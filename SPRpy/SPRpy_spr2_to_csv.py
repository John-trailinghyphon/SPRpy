import re
import os
import sys
import numpy as np
import tkinter
import tomllib
from tkinter.filedialog import askopenfilename, askopenfilenames

# Faster version of SPRpy_spr2_to_csv.py: same config.toml, same X-cal file, same .csv file names and format.
#
# Where the old version spent its time (Wet Scan.spr2, 170 MB, 4 channels: 89 s, this version: 7 s):
#  - ~98 %: the intensity calibration loop. For every angle it searched the calibration table from the start,
#    calling np.isclose() on two scalars each time (~6 million calls per channel for a wet scan).
#    Here the same lookup is done for all angles at once with array operations.
#  - The whole file was pickled into every worker process (170 MB each) and searched again there. Here the file
#    is parsed once, in one process. With the calibration vectorised there is nothing heavy left to parallelise.
#
# Calibration lookup. The old version used, for each angle, the FIRST calibration point within a *relative*
# tolerance of 0.1 % (then 0.3 %) of the angle, i.e. 0.04-0.08 deg. The calibration points are ~0.014 deg apart,
# so this picks a point 3-5 steps BELOW the angle instead of the nearest one (up to ~0.4 % intensity error).
#   'interp':   linear interpolation of the calibration curve at each scan angle (default)
#   'original': same lookup as the old version, gives .csv files identical to it
CALIBRATION_MODE = 'interp'

# Read configuration parameters
with open('config.toml', 'r') as f:
    config = tomllib.loads(f.read())

default_poly_file = config["default_poly_file"]


def read_polycoeffs(poly_file, channels):
    with open(poly_file, 'r') as p_file:
        return [list(map(float, p_file.readline().split('\t'))) for _ in range(channels)]


def load_polycoeffs(channels):
    try:
        return read_polycoeffs(default_poly_file, channels)
    except FileNotFoundError:
        print('Error: Polynomial coefficients not found in default location')
    except ValueError:
        print('Error: Default polynomial coefficients not matching number of channels')

    poly_file = askopenfilename(title='Error! Select correct polynomial coefficients (.csv)')
    if not poly_file:
        raise FileNotFoundError('No X-calibration file (polynomial coefficients) selected, conversion cancelled.')
    return read_polycoeffs(poly_file, channels)


def extract_parameters(content):

    # The calibration is in the first few kB of the file, so only search there for it
    calibration_end = content.find('</calibration>')
    head = content[:calibration_end] if calibration_end != -1 else content

    # Calibration data
    calibration_param_line = re.search(r'<calibration device_serial=.*>', head).group()
    channels = int(re.search(r'channels="(\d+)"', calibration_param_line).group(1))
    cal_scanspeed = int(re.search(r'step_len="(\d+)"', calibration_param_line).group(1))
    cal_points = int(re.search(r'points="(\d+)"', calibration_param_line).group(1))
    cal_start_pos = int(re.search(r'start_pos="(\d+)"', calibration_param_line).group(1))

    # set P1500 - Laser channels start at index pos 15 (starting from 0), ends at 22
    TIR_line = re.search(r'<prm_set set="P1500">.*</prm_set>', head).group()
    laser_channels = list(map(int, TIR_line.split(';')[15+(8-channels):23]))

    # Calibration intensity of each channel (the <data> line following <number>N)
    calib_data = []
    for c_ind in range(channels):
        calib_channel_string = head[head.index('<number>' + str(c_ind)):]
        calib_data_string = re.search(r'a>.*<', calib_channel_string).group().strip('a>').strip('<')
        calib_data.append(np.array(list(map(float, calib_data_string.split(';'))))/10000)

    # Spectra of all channels in one pass over the file: <ch number="N" start_pos="S">values</ch>
    spectra_strings = [[] for _ in range(channels)]
    start_pos = None
    for match in re.finditer(r'<ch number="(\d)" start_pos="(\d+)">(.*)</ch>', content):
        c_ind = int(match.group(1))
        if c_ind < channels:
            spectra_strings[c_ind].append(match.group(3))
        start_pos = int(match.group(2))

    # Remove init_scan if it is not the only measurement
    if re.search(r'<init_scan rtime="[1-9]', content):
        for channel in spectra_strings:
            if len(channel) > 1:
                channel.pop(0)

    #  List of measurement time for each point
    time_value_list = [int(rtime)/1000/60 for rtime in re.findall(r'<scan rtime="(\d+)"', content)]

    #  Determine the scanning speed/step length
    step_length = int(re.search(r'<scan rtime="\d+" step_len="(\d{1,2})" dir="Forward">', content).group(1))

    print('Start position: ', start_pos)
    print('Scan speed: ', step_length)

    return start_pos, step_length, channels, laser_channels, time_value_list, cal_scanspeed, cal_points, cal_start_pos, calib_data, spectra_strings


def first_close_index(calib_angles, angles, rtol):
    # Vectorised version of the old search, for all angles at once:
    #     next(ind1 for ind1, cal_angle in enumerate(calib_angles) if np.isclose(cal_angle, angle, rtol))
    # np.isclose(a, b, rtol) is abs(a - b) <= 1e-8 + rtol*abs(b). Gives -1 where no calibration angle is close enough.
    cal_ind = np.full(len(angles), -1)
    for chunk in range(0, len(angles), 500):  # 500 angles at a time to limit memory use
        angle = angles[chunk:chunk+500, np.newaxis]
        is_close = np.abs(calib_angles - angle) <= 1e-8 + rtol*np.abs(angle)
        found = is_close.any(axis=1)
        cal_ind[chunk:chunk+500][found] = is_close.argmax(axis=1)[found]
    return cal_ind


def extract_spectra(spectra_strings, c_ind, polycoff, start_pos, scanspeed, cal_scanspeed, cal_points, cal_start_pos, calib_data, time_values, laser_channels, spr2_file):
    #  Converts and calibrates the spectra of one channel, then saves it as .csv

    points = len(spectra_strings[-1].split(';'))
    spectra_array = np.ones((len(spectra_strings), points))

    for row, spectra in enumerate(spectra_strings):
        try:
            spectra_array[row, :] = list(map(float, spectra.split(';')))
        except ValueError:
            print('Row ', str(row), ': Mismatch in angular resolution compared to last scan. This scan will be skipped.')
            continue

    # Generation of angles for spectra and calibration
    spectra_steps = np.arange(start_pos, (scanspeed * points) + start_pos, scanspeed)
    spectra_angles = np.polyval(polycoff, spectra_steps)

    calib_steps = np.arange(float(cal_start_pos), float(cal_scanspeed*cal_points)+float(cal_start_pos), cal_scanspeed)
    calib_angles = np.polyval(polycoff, calib_steps)

    #  Intensity calibration of all angles at once
    if CALIBRATION_MODE == 'interp':
        if np.any(np.diff(calib_angles) <= 0):
            raise ValueError('Calibration angles are not increasing, cannot interpolate')
        spectra_array /= np.interp(spectra_angles, calib_angles, calib_data)
    else:
        cal_ind = first_close_index(calib_angles, spectra_angles, 0.001)
        not_found = cal_ind == -1
        cal_ind[not_found] = first_close_index(calib_angles, spectra_angles[not_found], 0.003)  # Try wider tolerance
        for angle in spectra_angles[cal_ind == -1]:
            print('WARNING: Angle ', str(angle), ' was not calibrated')
        found = cal_ind != -1
        spectra_array[:, found] /= calib_data[cal_ind[found]]

    #  Add angles as first row and time values as first column
    spectra_full_array = np.vstack((spectra_angles, spectra_array))
    spectra_full_array = np.column_stack((np.insert(np.array(time_values), 0, 0), spectra_full_array))

    #  Save data as .csv
    spr2_path, spr2_file_name = os.path.split(spr2_file)
    save_name = csv_file_name(spr2_file, c_ind, laser_channels[c_ind])
    header_string = 'Left most column is Time (min), First row is Angles (deg), Scanspeed=' + str(scanspeed)

    np.savetxt(save_name, spectra_full_array, fmt='%1.6f', delimiter=';', header=header_string)


def csv_file_name(spr2_file, channel_index, wavelength):
    """
    Name of the .csv file of one laser channel of a converted .spr2 file (same folder as the .spr2 file).
    """
    spr2_path, spr2_file_name = os.path.split(spr2_file)
    return os.path.join(spr2_path, '{head}-L{index}_{wavelength}nm.csv'.format(head=spr2_file_name[:-5], index=str(channel_index + 1), wavelength=str(wavelength)))


def convert_spr2_file(spr2_file, polycoeffs=None, progress=None):
    """
    Converts and calibrates one .spr2 file to one .csv file per laser channel (saved next to the .spr2 file).

    :param polycoeffs: X-calibration polynomial coefficients of each channel, read from the default file if None
    :param progress: optional function that is called with a progress message (e.g. to show it in the start menu)
    :return: the polynomial coefficients used (to reuse them for further files from the same instrument)
    """
    spr2_path, spr2_file_name = os.path.split(spr2_file)
    if progress:
        progress('Reading ' + spr2_file_name + '...')

    #  Read sp2 file
    with open(spr2_file, 'r') as f:
        content = f.read()

    #  Get various parameters and the raw spectra from file
    start_pos, scan_speed, channels, laser_channels, time_value_list, cal_scanspeed, cal_points, cal_start_pos, calib_data, spectra_strings = extract_parameters(content)
    del content

    # Assume all selected files are from the same instrument and has the same polynomial coefficients
    if polycoeffs is None:
        polycoeffs = load_polycoeffs(channels)

    #  Calibrate and save the spectra of each laser
    for c_ind in range(channels):
        print('Working on channel L' + str(c_ind + 1) + ' ' + str(laser_channels[c_ind]) + 'nm')
        if progress:
            progress('Converting {file}: channel L{index} ({wavelength} nm)...'.format(file=spr2_file_name, index=c_ind + 1, wavelength=laser_channels[c_ind]))
        extract_spectra(spectra_strings[c_ind], c_ind, polycoeffs[c_ind], start_pos, scan_speed, cal_scanspeed, cal_points, cal_start_pos, calib_data[c_ind], time_value_list, laser_channels, spr2_file)

    print('File: ' + spr2_file_name + ' is done.')
    return polycoeffs


if __name__ == '__main__':
    tkinter.Tk().withdraw()

    spr2_files = sys.argv[1:]  # Files can also be given on the command line
    if not spr2_files:
        spr2_files = askopenfilenames(title='Select spr2 files')

    polycoeffs = None
    for spr2_file in spr2_files:
        polycoeffs = convert_spr2_file(spr2_file, polycoeffs)
