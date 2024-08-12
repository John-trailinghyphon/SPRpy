import re
import os
import numpy as np
import tkinter
import multiprocessing as mp
from tkinter.filedialog import askopenfilename, askopenfilenames

# Default .csv file for polynomial coefficents of fit for motor steps vs angle (see SPR_poly_coefficients_generator.py).
# This should be regenerated when used in a new instrument. Does generally not need to be regenerated after normal calibrations,
# unless the stepper motor calibration in the instrument has changed for some reason (e.g. may happen from homing errors or similar).
# Assumes all selected measurement files comes from the same instrument. Use multiple copies of spr2_to_csv_v2.py and
# their associated SPR_poly_coeff_YY-MM-DD.csv file if working with multiple instruments.

default_poly_file = r'DEFAULT_X_cal_values.csv'

prism_n_635 = 1.5215  # Approximate prism RI based on Soda lime glass
prism_n_650 = 1.5210  # Approximate prism RI based on Soda lime glass
prism_n_670 = 1.5202  # From Bionavis
prism_n_785 = 1.5162  # From Bionavis
prism_n_850 = 1.5149  # From Bionavis
prism_n_980 = 1.5129  # From Bionavis


def extract_parameters(content):

    # Calibration data
    calibration_device_pattern = re.compile(r'<calibration device_serial=.*>')
    cal_channels_pattern = re.compile(r'channels="\d+"')
    cal_step_pattern = re.compile(r'step_len="\d+"')
    cal_points_pattern = re.compile(r'points="\d+"')
    cal_start_pos_pattern = re.compile(r'start_pos="\d+"')
    digit_pattern = re.compile(r'\d+')

    calibration_param_line = calibration_device_pattern.findall(content)[0]
    channels = int(digit_pattern.findall(cal_channels_pattern.findall(calibration_param_line)[0])[0])
    cal_scanspeed = int(digit_pattern.findall(cal_step_pattern.findall(calibration_param_line)[0])[0])
    cal_points = int(digit_pattern.findall(cal_points_pattern.findall(calibration_param_line)[0])[0])
    cal_start_pos = int(digit_pattern.findall(cal_start_pos_pattern.findall(calibration_param_line)[0])[0])

    #  TIR steps for each laser in calibration
    TIR_pattern = re.compile(r'<prm_set set="P1500">.*</prm_set>')
    TIR_line = TIR_pattern.findall(content)[0]

    # set P1500 - Laser channels start at index pos 15 (starting from 0), ends at 22
    # TIR angle calibration step values starts at 26 and ends at 33
    laser_channels = list(map(int, TIR_line.split(';')[15+(8-channels):23]))

    #  Starting position for measurement scans
    pos_pattern1 = re.compile(r'<ch number="\d" start_pos="\d+">')
    pos_line = pos_pattern1.findall(content)
    pos_pattern2 = re.compile(r'"\d+">')
    pos_value = pos_pattern2.search(pos_line[-1]).group()
    start_pos = int(pos_value.strip('">'))

    #  List of measurement time for each point
    time_pattern = re.compile(r'<scan rtime="\d+"')
    time_matches = time_pattern.findall(content)
    time_value_list = [0]*len(time_matches)
    time_value_match = re.compile(r'\d+')

    for t_ind, time_string in enumerate(time_matches):
        time_value_list[t_ind] = int(time_value_match.search(time_string).group())/1000/60

    #  Determine the scanning speed/step length
    step_pattern = re.compile(r'<scan rtime="\d+" step_len="\d{1,2}" dir="Forward">')
    step_string = step_pattern.search(content).group()
    step_length_pattern = re.compile(r'"\d{1,2}"')
    step_length = int(step_length_pattern.search(step_string).group().strip('"'))

    print('Start position: ', start_pos)
    print('Scan speed: ', step_length)

    return start_pos, step_length, channels, laser_channels, time_value_list, cal_scanspeed, cal_points, cal_start_pos


def extract_spectra(content, c_ind, polycoff, start_pos, scanspeed, cal_scanspeed, cal_points, cal_start_pos, time_values, laser_channels, spr2_file):
    #  Extracts and calibrates spectra from .sp2 file, then saves it as .csv

    #  Get the spectra data (angles and intensity)
    spectra_pattern = re.compile(r'<ch number="' + str(c_ind) + r'" start_pos="\d+">.*</ch>')
    channel = spectra_pattern.findall(content)

    # Remove init_scan if it is not the only measurement
    init_scan_pattern = re.compile(r'<init_scan rtime="[1-9]')
    init_scan_match = re.search(init_scan_pattern, content)
    if init_scan_match and len(channel) > 1:
        channel.pop(0)

    data_pattern = re.compile(r'>.*<')
    point_string = data_pattern.search(channel[-1]).group().strip('><')
    points = len(point_string.split(';'))
    spectra_array = np.ones((len(channel), points))

    for row, match in enumerate(channel):
        spectra = data_pattern.search(match).group().strip('><')
        try:
            spectra_array[row, :] = list(map(float, spectra.split(';')))
        except ValueError:
            print('Row ', str(row), ': Mismatch in angle range compared to last scan. This scan will be skipped.')
            continue

    # Generation of angles and combining with spectra
    spectra_steps = np.arange(start_pos, (scanspeed * points) + start_pos, scanspeed)
    spectra_angles = np.polyval(polycoff, spectra_steps)
    spectra_full_array = np.vstack((spectra_angles, spectra_array))

    #  Get the calibration data
    calib_channel_pattern = re.compile(r'<number>' + str(c_ind) + r'.*</data>', flags=re.DOTALL)
    calib_channel_string = calib_channel_pattern.search(content).group()

    calib_data_pattern = re.compile(r'a>.*<')
    calib_data_string = calib_data_pattern.search(calib_channel_string).group()
    
    calib_data_string = calib_data_string.replace('a>', '0.')
    calib_data_string = calib_data_string.strip('<')
    calib_data_string = calib_data_string.replace(';', ';0.')
    
    calib_data = list(map(float, calib_data_string.split(';')))
    calib_steps = np.arange(float(cal_start_pos), float(cal_scanspeed*cal_points)+float(cal_start_pos), cal_scanspeed)
    calib_angles = np.polyval(polycoff, calib_steps)
    calib_array = np.vstack((calib_angles, calib_data))

    #  Start intensity calibration
    for a_ind, angle in enumerate(spectra_full_array[0, :]):
        try:
            cal_ind = next(ind1 for ind1, cal_angle in enumerate(calib_array[0, :]) if np.isclose(cal_angle, angle, 0.001))  # Search for calibration value close to its respective spectra value at a given angle
            spectra_full_array[1:, a_ind] = np.true_divide(spectra_full_array[1:, a_ind], calib_array[1, cal_ind])   # Divide all values in  spectra_full_array[1:, a_ind] with calib_array[1, cal_ind]
        except:
            try:
                # Try wider tolerance
                cal_ind = next(ind1 for ind1, cal_angle in enumerate(calib_array[0, :]) if np.isclose(cal_angle, angle, 0.003))  # Search for calibration value close to its respective spectra value at a given angle
                spectra_full_array[1:, a_ind] = np.true_divide(spectra_full_array[1:, a_ind], calib_array[1, cal_ind])  # Divide all values in  spectra_full_array[1:, a_ind] with calib_array[1, cal_ind]
            except:
                print('WARNING: Angle ', str(angle), ' was not calibrated')
                continue

    #  Add time values as first column
    time_values_np = np.array(time_values)
    spectra_full_array = np.column_stack((np.insert(time_values_np, 0, 0), spectra_full_array))

    #  Save data as .csv
    spr2_path, spr2_file_name = os.path.split(spr2_file)
    file_identifier = '{head}-L{index}_{wavelength}nm.csv'.format(head=spr2_file_name[:-5], index=str(c_ind+1), wavelength=str(laser_channels[c_ind]))
    save_name = os.path.join(spr2_path, file_identifier)
    header_string = 'Left most column is Time (min), First row is Angles (deg), Scanspeed=' + str(scanspeed)

    np.savetxt(save_name, spectra_full_array, fmt='%1.6f', delimiter=';', header=header_string)


if __name__ == '__main__':  # This is important since mp.Process goes through this file for extract_spectra()
    tkinter.Tk().withdraw()

    spr2_files = askopenfilenames(title='Select spr2 files')

    for file_ind, spr2_file in enumerate(spr2_files):

        spr2_path, spr2_file_name = os.path.split(spr2_file)

        #  Read sp2 file
        with open(spr2_file, 'r') as f:
            content = f.read()

        #  Get various parameters from file
        start_pos, scan_speed, channels, laser_channels, time_value_list, cal_scanspeed, cal_points, cal_start_pos = extract_parameters(content)

        # Assume all selected files are from the same instrument and has the same polynomial coefficients
        if file_ind == 0:
            try:
                #  Read default polynomial file
                with open(default_poly_file, 'r') as p_file:
                    polycoeffs = [0] * channels
                    for p_ind in range(channels):
                        coeff = p_file.readline().split('\t')
                        polycoeffs[p_ind] = list(map(float, coeff))

                poly_path, poly_file_name = os.path.split(default_poly_file)

            except FileNotFoundError:
                print('Polynomial coefficients not found in default location')

                poly_file = askopenfilename(title='Select polynomial coefficients (.csv)')

                #  Read selected polynomial file
                with open(poly_file, 'r') as p_file:
                    polycoeffs = [0] * channels
                    for p_ind in range(channels):
                        coeff = p_file.readline().split('\t')
                        polycoeffs[p_ind] = list(map(float, coeff))

                poly_path, poly_file_name = os.path.split(poly_file)

        #  Extract and calibrate spectra for each laser
        jobs = []
        for i in range(channels):
            pr = mp.Process(target=extract_spectra, args=(content, i, polycoeffs[i], start_pos, scan_speed, cal_scanspeed, cal_points, cal_start_pos, time_value_list, laser_channels, spr2_file))
            jobs.append(pr)
            pr.start()
        print('Working...')

        # Wait for the first file to finish
        for job in jobs:
            job.join()

        print('File: ' + spr2_file_name + ' is done.')


