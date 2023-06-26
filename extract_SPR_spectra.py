import re
import os
import numpy as np
import tkinter
import multiprocessing as mp
from tkinter.filedialog import askopenfilename, askopenfilenames

def extract_parameters(content):

    #  Starting position for scans
    pos_pattern1 = re.compile(r'<ch number="\d" start_pos="\d+">')
    pos_line = pos_pattern1.findall(content)
    pos_pattern2 = re.compile(r'"\d+">')
    pos_value = pos_pattern2.search(pos_line[-1]).group()
    start_pos = int(pos_value.strip('">'))

    #  TIR steps for each laser in calibration
    TIR_pattern = re.compile(r'111;111;111;\d+;\d+;\d+;\d+;11111;11111;111;111')
    TIR_line = TIR_pattern.search(content).group()
    TIR_values = list(map(float, TIR_line.split(';')[3:7]))

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

    return start_pos, step_length, TIR_values, time_value_list


def extract_spectra(content, c_ind, polycoff, TIR_offset, start_pos, step_length, time_values, spr2_file):
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
            print('Row ', str(row), ': Mismatch in angle range compared to last measurement. This scan will be skipped.')
            continue

    # Generation of angles and combining with spectra
    spectra_steps = np.arange(start_pos, (step_length*points) + start_pos, step_length)
    spectra_angles = np.polyval(polycoff, spectra_steps) - TIR_offset
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
    calib_steps = np.arange(77.0, 27707.0, 10)  # Step length is always 10 for calibration
    calib_angles = np.polyval(polycoff, calib_steps) - TIR_offset
    calib_array = np.vstack((calib_angles, calib_data))

    #  Start intensity calibration
    for a_ind, angle in enumerate(spectra_full_array[0, :]):
        try:
            cal_ind = next(ind1 for ind1, cal_angle in enumerate(calib_array[0, :]) if np.isclose(cal_angle, angle, 0.001))  # Search for calibration value close to its respective spectra value at a given angle
            spectra_full_array[1:, a_ind] = np.true_divide(spectra_full_array[1:, a_ind], calib_array[1, cal_ind])   # Divide all values in  spectra_full_array[1:, a_ind] with calib_array[1, cal_ind]
        except:
            print('WARNING: Angle ', str(angle), ' was not calibrated')
            continue
    #  TODO: Change the calibration procedure to instead of looking for each value, find the starting and stopping steps in
    #   the calibration step array, then divide this into blocks where all angles between two blocks are calibrated
    #   against the same calibration values

    #  Add time values as first column
    time_values_np = np.array(time_values)
    spectra_full_array = np.column_stack((np.insert(time_values_np, 0, 0), spectra_full_array))

    #  Save data as .csv
    spr2_path, spr2_file_name = os.path.split(spr2_file)

    if c_ind == 0:
        file_identifier = spr2_file_name[:-5] + '-L1_670nm.csv'
    elif c_ind == 1:
        file_identifier = spr2_file_name[:-5] + '-L2_980nm.csv'
    elif c_ind == 2:
        file_identifier = spr2_file_name[:-5] + '-L3_670nm.csv'
    elif c_ind == 3:
        file_identifier = spr2_file_name[:-5] + '-L4_785nm.csv'

    save_name = os.path.join(spr2_path, file_identifier)
    header_string = 'Left most column is Time (min), First row is Angles (deg)'

    np.savetxt(save_name, spectra_full_array, fmt='%1.6f', delimiter=';', header=header_string)


if __name__ == '__main__':  # This is important since mp.Process goes through this file for extract_spectra()
    tkinter.Tk().withdraw()

    spr2_files = askopenfilenames(title='Select spr2 files')

    try:
        #  Read default polynomial file
        poly_file = r'C:\Users\anjohn\OneDrive - Chalmers\Dahlin group\Python\PhD-Projects\SPR_poly_coeff_22-01-29.csv'

        with open(poly_file, 'r') as p_file:
            polycoeffs = [0]*4
            for p_ind in range(4):
                coeff = p_file.readline().split('\t')
                polycoeffs[p_ind] = list(map(float, coeff))

        poly_path, poly_file_name = os.path.split(poly_file)

    except FileNotFoundError:
        print('Polynomial coefficients not found in default location')

        poly_file = askopenfilename(title='Select polynomial coefficients (.csv)')

        #  Read selected polynomial file
        with open(poly_file, 'r') as p_file:
            polycoeffs = [0]*4
            for p_ind in range(4):
                coeff = p_file.readline().split('\t')
                polycoeffs[p_ind] = list(map(float, coeff))

        poly_path, poly_file_name = os.path.split(poly_file)

    for spr2_file in spr2_files:

        spr2_path, spr2_file_name = os.path.split(spr2_file)

        #  Read sp2 file
        with open(spr2_file, 'r') as f:
            content = f.read()

        #  Get starting position for scan, TIR value steps and time points
        start_pos, scan_speed, TIR_steps, time_value_list = extract_parameters(content)

        #  Calculate angle for the TIR for each laser
        TIR_theoretical = [41.1479, 41.3866, 41.1479, 41.2802]
        angle_offsets = [0]*4
        for ind in range(4):
            angle_offsets[ind] = np.polyval(polycoeffs[ind], TIR_steps[ind]) - TIR_theoretical[ind]

        # extract_spectra(content, 0, polycoeffs[0], angle_offsets[0], start_pos, 10, time_value_list, sp2_file)

        #  Extract and calibrate spectra for each laser
       
        jobs = []
        for i in range(4):
            pr = mp.Process(target=extract_spectra, args=(content, i, polycoeffs[i], angle_offsets[i], start_pos, scan_speed, time_value_list, spr2_file))
            jobs.append(pr)
            pr.start()
            print('TIR Calib. Offset', i, ': ', angle_offsets[i], ' deg')
            print('Working...')

        # Wait for the first file to finish
        for job in jobs:
            job.join()

        print('File: ' + spr2_file_name + ' is done.')


