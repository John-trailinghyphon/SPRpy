import os
import numpy as np
from datetime import datetime
import SPRpy_functions
import SPRpy_spr2_to_csv
# Run this script to get polynomial coefficients translating motor step values to angles (needed for spr2_to_csv.py)
# Cannot take angle values directly as they may differ slightly in number of points between instrument reboots (step motor homing sequence)
# Polyfit "over fitting" on steps vs angles from .spr2 for each wavelength.
# This gives absolute angle value accuracy within the noise limit of the homing sequence.
# Use spectra from full angular scans
# Check from time to time that the angles obtained from spr2_to_csv conversion matches the angles from .dto export for a given scan.
# Small deviations of < 10-20 mdeg in absolute angle is likely within normal instrument uncertainty.
# If the step motor has been significantly affected or worn, it could be time to run this script again and replace this file for the spr2_to_csv.py script

# Prompt for required files (they are sorted by default)
spr2_file = SPRpy_functions.select_file('Select .spr2 file of full range slow speed angular scan.', file_types=[('SPR2 file', '*.spr2')])
dto_files = SPRpy_functions.select_files('Select .dto Bionavs viewer exported files for all wavelengths of full range slow speed angular scan.', file_types=[('DTO file', '*.dto')])
dto_wavelengths = [int(dto_file[-9:-6]) for dto_file in dto_files]

# Read .spr2 file and get calibration values
with open(spr2_file, 'r') as file:
    spr2_content = file.read()
start_pos, scanspeed, _, _, _, _, _, _ = SPRpy_spr2_to_csv.extract_parameters(spr2_content)

# Get polynomial coefficients
polycoffs = []
for dto_file in dto_files:
    data = np.loadtxt(dto_file, delimiter='\t', converters=lambda s: s.replace(',', '.')).T
    angles = data[0, :]
    step_max = start_pos + scanspeed * len(angles)
    steps = np.arange(start_pos, step_max, scanspeed)
    p = np.polyfit(steps, angles, 6)
    polycoffs.append(p)
    y_fit = np.polyval(p, steps)
    Chi_squared = np.sum(((angles - y_fit) / np.std(angles)) ** 2)
    print(Chi_squared)

# Write result to .csv file
polycoff_matrix = np.vstack(tuple(polycoffs))
cal_save_folder = os.getcwd()
file_name = f"{cal_save_folder}/SPRpy_X_cal_values_{datetime.now().strftime('%y-%m-%d')}.csv"
np.savetxt(file_name, polycoff_matrix, delimiter='\t')

