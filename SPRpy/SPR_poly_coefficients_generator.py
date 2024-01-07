import numpy as np
from datetime import datetime

# Run this script to get polynomial coefficients translating motor step values to angles (needed for spr2_to_csv.py)
# Cannot take angle values directly as they may differ slightly in number of points between instrument reboots (step motor homing sequence)
# Polyfit "over fitting" on steps vs angles from .spr2 for each wavelength.
# This gives absolute angle value accuracy within the noise limit of the homing sequence.
# Use spectra from full angular scans
# Check from time to time that the angles obtained from spr2_to_csv conversion matches the angles from .dto export for a given scan.
# If the step motor has been significantly affected or worn, it could be time to run this script again and replace this

cal_save_folder = ''
file_name = f"{cal_save_folder}/SPR_poly_coeff_{datetime.now().strftime('%y-%m-%d')}.csv"

data_L1 = np.loadtxt('X-cal_Full_range_L1 670nm.dto', delimiter='\t').T
data_L2 = np.loadtxt('X-cal_Full_range_L2 980nm.dto', delimiter='\t').T
data_L3 = np.loadtxt('X-cal_Full_range_L3 670nm.dto', delimiter='\t').T
data_L4 = np.loadtxt('X-cal_Full_range_L4 785nm.dto', delimiter='\t').T

y_L1 = data_L1[0, :]
y_L2 = data_L2[0, :]
y_L3 = data_L3[0, :]
y_L4 = data_L4[0, :]

scanspeed = 1  # Recommended to go for full range (slow) scans (1, 5 or 10 for slow, medium or fast scanspeed)
step_max = 10 + scanspeed * len(data_L1[0])  # This number comes from the number of measurement points in the calibration section of the .spr2 file
steps = np.arange(10, step_max)

p_L1 = np.polyfit(steps, y_L1, 20)
p_L2 = np.polyfit(steps, y_L2, 20)
p_L3 = np.polyfit(steps, y_L3, 20)
p_L4 = np.polyfit(steps, y_L4, 20)
polycoff_matrix = np.vstack((p_L1, p_L2, p_L3, p_L4))

y_fit_FC1 = np.polyval(p_L1, steps)
y_fit_FC2 = np.polyval(p_L2, steps)
y_fit_FC3 = np.polyval(p_L3, steps)
y_fit_FC4 = np.polyval(p_L4, steps)

Chi_squared_L1_c20 = np.sum(((y_L1 - y_fit_FC1) / np.std(y_L1)) ** 2)
Chi_squared_L2_c20 = np.sum(((y_L2 - y_fit_FC2) / np.std(y_L2)) ** 2)
Chi_squared_L3_c20 = np.sum(((y_L3 - y_fit_FC3) / np.std(y_L3)) ** 2)
Chi_squared_L4_c20 = np.sum(((y_L4 - y_fit_FC4) / np.std(y_L4)) ** 2)

print(Chi_squared_L1_c20)
print(Chi_squared_L2_c20)
print(Chi_squared_L3_c20)
print(Chi_squared_L4_c20)

np.savetxt(file_name, polycoff_matrix, delimiter='\t')

