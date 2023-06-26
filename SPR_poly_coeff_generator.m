% Polyfit on steps vs angles from .spr2 for each wavelength (2019 VT)
% Uses spectra from full angular scans

cal_save_folder = 'C:\Users\John\OneDrive - Chalmers University of Technology\Dahlin group\Python\PhD-Projects';
file_name = strcat(cal_save_folder, '\', 'SPR_poly_coeff_', datestr(now, 'yy-mm-dd'), '.csv');

data1_FC1 = dlmread('X-cal_210519scan1_L1 670nm.dto', '\t').';
data1_FC2 = dlmread('X-cal_210519scan1_L2 980nm.dto', '\t').';
data1_FC3 = dlmread('X-cal_210519scan1_L3 670nm.dto', '\t').';
data1_FC4 = dlmread('X-cal_210519scan1_L4 785nm.dto', '\t').';

y_FC1 = data1_FC1(1,:);
y_FC2 = data1_FC2(1,:);
y_FC3 = data1_FC3(1,:);
y_FC4 = data1_FC4(1,:);

step_max = 10 + 1*(length(data1_FC1)); % This number comes from number of measurement points in the calibration section of the .spr2 file
steps =  10:1:step_max-1;

p_FC1 = polyfit(steps, y_FC1, 20);
p_FC2 = polyfit(steps, y_FC2, 20);
p_FC3 = polyfit(steps, y_FC3, 20);
p_FC4 = polyfit(steps, y_FC4, 20);
polycoff_matrix = [p_FC1; p_FC2; p_FC3 ;p_FC4];

y_fit_FC1 = polyval(p_FC1, steps);
y_fit_FC2 = polyval(p_FC2, steps);
y_fit_FC3 = polyval(p_FC3, steps);
y_fit_FC4 = polyval(p_FC4, steps);

Chi_squared_FC1_20 = sum(((y_FC1-y_fit_FC1)./std(y_FC1)).^2);
Chi_squared_FC2_20 = sum(((y_FC2-y_fit_FC2)./std(y_FC2)).^2); 
Chi_squared_FC3_20 = sum(((y_FC3-y_fit_FC3)./std(y_FC3)).^2); 
Chi_squared_FC4_20 = sum(((y_FC4-y_fit_FC4)./std(y_FC4)).^2); 

% old_fit_FC1 = polyval(old_polys(1,:), steps);
% old_fit_FC2 = polyval(old_polys(2,:), steps);
% old_fit_FC3 = polyval(old_polys(3,:), steps);
% old_fit_FC4 = polyval(old_polys(4,:), steps);

writematrix(polycoff_matrix, file_name, 'Delimiter', 'tab')
