# This file contains utility functions

import numpy as np
import tkinter
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory, asksaveasfilename
import pandas as pd
import re
from fresnel_transfer_matrix import TIR_determination


def select_folder(prompt, prompt_folder=None):
    root = tkinter.Tk()
    root.attributes("-topmost", 1)
    root.withdraw()
    selected_folder = askdirectory(title=prompt, parent=root, initialdir=prompt_folder)
    root.destroy()
    return selected_folder


def select_file(prompt, prompt_folder=None, file_types=[('Pickle files', '*.pickle')]):
    root = tkinter.Tk()
    root.attributes("-topmost", 1)
    root.withdraw()
    selected_file = askopenfilename(title=prompt, filetypes=file_types, initialdir=prompt_folder, parent=root)
    root.destroy()
    return selected_file


def select_files(prompt, prompt_folder=None, file_types=[('Pickle files', '*.pickle')]):
    root = tkinter.Tk()
    root.attributes("-topmost", 1)
    root.withdraw()
    selected_files = askopenfilenames(title=prompt, filetypes=file_types, initialdir=prompt_folder, parent=root)
    root.destroy()
    return selected_files


def save_file(prompt, prompt_folder=None, file_types=[('CSV files', '*.csv')], default_extension='.csv'):
    root = tkinter.Tk()
    root.attributes("-topmost", 1)
    root.withdraw()
    save_file = asksaveasfilename(title=prompt, filetypes=file_types, defaultextension=default_extension, initialdir=prompt_folder, parent=root)
    root.destroy()
    return save_file


def load_csv_data(path=False, default_data_folder=None, prompt='Select the measurement data file (.csv)'):
    if not path:
        print(prompt)
        data_path_ = select_file('Select the measurement data file (.csv)', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')])
    else:
        data_path_ = path

    #  Determine the scanning speed/step length if present in the file
    try:
        with open(data_path_, 'r') as file:
            step_length_pattern = re.compile(r'=\d{1,2}')
            scanspeed = int(step_length_pattern.search(file.readline()).group().strip('='))

    except AttributeError:  # I think .group().strip() should return AttributeError if .search() returns None
        scanspeed = 5  # Assuming medium scanspeed if legacy spr2 to csv converter was used


    # Load in the measurement data from a .csv file
    data_frame_ = pd.read_csv(data_path_, delimiter=';', skiprows=1, header=None)
    time_df = data_frame_.iloc[1:, 0]
    angles_df = data_frame_.iloc[0, 1:]
    ydata_df = data_frame_.iloc[1:, 1:]

    # Select last scan as default reflectivity plot
    reflectivity_df_ = pd.DataFrame(data={'angles': angles_df, 'ydata': ydata_df.iloc[-1, :]})

    return data_path_, scanspeed, time_df, angles_df, ydata_df, reflectivity_df_

# TODO: The new fitting parameters need to be supplied as pandas dataframes with each data point as a pandas series object
def calculate_sensorgram(time, angles, ydata, TIR_range, scanspeed, SPR_boundary_points=(70, 70), SPR_fit_points=4000, TIR_window_count=None, TIR_fit_lower_ind=None, TIR_fit_higher_ind=None):

    # Convert dataframes to numpy ndarrays
    time = time.to_numpy()
    angles = angles.to_numpy()
    ydata = ydata.to_numpy()

    # TODO: The new fitting parameters need to be supplied as pandas dataframes with each data point as a pandas series object
    # Calculating SPR and TIR angles
    sensorgram_SPR_angles = np.empty(len(ydata))
    sensorgram_SPR_fit_y = np.empty(len(ydata))
    sensorgram_SPR_fit_x = np.empty(len(ydata))
    sensorgram_SPR_angles.fill(np.nan)
    sensorgram_SPR_fit_x.fill(np.nan)
    sensorgram_SPR_fit_y.fill(np.nan)
    sensorgram_TIR_angles = np.empty(len(ydata))
    sensorgram_TIR_deriv_x = np.empty(len(ydata))
    sensorgram_TIR_deriv_y = np.empty(len(ydata))
    sensorgram_TIR_deriv_fit_x = np.empty(len(ydata))
    sensorgram_TIR_deriv_fit_y = np.empty(len(ydata))
    sensorgram_TIR_angles.fill(np.nan)
    sensorgram_TIR_deriv_x.fill(np.nan)
    sensorgram_TIR_deriv_y.fill(np.nan)
    sensorgram_TIR_deriv_fit_x.fill(np.nan)
    sensorgram_TIR_deriv_fit_y.fill(np.nan)

    for ind, val in enumerate(time):
        reflectivity_spectrum = ydata[ind-1, :]
        min_index = np.argmin(reflectivity_spectrum)

        # SPR angles
        try:
            y_selection = reflectivity_spectrum[min_index - SPR_boundary_points[0]:min_index + SPR_boundary_points[1]]

            polynomial = np.polyfit(angles[min_index - SPR_boundary_points[0]:min_index + SPR_boundary_points[1]],
                                    y_selection, 3)
            x_selection = np.linspace(angles[min_index - SPR_boundary_points[0]],
                                      angles[min_index + SPR_boundary_points[1]], SPR_fit_points)
            y_polyfit = np.polyval(polynomial, x_selection)
            y_fit_min_ind = np.argmin(y_polyfit)

            sensorgram_SPR_angles[ind-1] = x_selection[y_fit_min_ind]
            sensorgram_SPR_fit_x[ind-1] = x_selection
            sensorgram_SPR_fit_y[ind-1] = y_polyfit

        except:
            print('No SPR minimum found. Skipping measurement point...')
            sensorgram_SPR_angles[ind-1] = np.nan
            sensorgram_SPR_fit_x[ind-1] = np.nan
            sensorgram_SPR_fit_y[ind-1] = np.nan

        # TIR angles
        try:
            TIR_theta, TIR_xdata_filtered, deriv_ydata, TIR_theta_fit_x, TIR_theta_fit_y  = TIR_determination(angles, reflectivity_spectrum, TIR_range, scanspeed, _window_count=TIR_window_count, _fit_lower_ind=TIR_fit_lower_ind, _fit_higher_ind=TIR_fit_higher_ind)
            sensorgram_TIR_angles[ind-1] = TIR_theta
            sensorgram_TIR_deriv_x[ind-1] = TIR_xdata_filtered
            sensorgram_TIR_deriv_y[ind-1] = deriv_ydata
            sensorgram_TIR_deriv_fit_x[ind-1] = TIR_theta_fit_x
            sensorgram_TIR_deriv_fit_y[ind-1] = TIR_theta_fit_y

        except:
            print('No TIR found. Skipping measurement point...')
            sensorgram_TIR_angles[ind-1] = np.nan
            sensorgram_TIR_deriv_x[ind - 1] = np.nan
            sensorgram_TIR_deriv_y[ind - 1] = np.nan
            sensorgram_TIR_deriv_fit_x[ind-1] = np.nan
            sensorgram_TIR_deriv_fit_y[ind-1] = np.nan


    sensorgram_df = pd.DataFrame(data={'time': time, 'SPR angle': sensorgram_SPR_angles, 'TIR angle': sensorgram_TIR_angles, 'SPR fit x': sensorgram_SPR_fit_x, 'SPR fit y': sensorgram_SPR_fit_y, 'TIR deriv x': sensorgram_TIR_deriv_x, 'TIR deriv y': sensorgram_TIR_deriv_y, 'TIR deriv fit x': sensorgram_TIR_deriv_fit_x, 'TIR deriv fit y': sensorgram_TIR_deriv_fit_y})

    return sensorgram_df

