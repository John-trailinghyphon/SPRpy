import numpy as np
import tkinter
from tkinter.filedialog import askopenfilename, askdirectory
import pandas as pd
import re
from fresnel_transfer_matrix import TIR_determination


# def generate_id():
#     """
#     Each time this function is called within an object it can return a new ID.
#     :yield:
#     """
#     new_id = 0
#     while True:
#         yield new_id
#         new_id += 1


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


def load_csv_data(path=False):
    if not path:
        print('Select the measurement data file (.csv)')
        data_path_ = select_file(prompt='Select the measurement data file', file_types=[('CSV files', '*.csv')])
    else:
        data_path_ = path

    #  Determine the scanning speed/step length if present in the file
    try:
        with open(data_path_, 'r') as file:
            step_length_pattern = re.compile(r'=\d{1,2}')
            scanspeed = int(step_length_pattern.search(file.readline()).group().strip('='))

    except AttributeError:  # I think .group().strip() should return attribute error if .search() returns None
        scanspeed = 5  # Assuming medium

    # Load in the measurement data from a .csv file
    data_frame_ = pd.read_csv(data_path_, delimiter=';', skiprows=1, header=None)
    time_df = data_frame_.iloc[1:, 0]
    angles_df = data_frame_.iloc[0, 1:]
    ydata_df = data_frame_.iloc[1:, 1:]

    # Select last scan as default reflectivity plot
    reflectivity_df_ = pd.DataFrame(data={'angles': angles_df, 'ydata': ydata_df.iloc[-1, :]})

    return data_path_, scanspeed, time_df, angles_df, ydata_df, reflectivity_df_


def calculate_sensorgram(time, angles, ydata, TIR_range, scanspeed, SPR_points=(70, 70)):

    # Convert dataframes to numpy ndarrays
    time = time.to_numpy()
    angles = angles.to_numpy()
    ydata = ydata.to_numpy()

    # Calculating SPR and TIR angles
    sensorgram_SPR_angles = np.empty(len(ydata)) * np.nan
    sensorgram_TIR_angles = np.empty(len(ydata)) * np.nan
    for ind, val in enumerate(time):
        reflectivity_spectrum = ydata[ind-1, :]
        min_index = np.argmin(reflectivity_spectrum)

        # SPR angles
        try:
            y_selection = reflectivity_spectrum[min_index-SPR_points[0]:min_index+SPR_points[1]]

            polynomial = np.polyfit(angles[min_index - SPR_points[0]:min_index + SPR_points[1]],
                                    y_selection, 3)
            x_selection = np.linspace(angles[min_index - SPR_points[0]],
                                      angles[min_index + SPR_points[1]], 4000)
            y_polyfit = np.polyval(polynomial, x_selection)
            y_fit_min_ind = np.argmin(y_polyfit)

            sensorgram_SPR_angles[ind-1] = x_selection[y_fit_min_ind]

        except:
            print('No SPR minimum found. Skipping measurement point...')
            sensorgram_SPR_angles[ind-1] = np.NaN

        # TIR angles
        try:
            TIR_theta, _, _ = TIR_determination(angles, reflectivity_spectrum, TIR_range, scanspeed)
            sensorgram_TIR_angles[ind-1] = TIR_theta

        except:
            print('No TIR found. Skipping measurement point...')
            sensorgram_TIR_angles[ind-1] = np.NaN

    sensorgram_df = pd.DataFrame(data={'time': time, 'SPR angle': sensorgram_SPR_angles, 'TIR angle': sensorgram_TIR_angles})

    return sensorgram_df

