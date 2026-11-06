# This file contains utility functions

import copy
import math
import tomllib
import numpy as np
import tkinter
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory, asksaveasfilename
import pandas as pd
import re
from fresnel_transfer_matrix import TIR_determination, SPR_determination, select_TIR_range, angle_range_around_minimum, SPR_minimum_index, TIR_METHODS, SPR_METHODS


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


def load_config(config_path='config.toml', defaults_path='config.toml.defaults'):
    """
    Reads config.toml and fills in any entries that are missing (e.g. after upgrading SPRpy while keeping an older
    config.toml) from config.toml.defaults. Missing entries are reported, so that they can be added to config.toml.

    :return: config dictionary
    """
    with open(config_path, 'rb') as f:
        config = tomllib.load(f)
    try:
        with open(defaults_path, 'rb') as f:
            defaults = tomllib.load(f)
    except FileNotFoundError:
        return config

    missing = []

    def merge(target, source, prefix=''):
        for key, value in source.items():
            if key not in target:
                target[key] = copy.deepcopy(value)
                missing.append(prefix + key)
            elif isinstance(value, dict) and isinstance(target[key], dict):
                merge(target[key], value, prefix + key + '.')

    merge(config, defaults)
    if missing:
        print('Note: the following settings are missing in config.toml, using the default values from '
              'config.toml.defaults instead: ' + ', '.join(missing))
    return config


class FileSelectionCancelled(Exception):
    """
    Raised when the user cancels a file dialog.
    """


def load_csv_data(path=False, default_data_folder=None, prompt='Select the measurement data file (.csv)'):
    if not path:
        print(prompt)
        data_path_ = select_file(prompt, prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')])
        if not data_path_:
            raise FileSelectionCancelled('No measurement file selected.')
    else:
        data_path_ = path

    #  Determine the scanning speed/step length if present in the file
    try:
        with open(data_path_, 'r') as file:
            step_length_pattern = re.compile(r'=\d{1,2}')
            scanspeed = int(step_length_pattern.search(file.readline()).group().strip('='))

    except AttributeError:  # .search() returns None if no scanspeed is present in the file
        scanspeed = 5  # Assuming medium scanspeed if legacy spr2 to csv converter was used

    # Load in the measurement data from a .csv file
    data_frame_ = pd.read_csv(data_path_, delimiter=';', skiprows=1, header=None)

    # Remove empty columns and rows (e.g. from trailing delimiters or line breaks added by other programs like Excel)
    data_frame_ = data_frame_.dropna(axis=1, how='all').dropna(axis=0, how='all').reset_index(drop=True)
    data_frame_.columns = range(data_frame_.shape[1])

    time_df = data_frame_.iloc[1:, 0]
    angles_df = data_frame_.iloc[0, 1:]
    ydata_df = data_frame_.iloc[1:, 1:]

    # Select last scan as default reflectivity plot
    reflectivity_df_ = pd.DataFrame(data={'angles': angles_df, 'ydata': ydata_df.iloc[-1, :]})

    return data_path_, scanspeed, time_df, angles_df, ydata_df, reflectivity_df_


def default_SPR_TIR_fitting_parameters(TIR_default_parameters, SPR_default_parameters):
    """
    Builds the SPR and TIR fitting parameters of a new session from the [TIR_fitting_parameters] and
    [SPR_fitting_parameters] sections of config.toml.
    """
    parameters = {
        'TIR range': [float(i) for i in TIR_default_parameters['TIR_range_water_or_long_measurement']],
        'Fresnel_angle_range_points': [int(i) for i in SPR_default_parameters['Fresnel_angle_range_points']],
        'sensorgram_angle_range_points': [int(i) for i in SPR_default_parameters['sensorgram_angle_range_points']],
        'window_count_scanspeeds_1_5': int(TIR_default_parameters['window_count_scanspeeds_1_5']),
        'points_above_TIR_peak_scanspeed_1_5': int(TIR_default_parameters['points_above_TIR_peak_scanspeed_1_5']),
        'points_below_TIR_peak_scanspeed_1_5': int(TIR_default_parameters['points_below_TIR_peak_scanspeed_1_5']),
        'window_count_scanspeeds_10': int(TIR_default_parameters['window_count_scanspeeds_10']),
        'points_above_TIR_peak_scanspeed_10': int(TIR_default_parameters['points_above_TIR_peak_scanspeed_10']),
        'points_below_TIR_peak_scanspeed_10': int(TIR_default_parameters['points_below_TIR_peak_scanspeed_10']),
        'TIR fit points': int(TIR_default_parameters['TIR_fit_points']),
        'SPR fit points': int(SPR_default_parameters['SPR_fit_points']),
        'TIR method': str(TIR_default_parameters.get('TIR_method', 'savitzky_golay')),
        'TIR smooth width deg': float(TIR_default_parameters.get('TIR_smooth_width_deg', 0.15)),
        'TIR fit halfwidth deg': float(TIR_default_parameters.get('TIR_fit_halfwidth_deg', 0.045)),
        'SPR method': str(SPR_default_parameters.get('SPR_method', 'dip_level')),
        'SPR dip level': float(SPR_default_parameters.get('SPR_dip_level', 0.2)),
    }
    if parameters['TIR method'] not in TIR_METHODS:
        print('Warning! Unknown TIR_method "{m}" in config.toml, using "savitzky_golay".'.format(m=parameters['TIR method']))
        parameters['TIR method'] = 'savitzky_golay'
    if parameters['SPR method'] not in SPR_METHODS:
        print('Warning! Unknown SPR_method "{m}" in config.toml, using "dip_level".'.format(m=parameters['SPR method']))
        parameters['SPR method'] = 'dip_level'
    return parameters


def ensure_SPR_TIR_fitting_parameters(parameters, TIR_default_parameters):
    """
    Adds fitting parameters that are missing in sessions created with older SPRpy versions. Old sessions keep the
    legacy TIR and SPR algorithms by default so that previous results are reproduced (they can be switched in the GUI).
    """
    parameters.setdefault('TIR method', 'legacy')
    parameters.setdefault('SPR method', 'legacy')
    parameters.setdefault('TIR smooth width deg', float(TIR_default_parameters.get('TIR_smooth_width_deg', 0.15)))
    parameters.setdefault('TIR fit halfwidth deg', float(TIR_default_parameters.get('TIR_fit_halfwidth_deg', 0.045)))
    parameters.setdefault('SPR dip level', 0.2)
    return parameters


def update_TIR_settings_for_measurement(parameters, scanspeed, angles_df, ydata_df, TIR_default_parameters):
    """
    Chooses the TIR range (liquid or air) and the scan speed dependent legacy TIR settings for a newly loaded
    measurement. The TIR range is chosen from where the TIR edge actually is in the last spectrum (the range with the
    steepest increase in reflectivity). If that fails, more than 50 scans is assumed to be a liquid measurement.
    """
    liquid_range = [float(i) for i in TIR_default_parameters['TIR_range_water_or_long_measurement']]
    air_range = [float(i) for i in TIR_default_parameters['TIR_range_air_or_few_scans']]
    try:
        selected_range = select_TIR_range(angles_df, ydata_df.iloc[-1, :], [liquid_range, air_range])
    except Exception:
        selected_range = None
    if selected_range is None:
        selected_range = liquid_range if ydata_df.shape[0] > 50 else air_range
    parameters['TIR range'] = selected_range

    if scanspeed <= 5:
        parameters['TIR window count'] = parameters['window_count_scanspeeds_1_5']
        parameters['points_above_TIR_peak'] = parameters['points_above_TIR_peak_scanspeed_1_5']
        parameters['points_below_TIR_peak'] = parameters['points_below_TIR_peak_scanspeed_1_5']
    else:
        parameters['TIR window count'] = parameters['window_count_scanspeeds_10']
        parameters['points_above_TIR_peak'] = parameters['points_above_TIR_peak_scanspeed_10']
        parameters['points_below_TIR_peak'] = parameters['points_below_TIR_peak_scanspeed_10']
    return parameters


def calculate_sensorgram(time, angles, ydata, SPR_TIR_fitting_parameters):
    """
    Calculates the SPR and TIR angle of every spectrum of a measurement.

    Spectra where an angle cannot be determined get the value 0 (and a message is printed), instead of NaN, so that
    they stay accessible when hovering over the sensorgram (e.g. for Pd or Pt films).

    :return: DataFrame with columns 'time', 'SPR angle' and 'TIR angle'
    """

    # Convert dataframes to numpy ndarrays
    time = np.asarray(time, dtype=np.float64)
    angles = np.asarray(angles, dtype=np.float64)
    ydata = np.asarray(ydata, dtype=np.float64)

    sensorgram_SPR_angles = np.zeros(len(ydata))
    sensorgram_TIR_angles = np.zeros(len(ydata))
    failed_SPR = []
    failed_TIR = []

    for ind in range(len(ydata)):
        reflectivity_spectrum = ydata[ind, :]

        try:
            sensorgram_SPR_angles[ind] = SPR_determination(angles, reflectivity_spectrum, SPR_TIR_fitting_parameters)[0]
        except Exception:
            failed_SPR.append(time[ind])

        try:
            sensorgram_TIR_angles[ind] = TIR_determination(angles, reflectivity_spectrum, SPR_TIR_fitting_parameters)[0]
        except Exception:
            failed_TIR.append(time[ind])

    for name, failed in (('SPR', failed_SPR), ('TIR', failed_TIR)):
        if failed:
            print('No {name} angle found for {n} of {total} spectra (set to 0), first at t = {t:.2f} min. '
                  'Check the {name} fitting parameters.'.format(name=name, n=len(failed), total=len(ydata), t=failed[0]))

    sensorgram_df = pd.DataFrame(data={'time': time, 'SPR angle': sensorgram_SPR_angles, 'TIR angle': sensorgram_TIR_angles})

    return sensorgram_df


def sensorgram_fit_details(angles, reflectivity_spectrum, SPR_TIR_fitting_parameters):
    """
    Details of the TIR and SPR angle determination of one spectrum, used for plotting the fits in the GUI. NaN-filled
    placeholders are returned for a fit that fails.
    """
    angles = np.asarray(angles, dtype=np.float64)
    reflectivity_spectrum = np.asarray(reflectivity_spectrum, dtype=np.float64)
    details = {}
    try:
        TIR_theta, deriv_x, deriv_y, fit_x, fit_y = TIR_determination(angles, reflectivity_spectrum, SPR_TIR_fitting_parameters)
        details.update({'TIR angle': TIR_theta, 'TIR deriv x': deriv_x, 'TIR deriv y': deriv_y, 'TIR deriv fit x': fit_x, 'TIR deriv fit y': fit_y})
    except Exception:
        details.update({'TIR angle': np.nan, 'TIR deriv x': [np.nan], 'TIR deriv y': [np.nan], 'TIR deriv fit x': [np.nan], 'TIR deriv fit y': [np.nan]})
    try:
        SPR_theta, fit_x, fit_y, window_x, window_y = SPR_determination(angles, reflectivity_spectrum, SPR_TIR_fitting_parameters)

        # Show the data a bit outside the fit window as well
        margin = max(10, len(window_x) // 2)
        start = max(int(np.searchsorted(angles, window_x[0])) - margin, 0)
        stop = min(int(np.searchsorted(angles, window_x[-1])) + margin + 1, len(angles))
        details.update({'SPR angle': SPR_theta, 'SPR fit x': fit_x, 'SPR fit y': fit_y,
                        'SPR data x': angles[start:stop], 'SPR data y': reflectivity_spectrum[start:stop]})
    except Exception:
        details.update({'SPR angle': np.nan, 'SPR fit x': [np.nan], 'SPR fit y': [np.nan], 'SPR data x': angles, 'SPR data y': reflectivity_spectrum})
    return details


SENSORGRAM_TRACES = ['SPR angle', 'TIR angle', 'Bulk corrected']  # Traces of the sensorgram figure, in plotting order


def sensorgram_traces(sensorgram, S_SPR, S_TIR, layer_thickness, decay_length, drift=None):
    """
    Traces of the response quantification sensorgram (deg, not offset): SPR angle (drift corrected if drift is given),
    TIR angle and bulk corrected SPR angle,
        bulk corrected = SPR - TIR * S_SPR / S_TIR * exp(-2 d / delta)

    :param drift: optional dict with the SPR drift 'slope' (deg/min) that is removed (relative to the first time point)
    :return: dict {trace name: numpy array}, in the order of SENSORGRAM_TRACES
    """
    time_values = np.asarray(sensorgram['time'], dtype=np.float64)
    SPR_angle = np.asarray(sensorgram['SPR angle'], dtype=np.float64)
    TIR_angle = np.asarray(sensorgram['TIR angle'], dtype=np.float64)
    if drift:
        SPR_angle = SPR_angle - drift['slope'] * (time_values - time_values[0])
    try:
        bulk_factor = S_SPR / S_TIR * math.exp(-2 * layer_thickness / decay_length)
    except (TypeError, ZeroDivisionError):
        bulk_factor = 0
    return {'SPR angle': SPR_angle, 'TIR angle': TIR_angle, 'Bulk corrected': SPR_angle - TIR_angle * bulk_factor}


def surface_coverage(delta_theta, decay_length, dn_dC, S_SPR, thickness=0):
    """
    Surface coverage (ng/cm^2) of an adsorbed layer from an SPR angular shift:
        d > 0:  m_S = d * delta_theta_SPR / (dn/dC * S_SPR * (1 - exp(-2 d / delta)))
        d = 0:  m_S = delta * delta_theta_SPR / (2 * dn/dC * S_SPR)   (thin film approximation, the limit d -> 0)

    :param delta_theta: SPR angular shift (deg)
    :param decay_length: evanescent field decay length delta (nm)
    :param dn_dC: refractive index increment of the adsorbed material (cm^3/g, e.g. ~0.18 for proteins)
    :param S_SPR: bulk sensitivity of the SPR angle (deg/RIU)
    :param thickness: estimated layer thickness d (nm), 0 (or None) for the thin film approximation
    """
    # nm -> cm (1e-7) and g -> ng (1e9) gives a factor 100
    if not thickness or thickness <= 0:
        return 100.0 * decay_length * delta_theta / (2.0 * dn_dC * S_SPR)
    return 100.0 * thickness * delta_theta / (dn_dC * S_SPR * (1.0 - math.exp(-2.0 * thickness / decay_length)))


def linear_drift(time, signal, interval):
    """
    Slope (per minute) of a straight line fitted to the signal within interval = [start, end] (minutes).

    :return: slope, or None if fewer than 2 valid points lie in the interval
    """
    time = np.asarray(time, dtype=np.float64)
    signal = np.asarray(signal, dtype=np.float64)
    selection = (time >= min(interval)) & (time <= max(interval)) & np.isfinite(signal)
    if np.count_nonzero(selection) < 2:
        return None
    slope, _ = np.polyfit(time[selection], signal[selection], 1)
    return float(slope)
