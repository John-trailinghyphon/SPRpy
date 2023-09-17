# This is the main file where the webapp is initiated and further selections are made. It should ask for a datafile and
# either load .csv files directly or run conversion of .spr2 files using 'spr2_to_csv.py' in a separate thread
# (preferably showing a progress bar if possible)

#  Start with rewriting the dry scan fitting code and function files into python. Use scipy.optimize.least_squares
#  to replace MATLABs lsqnonlin.

#  Should have a class for modelling dry scan reflectivity traces and a second one for fitting reflectivity traces.
#  Each class should have methods for running the calculations using the fresnel_calculation() function and for saving
#  the results in an "experiment workbook". The optical parameters (with many defaults) and some naming strings should
#  be provided when an object is instanced from the class.

#  The workbook can be its own class with methods that define how data is stored and loaded for the dash app. The
#  idea is that this can be loaded by the app if a user wants to redo some modelling without starting all over again.
#  It should save a path to the datafile that was used for the analysis so that it has access to the data.


# TODO: Next, it is time to populate the Modelled/FittedReflectivityTrace classes with methods and attributes for
#  performing the fresnel_calculation() function. Note that the calculation results should be saved in the object, but
#  it is better to avoid saving the raw data, unless that is all that was used (simply plotting reflectivity trace for
#  instance).

# TODO: For the non-interacting probe method where a previously calculated background is used in the calculations, there
#  is a need of some way to select a previous analysis. The easiest way is probably to make sure that each analysis is
#  named something? I imagine that the user selects from a list of names of each analysis (and there should be a default
#  name that auto-increments with the object ID.

# Regarding the dash app below

# TODO: "Main control DIV". Need buttons for controlling a lot of things:
#  - Which measurement file to load data from (.csv)
#  - Choosing analysis method (default should be fitting reflectivity traces, like background)
#  - Session control (this should be done lasts, as it requires to figure out how to reinitiate the whole Dash interface with new values)
#     * loading previous session ("Continue from previous session")
#     * importing sensor and analysis objects from previous session (like background)
#     * saving session automatically
#     * removing sensor objects and analysis objects from the active session (this can be done by deleting their pickled files)
#  - Adding new analysis object (starting an analysis DIV)
#     * Selecting between different available types, which will change the analysis interface DIV and its options
#     * run calculations, fitting, plotting, selecting values, etc.
#  - Exporting the finished analysis as a HTML file with retained interactivity (omitting control DIV elements). This
#  can then be added to obsidian notes for instance, allowing straight forward result documentation.

# Plotly default discrete colors
# 1 '#636EFA',
# 2 '#EF553B',
# 3 '#00CC96',
# 4 '#AB63FA',
# 5 '#FFA15A',
# 6 '#19D3F3',
# 7 '#FF6692',
# 8 '#B6E880',
# 9 '#FF97FF',
# 10 '#FECB52'

import numpy as np
import datetime
import os
import tkinter
from tkinter.filedialog import askopenfilename, askdirectory
import pandas as pd
import dash
import dash_bootstrap_components as dbc
import copy
import plotly
import plotly.express as px
import plotly.graph_objects as go
import scipy
import fresnel_transfer_matrix as ftm
import re

# Configuration parameters

TIR_range_water_or_long_measurement = (60.8, 63)  # TIR range for water --> Automatically used for 50 or more scans per file
TIR_range_air_or_few_scans = (40.9, 41.8)  # TIR range for dry scans --> Automatically used for less than 50 scans per file


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name='Initial session', directory=os.getcwd()):
        self.name = datetime.datetime.now().__str__()[0:16] + ' ' + name
        if not os.path.exists(directory + r'\\pySPR sessions'):
            os.mkdir(directory + r'\\pySPR sessions')
        self.location = directory + r'\\pySPR sessions'
        self.sensor_instances = {}  # NOTE: The sessions in this list are also updated when modified as current sensor object
        self.sensor_ID = generate_id()
        self.analysis_instances = {}
        self.analysis_ID = generate_id()
        self.log = datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Welcome to pySPR!' \
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Start your session by defining your SPR sensor layers.' \
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'You can load a previous session or import previous results under "File and session controls".'

    def remove_sensor(self, sensor_object_id):
        """
        Remove a sensor object from the session.
        :return:
        """
        removed = self.sensor_instances.pop(sensor_object_id)
        print('Removed the following sensor object: ' + str(removed))

    def remove_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.analysis_instances.pop(analysis_object_id)
        print('Removed the following analysis object: ' + str(removed))


class Sensor:

    """
      An SPR measurement typically have some things in common, such as the sensor layers, measured angles,
    measured reflectivity, measurement time, etc. This information can be shared between different analysis methods for
    one measurement. This class serves as a basis for describing the current sensor, containing information about its
    layers and their optical properties.
    """

    def __init__(self, data_path_, object_id_, sensor_metal='Au', data_type='R', polarization=1):
        """
        :param data_path_: string
        :param sensor_metal: string, see options in method "set_default_optical_properties"
        :param polarization: int, default 1 or "p-polarization"
        """
        # Load sensor's default optical properties
        self.object_id = object_id_
        self.data_path = data_path_
        self.polarization = polarization
        self.data_type = data_type
        self.wavelength = int(data_path_[-9:-6])
        self.channel = data_path_[-12:-4].replace('_', ' ')
        self.sensor_metal = sensor_metal
        self.set_default_optical_properties(self.sensor_metal)
        self.fitted_var = self.optical_parameters.iloc[self.fitted_layer_index]

    def set_default_optical_properties(self, sensor_metal):

        # These default parameters should be set based on material layer and wavelength from loaded .csv file
        match sensor_metal:
            case 'Au' | 'gold' | 'Gold' | 'GOLD':
                self.layer_thicknesses = np.array([np.NaN, 2, 50, np.NaN])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([1.5202, 3.3105, 0.2238, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 3.9259, 0])
                    case 785:
                        self.refractive_indices = np.array([1.5162, 3.3225, 0.2580, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.6148, 4.88, 0])
                    case 980:
                        self.refractive_indices = np.array([1.5130, 3.4052, 0.28, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.5678, 6.7406, 0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Au', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'sio2' | 'SiO2' | 'SIO2' | 'Glass' | 'glass' | 'silica':
                # Fused silica values source: L. V. Rodríguez-de Marcos, J. I. Larruquert, J. A. Méndez, J. A. Aznárez.
                # Self-consistent optical constants of SiO2 and Ta2O5 films
                # Opt. Mater. Express 6, 3622-3637 (2016) (Numerical data kindly provided by Juan Larruquert)
                self.layer_thicknesses = np.array([np.NaN, 2, 50, 14, np.NaN])
                self.fitted_layer_index = (3, 1)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([1.5202, 3.3105, 0.2238, 1.4628, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 3.9259, 0, 0])
                    case 785:
                        self.refractive_indices = np.array([1.5162, 3.3225, 0.2580, 1.4610, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.6148, 4.88, 0, 0])
                    case 980:
                        self.refractive_indices = np.array([1.5130, 3.4052, 0.28, 1.4592, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.5678, 6.7406, 0, 0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Au', 'SiO2', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})
            case 'Pd' | 'palladium' | 'Palladium' | 'PALLADIUM':
                self.layer_thicknesses = np.array([np.NaN, 2, 20, np.NaN])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([1.5202, 3.3105, 2.25, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 4.60, 0])
                    case 785:
                        self.refractive_indices = np.array([1.5162, 3.3225, 2.5467, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.6148, 5.1250, 0])
                    case 980:
                        self.refractive_indices = np.array([1.5130, 3.4052, 3.0331, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.5678, 6.1010, 0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Pd', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'Pt' | 'platinum' | 'Platinum' | 'PLATINUM':
                self.layer_thicknesses = np.array([np.NaN, 2, 20, np.NaN])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([1.5202, 3.3105, 2.4687, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 5.2774, 0])
                    case 785:
                        print('WARNING! Default values for Pt, platinum, at 785 and 980 nm not yet supported. Enter values manually')
                        self.refractive_indices = np.array([1.5202, 3.3105, 2.4687, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 5.2774, 0])
                    case 980:
                        print('WARNING! Default values for Pt, platinum, at 785 and 980 nm not yet supported. Enter values manually')
                        self.refractive_indices = np.array([1.5202, 3.3105, 2.4687, 1.0003])
                        self.extinction_coefficients = np.array([0, 3.4556, 5.2774, 0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Pt', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

    def add_material_layer(self, thickness, n_re, n_im, layer_index_=-1):
        """
        Add additional layers on top of the sensor (before bulk medium). not used by dash app
        :return:
        """
        # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
        self.layer_thicknesses = np.insert(self.layer_thicknesses, layer_index_, thickness)
        self.refractive_indices = np.insert(self.refractive_indices, layer_index_, n_re)
        self.extinction_coefficients = np.insert(self.extinction_coefficients, layer_index_, n_im)
        self.fitted_layer = 'h_surf'
        print('Sensor thicknesses: ', self.layer_thicknesses)
        print('Sensor refractive indices: ', self.refractive_indices)
        print('Sensor extinction coefficients: ', self.extinction_coefficients)

    def remove_material_layer(self, layer_index_):

        """
        Removes a layer from a sensor. (Not used by dash app UI.)

        :param layer_index_: int, which layer to remove (starting from 1)
        :return:
        """

        self.layer_thicknesses = np.delete(self.layer_thicknesses, layer_index_-1, axis=0)
        if len(self.layer_thicknesses) == 4:
            self.fitted_layer = 'n_im_metal'
        self.refractive_indices = np.delete(self.refractive_indices, layer_index_-1, axis=0)
        self.extinction_coefficients = np.delete(self.extinction_coefficients, layer_index_-1, axis=0)

        print('Sensor thicknesses: ', self.layer_thicknesses)
        print('Sensor refractive indices: ', self.refractive_indices)
        print('Sensor extinction coefficients: ', self.extinction_coefficients)


class ModelledReflectivityTrace:
    """
    This class defines how a modelled reflectivity trace behaves. Note that a different sensor object should be used for
    each layer added to the sensor!

    TODO: Each object should also have a .csv export function.
    TODO: IN dash app, add functionality to update current sensor object with the model object optical parameters
    """

    def __init__(self, sensor_object, data_path_, TIR_range_, angle_range_, scanspeed_, object_id_):
        self.object_id = object_id_
        self.sensor_id = sensor_object.object_id
        self.polarization = sensor_object.polarization
        self.data_type = sensor_object.data_type
        self.wavelength = sensor_object.wavelength
        self.layer_thicknesses = sensor_object.layer_thicknesses
        self.fitted_var = sensor_object.fitted_var
        self.fitted_layer_index = sensor_object.fitted_layer_index
        self.refractive_indices = sensor_object.refractive_indices
        self.optical_parameters = sensor_object.optical_parameters
        self.extinction_coefficients = sensor_object.extinction_coefficients
        self.data_path = data_path_
        self.fit_result = None
        self.y_offset = 0
        self.TIR_range = TIR_range_  # Default for air. For water: (60.8, 63)
        self.angle_range = angle_range_
        self.scanspeed = scanspeed_  # Scan speed from .spr2 file, 1 (slow), 5 (medium) or 10 (fast)

    def calculate_fresnel_trace(self, angles_=np.linspace(39, 50, 1567)):

        fresnel_coefficients_ = ftm.fresnel_calculation(None,
                                                        angles=angles_,
                                                        fitted_layer_index=self.fitted_layer_index,
                                                        wavelength=self.wavelength,
                                                        layer_thicknesses=self.layer_thicknesses,
                                                        n_re=self.refractive_indices,
                                                        n_im=self.extinction_coefficients,
                                                        ydata=None,
                                                        ydata_type=self.data_type,
                                                        polarization=self.polarization)
        return fresnel_coefficients_

    def model_reflectivity_trace(self, ini_guess, bounds):
        """

        :param ini_guess:
        :param bounds:
        :param TIR_range:
        :param scanspeed:
        :return:
        """

        global current_data_path

        # Check if current data path matches data_path when object was first initialized, otherwise load previous data
        if current_data_path == self.data_path:
            global reflectivity_df
            xdata_ = reflectivity_df['angles']
            ydata_ = reflectivity_df['ydata']

        else:
            _, _, _, _, _, reflectivity_df_ = load_csv_data(path=self.data_path)
            xdata_ = reflectivity_df_['angles']
            ydata_ = reflectivity_df_['ydata']

        # Calculate TIR angle and bulk refractive index
        TIR_angle, TIR_fitted_angles, TIR_fitted_ydata = ftm.TIR_determination(xdata_, ydata_, self.TIR_range, self.scanspeed)
        self.refractive_indices[-1] = self.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)  # TODO: Currently, the sensor bulk RI in optical parameters do not update according to the TIR angle

        # Selecting a range of measurement data to use for fitting, and including an offset in reflectivity
        selection_xdata_ = xdata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])]
        selection_ydata_ = ydata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])] - self.y_offset

        # Perform the fitting
        result = scipy.optimize.least_squares(ftm.fresnel_calculation,
                                              ini_guess,
                                              bounds=bounds,
                                              kwargs={'fitted_layer_index': self.fitted_layer_index,
                                                      'wavelength': self.wavelength,
                                                      'layer_thicknesses': self.layer_thicknesses,
                                                      'n_re': self.refractive_indices,
                                                      'n_im': self.extinction_coefficients,
                                                      'angles': selection_xdata_,
                                                      'ydata': selection_ydata_,
                                                      'ydata_type': self.data_type,
                                                      'polarization': self.polarization}
                                              )
        # Collect the results from least_squares object and calculate corresponding fresnel coefficients
        self.fit_result = result['x']
        fresnel_coefficients = ftm.fresnel_calculation(self.fit_result,
                                                       fitted_layer_index=self.fitted_layer_index,
                                                       angles=selection_xdata_,
                                                       wavelength=self.wavelength,
                                                       layer_thicknesses=self.layer_thicknesses,
                                                       n_re=self.refractive_indices,
                                                       n_im=self.extinction_coefficients,
                                                       ydata=None,
                                                       ydata_type='R',
                                                       polarization=1
                                                       )

        return self.fit_result, selection_xdata_, fresnel_coefficients

    def export_results(self):
        """
        Exporting the result (including parameters) of a particular analysis as a .csv file
        :return:
        """
        pass


class FittedReflectivityTrace(ModelledReflectivityTrace):

    """

    """

    def __init__(self, sensor_object, data_path_, object_id_, ydata_type='R'):
        super().__init__(sensor_object, data_path_, sensor_object.object_id)  # Initializes the same way as parent objects, to shorten code
        self.object_id = object_id_
        self.ydata_type = ydata_type

    def calculate_fit(self):
        pass


class NonInteractingProbe(ModelledReflectivityTrace):

    """

    """

    def __init__(self, sensor_object, data_path_, object_id_, ydata_type='R'):
        super().__init__(sensor_object, data_path_, sensor_object.object_id)  # Initializes the same way as parent objects, to shorten code
        self.object_id = object_id_
        self.ydata_type = ydata_type

    def calculate_fit(self):
        pass


def add_sensor_backend(session_object, data_path_, sensor_metal='Au', polarization=1):

    """
    Adds sensor objects to a session object.
    :return: a sensor object
    """
    id_ = next(session_object.sensor_ID)
    sensor_object = Sensor(data_path_, id_, sensor_metal=sensor_metal, polarization=polarization)
    session_object.sensor_instances[id_] = sensor_object

    return sensor_object


def copy_sensor_backend(session_object, sensor_object):

    """
    Copies sensor object to a session object.
    :return: a sensor object
    """
    id_ = next(session_object.sensor_ID)
    copied_sensor_object = copy.deepcopy(sensor_object)
    copied_sensor_object.object_id = id_
    session_object.sensor_instances[id_] = copied_sensor_object

    return copied_sensor_object



def add_modelled_reflectivity_trace(session_object, sensor_object, data_path_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_object.analysis_ID)
    analysis_object = ModelledReflectivityTrace(sensor_object, data_path_, id_)
    session_object.analysis_instances[id_] = analysis_object

    return analysis_object


def add_fitted_reflectivity_trace(session_object, sensor_object, data_path_, ydata_type='R'):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_object.analysis_ID)
    analysis_object = FittedReflectivityTrace(sensor_object, data_path_, id_, ydata_type=ydata_type)
    session_object.analysis_instances[id_] = analysis_object

    return analysis_object


def load_session(filename):
    """
    Loads a previous session.
    :return:
    """
    # TODO: ADD DASH CALLBACK FUNCTIONALITY. It should load a session object and load its dictionaries containing
    #  measurement_instances and analysis_instances.

    pass


def save_to_session(session_handle, object_handle):
    """
    Save a session to a binary pickle file.
    :return:
    """
    # TODO: Save a measurement or analysis object to a session
    pass


def import_to_session(past_session, current_session, sensor_import_id=None, analysis_import_id=None):
    """

    :return:
    """
    pass


def generate_id():
    """
    Each time this function is called within an object it can return a new ID.
    :yield:
    """
    new_id = 0
    while True:
        yield new_id
        new_id += 1


def load_csv_data(path=False):
    if not path:
        print('Select the measurement data file (.csv)')
        root = tkinter.Tk()
        root.attributes("-topmost", 1)
        root.withdraw()
        data_path_ = askopenfilename(title='Select the measurement data file', filetypes=[('CSV files', '*.csv')],
                                     initialdir=r'C:\Users\anjohn\OneDrive - Chalmers\Dahlin group\Data\SPR')
        root.destroy()
    else:
        data_path_ = path

    #  Determine the scanning speed/step length if present in the file
    try:
        with open(data_path_, 'r') as file:
            step_length_pattern = re.compile(r'=\d{1,2}')
            scanspeed = int(step_length_pattern.search(file.readline()).group().strip('='))

    except:
        scanspeed = 5

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
            TIR_theta, _, _ = ftm.TIR_determination(angles, reflectivity_spectrum, TIR_range, scanspeed)
            sensorgram_TIR_angles[ind-1] = TIR_theta

        except:
            print('No TIR found. Skipping measurement point...')
            sensorgram_TIR_angles[ind-1] = np.NaN

    sensorgram_df = pd.DataFrame(data={'time': time, 'SPR angle': sensorgram_SPR_angles, 'TIR angle': sensorgram_TIR_angles})

    return sensorgram_df

# def save_new_measurement(self):
#     """
#     Saving a new measurement instance for the session.
#     :return:
#     """
#     pass
#
# def overwrite_measurement(self):
#     """
#     Overwriting a measurement instance for the session (does not produce a new instance).
#     :return:
#     """
#     pass


if __name__ == '__main__':

    # Create initial session
    current_session = Session()

    # Prompt user for initial measurement data
    current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data()

    # Calculate sensorgram from loaded data (assume air or liquid medium for TIR calculation based on number of scans)
    if ydata_df.shape[0] > 50:
        TIR_range = TIR_range_water_or_long_measurement
    else:
        TIR_range = TIR_range_air_or_few_scans

    sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, TIR_range, scanspeed)

    # Offset to start at 0 degrees at 0 minutes
    sensorgram_df_selection = sensorgram_df
    sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - sensorgram_df_selection['SPR angle'][0]
    sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - sensorgram_df_selection['TIR angle'][0]

    # Add sensor object based on chosen measurement data
    current_sensor = add_sensor_backend(current_session, current_data_path)

    # Dash app
    app = dash.Dash(external_stylesheets=[dbc.themes.SPACELAB])

    # Dash figures
    reflectivity_fig = px.line(reflectivity_df, x='angles', y='ydata')
    reflectivity_fig.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                   yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                   font_family='Balto',
                                   font_size=19,
                                   margin_r=25,
                                   margin_l=60,
                                   margin_t=40,
                                   template='simple_white')
    reflectivity_fig.update_xaxes(mirror=True, showline=True)
    reflectivity_fig.update_yaxes(mirror=True, showline=True)

    sensorgram_fig = px.line(sensorgram_df_selection, x='time', y='SPR angle')
    sensorgram_fig['data'][0]['showlegend'] = True
    sensorgram_fig['data'][0]['name'] = 'SPR angle'
    sensorgram_fig.add_trace(go.Scatter(x=sensorgram_df_selection['time'],
                                        y=sensorgram_df_selection['TIR angle'],
                                        name='TIR angle'))
    sensorgram_fig.update_layout(xaxis_title=r'$\large{\text{Time [min]}}$',
                                 yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                                 font_family='Balto',
                                 font_size=19,
                                 margin_r=25,
                                 margin_l=60,
                                 margin_t=40,
                                 template='simple_white',
                                 clickmode='event+select')
    sensorgram_fig.update_xaxes(mirror=True, showline=True)
    sensorgram_fig.update_yaxes(mirror=True, showline=True)

    # Dash webapp layout
    app.layout = dash.html.Div([

        # Heading for page
        dbc.Container(
            [
                dbc.Card(
                    [
                        dbc.CardImg(src='static/images/SPR_principle.svg', top=True),
                        # dbc.CardBody([dash.html.H4('Surface plasmon resonance (SPR)', className='card-title')])
                    ], style={'width': '22rem'}
                ),
                dbc.Card(
                    [
                        dbc.CardImg(src='static/images/fresnel_material.svg', top=True),
                        # dbc.CardBody([dash.html.H4('Fresnel modelling', className='card-title')])
                    ], style={'width': '19rem', 'padding-top': '30px', 'margin-left': '2rem'}
                ),
                dash.dcc.Markdown('''
                # **#pySPR#**
                ''', className='dash-bootstrap', style={'margin-top': '6rem', 'margin-left': '5rem', 'margin-right': '5rem'}),
                dbc.Card(
                    [
                        dbc.CardImg(src='static/images/SPR_angular_spectrum.svg', top=True),
                        # dbc.CardBody([dash.html.H4('SPR sensor', className='card-title')])
                    ], style={'width': '23rem', 'padding-top': '18px', 'margin-right': '2rem'}
                ),
                dbc.Card(
                    [
                        dbc.CardImg(src='static/images/non-interacting_height_probe.PNG', top=True),
                        # dbc.CardBody([dash.html.H4('Non-interacting height probing', className='card-title')])
                    ], style={'width': '17rem', 'padding-top': '20px'}
                ),
            ], style={'margin-top': '20px', 'display': 'flex', 'justify-content': 'space-between'}
        ),

        # Session log div
        dash.html.Div([
            dash.html.H3("Session log", className='dash-bootstrap'),
            dash.dcc.Textarea(
                id='console',
                value=current_session.log,
                readOnly=True,
                className='dash-bootstrap',
                style={'width': '98%', 'height': '150px', 'margin-right': '2%'}
            )
        ], style={'margin-top': '40px', 'margin-left': '2%', 'text-align': 'left'}),

        # Button for adding note to session log
        dash.html.Div([
            dbc.InputGroup(
                [
                    dbc.Button('Add note to log', id='submit-button', n_clicks=0, color='info'),
                    dbc.Input(id='test-input', value='', type='text', style={'margin-right': '2%'})
                ]
            )

        ], style={'margin-left': '2%'}),

        # File and session control
        dash.html.H3("File and session controls", className='dash-bootstrap', style={'margin-top': '20px', 'text-align': 'center'}),
        dash.html.Div(['Current measurement file:    ', current_data_path.split('/')[-1]],
                      id='datapath-textfield',
                      style={'margin-right': '10px', 'textAlign': 'center'}),
        dbc.Container([
            dbc.ButtonGroup([
                dbc.Button('Load data',
                           id='load-data',
                           n_clicks=0,
                           color='primary',
                           title='Load data from another measurement. Analysis is always performed on this active measurement'),
                dbc.Button('Load session',  # TODO: It is probably best to prompt for this outside the dash app, if at all.
                           id='load-session',
                           n_clicks=0,
                           color='primary',
                           title='Load a previous session in its entirety'),
                dbc.Button('Import result',
                           id='import-from-session',
                           n_clicks=0,
                           color='primary',
                           title='Use this to import previous results from another session'),
                dbc.DropdownMenu(
                    id='create-new-sensor-dropdown',
                    label='Add new sensor',
                    color='primary',
                    children=[dbc.DropdownMenuItem('Gold', id='new-sensor-gold', n_clicks=0),
                              dbc.DropdownMenuItem('Glass', id='new-sensor-glass', n_clicks=0),
                              dbc.DropdownMenuItem('Palladium', id='new-sensor-palladium', n_clicks=0),
                              dbc.DropdownMenuItem('Platinum', id='new-sensor-platinum', n_clicks=0)], style={'margin-left': '-5px'}),
                dbc.Button('Copy current sensor',
                           id='copy-sensor',
                           n_clicks=0,
                           color='primary',
                           title='Use this to copy a sensor table\'s values into a new sensor'),
                dbc.DropdownMenu(
                    id='chosen-sensor-dropdown',
                    label='Sensors',
                    color='primary',
                    children=[
                        dbc.DropdownMenuItem('Sensor ' + str(sensor_id), id={'type': 'sensor', 'index': sensor_id},
                                             n_clicks=0) for sensor_id in current_session.sensor_instances], style={'margin-left': '-5px'})
            ])
        ], style={'margin-bottom': '20px', 'display': 'flex', 'justify-content': 'center'}),

        # Sensor datatable
        dash.html.Div([
            dash.html.Div([
                dash.html.H4(['Sensor {sensor_number} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                    sensor_number=current_sensor.object_id,
                    channel=current_sensor.channel,
                    fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                    fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])
                              ], id='sensor-table-title', style={'text-align': 'center'}),
                dash.html.Div([
                    dash.dash_table.DataTable(data=current_sensor.optical_parameters.to_dict('records'),
                                              columns=[{'name': col, 'id': col} for col in
                                                       current_sensor.optical_parameters.columns],
                                              editable=True,
                                              row_deletable=True,
                                              cell_selectable=True,
                                              id='sensor-table',
                                              style_header={
                                                  'backgroundColor': '#446e9b',
                                                  'color': 'white',
                                                  'fontWeight': 'bold'
                                              },
                                              style_cell={'textAlign': 'center'}),
                ], style={'margin-left': '6px'}),
                dbc.ButtonGroup([
                    dbc.Button('Add layer',
                               id='add-table-layer',
                               n_clicks=0,
                               color='primary',
                               title='Add a new layer on the sensor surface'),
                    dbc.Button('Update table',
                               id='table-update-values',
                               n_clicks=0,
                               color='danger',
                               title='Refresh the table after editing its values'),
                    dbc.Button('Select fitted variable',
                               id='table-select-fitted',
                               n_clicks=0,
                               color='success',
                               title='Click this button after selecting a different parameter to fit by clicking it such'
                                     ' that it is marked in red.'),
                    dbc.Button(
                        "Show default values",
                        id="show-default-param-button",
                        color="secondary",
                        n_clicks=0,
                        title='CTRL+Z not supported, check default values here.'
                    ),
                ], style={'width': '672px', 'margin-left': '4px', 'margin-top': '5px', 'margin-bottom': '20px'}),
            ], style={'width': '675px'}),
            dash.html.Div([
                dbc.Collapse(
                    dbc.Card(
                        dbc.CardBody(
                            dbc.Table.from_dataframe(pd.DataFrame(
                                {
                                    "Layer": ['Prism', 'Cr', 'Au', 'SiO2', 'Pd', 'Pt'],
                                    "d[nm]": ['', '2', '50', '14', '20', '20'],
                                    "n (670)": ['1.5202', '3.3105', '0.2238', '1.4628', '2.2500', '2.4687'],
                                    "n (785)": ['1.5162', '3.3225', '0.2580', '1.4610', '2.5467', '?'],
                                    "n (980)": ['1.5130', '3.4052', '0.2800', '1.4592', '3.0331', '?'],
                                    "k (670)": ['0', '3.4556', '3.9259', '0', '4.6000', '5.2774'],
                                    "k (785)": ['0', '3.6148', '4.8800', '0', '5.1250', '?'],
                                    "k (980)": ['0', '3.5678', '6.7406', '0', '6.7406', '?'],
                                }
                            ), size='sm', striped=True, bordered=True, hover=True)
                        ), style={'width': '650px'}),
                    id='default-values-collapse',
                    is_open=False)
            ], style={'margin-top': '40px', 'margin-left': '10px'}),
        ], style={'display': 'flex', 'justify-content': 'center'}),

        # Analysis tabs
        dash.html.Div([
            dash.html.H1(['Analysis options']),
            dbc.Tabs([

                # Data plotting tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.dcc.Graph(id='plotting-angular-reflectivity-graph',
                                           figure=reflectivity_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.Button('Add trace',
                                           id='plotting-reflectivity-add-trace',
                                           n_clicks=0,
                                           color='warning',
                                           title='Add a trace to the figure from an external dry scan .csv file. The most recent scan in the file is used.'),
                                dbc.DropdownMenu(
                                    id='reflectivity-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[
                                        dbc.DropdownMenuItem('.PNG', id='plotting-reflectivity-save-png', n_clicks=0),
                                        dbc.DropdownMenuItem('.SVG', id='plotting-reflectivity-save-svg', n_clicks=0),
                                        dbc.DropdownMenuItem('.HTML', id='plotting-reflectivity-save-html', n_clicks=0)],
                                    style={'margin-left': '-5px'})
                            ], style={'margin-left': '30%'}),
                        ], style={'width': '35%'}),
                        dash.html.Div([
                            dash.dcc.Graph(id='plotting-sensorgram-graph',
                                           figure=sensorgram_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.DropdownMenu(
                                    id='sensorgram-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[dbc.DropdownMenuItem('.PNG', id='plotting-sensorgram-save-png', n_clicks=0),
                                              dbc.DropdownMenuItem('.SVG', id='plotting-sensorgram-save-svg', n_clicks=0),
                                              dbc.DropdownMenuItem('.HTML', id='plotting-sensorgram-save-html', n_clicks=0)],
                                    style={'margin-left': '-5px'}),
                                ], style={'margin-left': '40%'}),
                        ], style={'width': '60%'})
                    ], id='plotting-tab-content', style={'display': 'flex', 'justify-content': 'center'})
                ], label='Data plotting', tab_id='plotting-tab', style={'margin-top': '10px'}),

                # Fresnel modelling tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.dcc.Graph(id='fresnel-angular-reflectivity-graph',
                                           figure=reflectivity_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.Button('Run modelling',
                                           id='fresnel-reflectivity-run-model',
                                           n_clicks=0,
                                           color='warning',
                                           title='Run the fresnel model'),
                                dbc.DropdownMenu(
                                    id='fresnel-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[
                                        dbc.DropdownMenuItem('.PNG', id='fresnel-reflectivity-save-png', n_clicks=0),
                                        dbc.DropdownMenuItem('.SVG', id='fresnel-reflectivity-save-svg', n_clicks=0),
                                        dbc.DropdownMenuItem('.HTML', id='fresnel-reflectivity-save-html', n_clicks=0)],
                                    style={'margin-left': '-5px'})
                            ], style={'margin-left': '30%'}),
                        ], style={'width': '35%'}),
                        dash.html.Div([
                            dash.html.H3(['Fitting options']),
                            dbc.Form([
                                dash.html.Div([
                                    dbc.ButtonGroup([
                                        dbc.Button('Add new analysis',
                                                   id='fresnel-add-analysis-button',
                                                   n_clicks=0,
                                                   color='primary',
                                                   title='Add a new analysis object for the current sensor.'),
                                        dbc.DropdownMenu(id='fresnel-analysis-dropdown',
                                                         label='Chose analysis',
                                                         color='primary',
                                                         children=[])
                                    ])
                                ]),
                                dash.html.Div([
                                    dbc.Collapse(
                                        dbc.Card(
                                            dbc.CardBody(
                                                dbc.Form([
                                                    #  TIR_range, scanspeed (for TIR)
                                                    #  polarization, angle_range, ini_guess, upper_lowerbound, offset, weight_factor
                                                    dbc.Row([
                                                        dbc.Label('Initial guess', width='auto'),
                                                        dbc.Col([
                                                            dbc.Input(id='fresnel-fit-option-iniguess',
                                                                      value=current_sensor.fitted_var, type='number')
                                                        ], width=2),
                                                        dbc.Label('Bounds', width='auto'),
                                                        dbc.Col([
                                                            dbc.InputGroup([
                                                                dbc.Input(id='fresnel-fit-option-lowerbound',
                                                                          value=current_sensor.fitted_var - current_sensor.fitted_var / 2,
                                                                          type='number'),
                                                                dbc.Input(id='fresnel-fit-option-upperbound',
                                                                          value=current_sensor.fitted_var + current_sensor.fitted_var / 2,
                                                                          type='number')
                                                            ])
                                                        ], width=4)
                                                    ], style={'margin-bottom': '10px'}),
                                                    dbc.Row([
                                                        dbc.Label('Angle range', width='auto'),
                                                        dbc.Col([
                                                            dash.dcc.RangeSlider(40, 80,
                                                                                 marks={'40': '40', '45': '45',
                                                                                        '50': '50', '55': '55',
                                                                                        '60': '60', '65': '65',
                                                                                        '70': '70', '75': '75',
                                                                                        '80': '80'},
                                                                                 step=0.1,
                                                                                 allowCross=False,
                                                                                 tooltip={"placement": "top",
                                                                                          "always_visible": True},
                                                                                 id='fresnel-fit-options-rangeslider')
                                                        ])
                                                    ], style={'margin-bottom': '10px'}),
                                                    dbc.Row([
                                                        dbc.Label('Extinction correction', width='auto'),
                                                        dbc.Col([
                                                            dash.dcc.Slider(min=0, max=0.1,
                                                                            step=0.005,
                                                                            marks={0: '0', 0.01: '0.01',
                                                                                   0.02: '0.02', 0.03: '0.03',
                                                                                   0.04: '0.04', 0.05: '0.05',
                                                                                   0.06: '0.06', 0.07: '0.07',
                                                                                   0.08: '0.08', 0.09: '0.09',
                                                                                   0.1: '0.1'},
                                                                            tooltip={"placement": "top",
                                                                                     "always_visible": True},
                                                                            id='fresnel-fit-options-extinctionslider')
                                                        ])
                                                    ], style={'margin-bottom': '10px'}),
                                                ])
                                            )
                                        ), id='fresnel-analysis-options-collapse', is_open=True)
                                ])


                            ], id='fresnel-fit-options-form')
                        ], style={'margin-top': '1.9rem'})
                    ], id='fresnel-tab-content', style={'display': 'flex', 'justify-content': 'center'})
                ], label='Fresnel modelling', tab_id='fresnel-tab', style={'margin-top': '10px'}),

                # Non-interacting height probe tab
                dbc.Tab([
                    dash.html.Div(
                        ['Non-interacting height probe'],
                        id='probe-tab-content')
                ], label='Exclusion height probing', tab_id='probe-tab', style={'margin-top': '10px'}),

                # Result summary tab
                dbc.Tab([
                    dash.html.Div(
                        ['Summary'],
                        id='summary-tab-content')
                ], label='Result summary', tab_id='summary-tab', style={'margin-top': '10px'}),
            ], id='analysis-tabs', active_tab='fresnel-tab'),

        ], style={'margin-left': '2%', 'margin-right': '2%'})

    ])

    # Adding note to session log
    @dash.callback(
        dash.Output('console', 'value'),
        dash.Input('submit-button', 'n_clicks'),
        dash.State('console', 'value'),
        dash.State('test-input', 'value')
    )
    def update_session_log(input1, state1, state2):
        new_message = state1 + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + state2
        current_session.log = new_message
        return new_message

    # Adding new sensor
    @dash.callback(
        dash.Output('chosen-sensor-dropdown', 'children'),  # Update chosen sensor dropdown
        dash.Input('new-sensor-gold', 'n_clicks'),
        dash.Input('new-sensor-glass', 'n_clicks'),
        dash.Input('new-sensor-palladium', 'n_clicks'),
        dash.Input('new-sensor-platinum', 'n_clicks'),
        dash.Input('copy-sensor', 'n_clicks'),
    )
    def add_sensor_UI(input1, input2, input3, input4, copy):
        """
        Dictates what happens when creating a new sensor under the "Add new sensor" dropdown menu.

        :param input1: Adding gold sensor
        :param input2: Adding glass sensor
        :param input3: Adding palladium sensor
        :param input4: Adding platinum sensor
        :return: New options to the sensor list
        """

        global current_sensor
        global current_session
        global current_data_path

        if 'new-sensor-gold' == dash.ctx.triggered_id:
            add_sensor_backend(current_session, current_data_path, sensor_metal='Au')
        elif 'new-sensor-glass' == dash.ctx.triggered_id:
            add_sensor_backend(current_session, current_data_path, sensor_metal='SiO2')
        elif 'new-sensor-palladium' == dash.ctx.triggered_id:
            add_sensor_backend(current_session, current_data_path, sensor_metal='Pd')
        elif 'new-sensor-platinum' == dash.ctx.triggered_id:
            add_sensor_backend(current_session, current_data_path, sensor_metal='Pt')
        elif 'copy-sensor' == dash.ctx.triggered_id:
            copy_sensor_backend(current_session, current_sensor)

        sensor_options = [dbc.DropdownMenuItem('Sensor ' + str(sensor_id), id={'type': 'sensor-list', 'index': sensor_id},
                                               n_clicks=0) for sensor_id in current_session.sensor_instances]

        return sensor_options

    # Updating the sensor table with new values and properties
    @dash.callback(
        dash.Output('sensor-table', 'data'),  # Update sensor table data
        dash.Output('sensor-table-title', 'children'),  # Update sensor table title
        dash.Input({'type': 'sensor-list', 'index': dash.ALL}, 'n_clicks'),
        dash.Input('add-table-layer', 'n_clicks'),
        dash.Input('table-update-values', 'n_clicks'),
        dash.Input('table-select-fitted', 'n_clicks'),
        dash.State('sensor-table', 'data'),
        dash.State('sensor-table', 'columns'),
        dash.State('sensor-table', 'active_cell'),
        prevent_initial_call=True)
    def update_sensor_table(n_clicks_sensor_list, n_clicks_add_row, n_clicks_update, n_clicks_fitted, table_rows, table_columns, active_cell):
        """
        This callback function controls all updates to the sensor table.

        :param n_clicks_sensor_list: Choose sensor dropdown menu
        :param n_clicks_add_row: Add layers button
        :param n_clicks_update: Update table values button
        :param n_clicks_fitted: Update fitted variable
        :param table_rows: Data rows (state)
        :param table_columns: Column names (state)
        :param active_cell: Dict with columns and rows of highlighted cell (state)

        :return: Updated data rows in sensor table and the sensor table title
        """

        global current_sensor

        if 'table-update-values' == dash.ctx.triggered_id:

            # Update background sensor object with new row
            current_sensor.optical_parameters = pd.DataFrame.from_records(table_rows)
            current_sensor.layer_thicknesses = current_sensor.optical_parameters['d [nm]'].to_numpy()
            current_sensor.refractive_indices = current_sensor.optical_parameters['n'].to_numpy()
            current_sensor.extinction_coefficients = current_sensor.optical_parameters['k'].to_numpy()
            current_sensor.fitted_var = current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index]

            return table_rows, dash.no_update

        elif 'add-table-layer' == dash.ctx.triggered_id:
            table_rows.insert(-1, {c['id']: '' for c in table_columns})

            return table_rows, dash.no_update

        elif 'table-select-fitted' == dash.ctx.triggered_id:

            current_sensor.fitted_layer_index = (active_cell['row'], active_cell['column'])
            current_sensor.fitted_var = current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index]
            sensor_table_title = 'Sensor {sensor_number} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[active_cell['row'], 0],
                fitted_param=current_sensor.optical_parameters.columns[active_cell['column']])

            return dash.no_update, sensor_table_title

        else:
            current_sensor = current_session.sensor_instances[dash.callback_context.triggered_id.index]

            data_rows = current_sensor.optical_parameters.to_dict('records')
            sensor_table_title = 'Sensor {sensor_number} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            return data_rows, sensor_table_title

    # Toggle view of default optical parameters for different materials
    @dash.callback(
        dash.Output('default-values-collapse', 'is_open'),
        dash.Input('show-default-param-button', 'n_clicks'),
        dash.State('default-values-collapse', 'is_open')
    )
    def show_default_parameters(n_clicks, is_open):

        if n_clicks:
            return not is_open

        return is_open

    # Update the reflectivity plot in the Data plotting tab
    @dash.callback(
        dash.Output('plotting-angular-reflectivity-graph', 'figure'),
        dash.Input('plotting-reflectivity-add-trace', 'n_clicks'),
        dash.Input('plotting-reflectivity-save-png', 'n_clicks'),
        dash.Input('plotting-reflectivity-save-svg', 'n_clicks'),
        dash.Input('plotting-reflectivity-save-html', 'n_clicks'),
        dash.Input('plotting-sensorgram-graph', 'hoverData'),
        dash.State('plotting-angular-reflectivity-graph', 'figure'),
    )
    def update_reflectivity_plotting_graph(add_trace, save_png, save_svg, save_html, hoverData, figure_JSON):

        figure_object = go.Figure(figure_JSON)

        # Update based on hover over sensorgram figure
        if 'plotting-sensorgram-graph' == dash.ctx.triggered_id:

            # First make sure no other traces has been added and the very first value is ignored
            if figure_object.data.__len__() == 1:

                time_index = hoverData['points'][0]['pointIndex']

                new_trace_data = ydata_df.loc[time_index+1]

                new_figure = go.Figure(go.Scatter(x=angles_df,
                                                  y=new_trace_data,
                                                  mode='lines',
                                                  showlegend=False,
                                                  line_color='#636efa'))
                new_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                         yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                         font_family='Balto',
                                         font_size=19,
                                         margin_r=25,
                                         margin_l=60,
                                         margin_t=40,
                                         template='simple_white',
                                         uirevision=True)
                new_figure.update_xaxes(mirror=True,
                                        showline=True)
                new_figure.update_yaxes(mirror=True,
                                        showline=True)

                return new_figure

            else:
                return dash.no_update

        # This adds a trace to the reflectivity plot from a separate measurement file. The trace data is not stored.
        elif 'plotting-reflectivity-add-trace' == dash.ctx.triggered_id:

            _, _, _, _, _, trace_reflectivity_df = load_csv_data()
            figure_object.add_trace(go.Scatter(x=trace_reflectivity_df['angles'],
                                               y=trace_reflectivity_df['ydata'],
                                               mode='lines',
                                               showlegend=False))

        elif 'plotting-reflectivity-save-html' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_html(figure_object, save_folder+r'\reflectivity_plot.html', include_mathjax='cdn')

        elif 'plotting-reflectivity-save-svg' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder+r'\reflectivity_plot.svg', format='svg')

        elif 'plotting-reflectivity-save-png' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder+r'\reflectivity_plot.png', format='png')

        return figure_object

    # Update the sensorgram plot in the Data plotting tab
    @dash.callback(
        dash.Output('plotting-sensorgram-graph', 'figure'),
        dash.Input('plotting-sensorgram-save-png', 'n_clicks'),
        dash.Input('plotting-sensorgram-save-svg', 'n_clicks'),
        dash.Input('plotting-sensorgram-save-html', 'n_clicks'),
        dash.Input('plotting-sensorgram-graph', 'clickData'),
        dash.State('plotting-sensorgram-graph', 'figure'),
        prevent_initial_call=True)  # Adding this fixed a weird bug with graph not updating after firing clickData callbacks
    def update_sensorgram_plotting_tab(save_png, save_svg, save_html, clickData, figure_JSON):

        figure_object = go.Figure(figure_JSON)

        if 'plotting-sensorgram-graph' == dash.ctx.triggered_id:
            global sensorgram_df_selection

            offset_index = clickData['points'][0]['pointIndex']

            new_sensorgram_fig = go.Figure(go.Scatter(x=sensorgram_df_selection['time'],
                                                      y=sensorgram_df_selection['SPR angle']-sensorgram_df_selection['SPR angle'].loc[offset_index],
                                                      name='SPR angle',
                                                      line_color='#636efa'))

            new_sensorgram_fig.add_trace(go.Scatter(x=sensorgram_df_selection['time'],
                                                    y=sensorgram_df_selection['TIR angle']-sensorgram_df_selection['TIR angle'].loc[offset_index],
                                                    name='TIR angle',
                                                    line_color='#ef553b'))

            new_sensorgram_fig.update_layout(xaxis_title=r'$\large{\text{Time [min]}}$',
                                             yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                                             font_family='Balto',
                                             font_size=19,
                                             margin_r=25,
                                             margin_l=60,
                                             margin_t=40,
                                             template='simple_white',
                                             uirevision=True)
            new_sensorgram_fig.update_xaxes(mirror=True, showline=True)
            new_sensorgram_fig.update_yaxes(mirror=True, showline=True)

            return new_sensorgram_fig

        elif 'plotting-sensorgram-save-html' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_html(figure_object, save_folder + r'\reflectivity_plot.html', include_mathjax='cdn')

            return figure_object

        elif 'plotting-sensorgram-save-svg' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder + r'\reflectivity_plot.svg', format='svg')

            return figure_object

        elif 'plotting-sensorgram-save-png' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder + r'\reflectivity_plot.png', format='png')

            return figure_object

    # Update the reflectivity plot in the Fresnel fitting tab
    @dash.callback(
        dash.Output('fresnel-angular-reflectivity-graph', 'figure'),
        dash.Input('fresnel-reflectivity-run-model', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-png', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-svg', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-html', 'n_clicks'),
        dash.Input('fresnel-fit-options-rangeslider', 'value'),
        dash.State('fresnel-angular-reflectivity-graph', 'figure'),
        dash.State('fresnel-fit-options-rangeslider', 'value'),
        dash.State('fresnel-fit-options-iniguess', 'value'),
        dash.State('fresnel-fit-options-lowerbound', 'value'),
        dash.State('fresnel-fit-options-upperbound', 'value'),
        dash.State('fresnel-fit-options-extinctionslider', 'value'),
        )
    def update_reflectivity_fresnel_graph(run_model, save_png, save_svg, save_html, rangeslider_inp,
                                          figure_JSON, rangeslider_state, ini_guess, lower_bound, upper_bound,
                                          extinction_correction):

        figure_object = go.Figure(figure_JSON)

        if 'fresnel-fit-options-rangeslider' == dash.ctx.triggered_id:

            # First check if model has been run previously, then include model data before adding angle range lines
            if figure_object.data.__len__() > 3:
                new_figure = go.Figure(go.Scatter(x=figure_object.data[0]['x'],
                                                  y=figure_object.data[0]['y'],
                                                  mode='lines',
                                                  showlegend=False,
                                                  line_color='#636efa'
                                                  ))
                new_figure.add_trace(go.Scatter(x=figure_object.data[1]['x'],
                                                y=figure_object.data[1]['y'],
                                                mode='lines',
                                                showlegend=False,
                                                line_color='#ef553b'
                                                ))
            else:
                new_figure = go.Figure(go.Scatter(x=figure_object.data[0]['x'],
                                                  y=figure_object.data[0]['y'],
                                                  mode='lines',
                                                  showlegend=False,
                                                  line_color='#636efa'
                                                  ))
            # Adding angle range lines
            new_figure.add_trace(go.Scatter(x=[rangeslider_inp[0], rangeslider_inp[0]],
                                            y=[min(figure_object.data[0]['y']), max(figure_object.data[0]['y'])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            new_figure.add_trace(go.Scatter(x=[rangeslider_inp[1], rangeslider_inp[1]],
                                            y=[min(figure_object.data[0]['y']), max(figure_object.data[0]['y'])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            # Updating layout
            new_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                     yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                     font_family='Balto',
                                     font_size=19,
                                     margin_r=25,
                                     margin_l=60,
                                     margin_t=40,
                                     template='simple_white',
                                     uirevision=True)
            new_figure.update_xaxes(mirror=True,
                                    showline=True)
            new_figure.update_yaxes(mirror=True,
                                    showline=True)
            return new_figure

        # TODO: Incorporate fresnel fitting backend functionality into the "Run modelling" button
        elif 'fresnel-reflectivity-run-model' == dash.ctx.triggered_id:
            fresnel_figure = None
            return fresnel_figure

        elif 'fresnel-reflectivity-save-html' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_html(figure_object, save_folder + r'\fresnel_plot.html', include_mathjax='cdn')

        elif 'fresnel-reflectivity-save-svg' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder + r'\fresnel_plot.svg', format='svg')

        elif 'fresnel-reflectivity-save-png' == dash.ctx.triggered_id:
            root = tkinter.Tk()
            root.attributes("-topmost", 1)
            root.withdraw()
            save_folder = askdirectory(title='Choose folder', parent=root)
            root.destroy()
            plotly.io.write_image(figure_object, save_folder + r'\fresnel_plot.png', format='png')

        return figure_object


    app.run_server(debug=True, use_reloader=False)
