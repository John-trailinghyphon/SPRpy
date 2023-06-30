# This is the main file where the webapp is initiated and further selections are made. It should ask for a datafile and
# either load .csv files directly or run conversion of .spr2 files using 'extract_SPR_spectra.py' in a separate thread
# (preferably showing a progress bar if possible)

# TODO: Start with rewriting the dry scan fitting code and function files into python. Use scipy.optimize.least_squares
#  to replace MATLABs lsqnonlin.

# TODO: Should have a class for modelling dry scan reflectivity traces and a second one for fitting reflectivity traces.
#  Each class should have methods for running the calculations using the fresnel_calculation() function and for saving
#  the results in an "experiment workbook". The optical parameters (with many defaults) and some naming strings should
#  be provided when an object is instanced from the class.

# TODO: The workbook can be its own class with methods that define how data is stored and loaded for the dash app. The
#  idea is that this can be loaded by the app if a user wants to redo some modelling without starting all over again.
#  It should save a path to the datafile that was used for the analysis so that it has access to the data.

import fresnel_transfer_matrix as ftm
import numpy as np
import datetime
import os
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import pandas as pd
import plotly.express as px
import dash
import extract_SPR_spectra
import copy


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name='Session', directory=os.getcwd()):
        self.name = datetime.datetime.now().__str__()[0:16] + ' ' + name
        if not os.path.exists(directory + r'\\sessions'):
            os.mkdir(directory + r'\\pySPR sessions')
        self.location = directory + r'\\pySPR sessions'
        self.measurement_instances = []
        self.session_ID_generator = generate_id()

    def remove_experiment(self):
        """
        Remove a measurement from the session.
        :return:
        """
        pass


class SPRMeasurement:

    """
    An SPR measurement typically have some things in common, such as the sensor layers, measured angles, measured reflectivity, measurement
    time, etc. This information can be shared between different analysis methods for one measurement. This class serves
    as a base class for inheritance between different measurements within a session.

    NOTE: Create a new SPRMeasuerment object/instance every time before an analysis run is performed after the sensor properties has changed.
    Otherwise, previous runs will be broken. Essentially,


    """

    def __init__(self, data_path_, sensor_metal='Au', polarization=1):

        self.polarization = polarization

        # Load in the measurement data from a .csv file
        self.data_table = pd.read_csv(data_path_, delimiter=';')
        self.time = self.data_table.iloc[:, 0]
        self.angles = self.data_table.iloc[0, 1:]
        self.ydata = self.data_table.iloc[1, 1:]

        # Loading sensor's default optical properties
        self.wavelength = int(data_path_[-9:-6])
        self.sensor_metal = sensor_metal
        self.set_default_sensor_properties(self.sensor_metal)

    def set_default_sensor_properties(self, sensor_metal):

        # These default parameters should be set based on material layer and wavelength from loaded .csv file
        match sensor_metal:
            case 'Au' | 'gold' | 'Gold' | 'GOLD':
                self.layer_thicknesses = np.array([np.NaN, 2, 50, np.NaN])
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

            case 'sio2' | 'SiO2' | 'SIO2' | 'glass' | 'silica':
                # Fused silica values source: L. V. Rodríguez-de Marcos, J. I. Larruquert, J. A. Méndez, J. A. Aznárez.
                # Self-consistent optical constants of SiO2 and Ta2O5 films
                # Opt. Mater. Express 6, 3622-3637 (2016) (Numerical data kindly provided by Juan Larruquert)
                self.layer_thicknesses = np.array([np.NaN, 2, 50, 14, np.NaN])
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

            case 'Pd' | 'palladium' | 'Palladium' | 'PALLADIUM':
                self.layer_thicknesses = np.array([np.NaN, 2, 20, np.NaN])
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

            case 'Pt' | 'platinum' | 'Platinum' | 'PLATINUM':
                self.layer_thicknesses = np.array([np.NaN, 2, 20, np.NaN])
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

    def add_sensor_layer(self, thickness, n_re, n_im):
        """
        Add additional layers on top of the sensor (before bulk medium).
        :return:
        """
        # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
        self.layer_thicknesses = np.insert(self.layer_thicknesses, -1, thickness)
        self.refractive_indices = np.insert(self.refractive_indices, -1, n_re)
        self.extinction_coefficients = np.insert(self.extinction_coefficients, -1, n_im)
        print('Current layer thicknesses: ', self.layer_thicknesses)
        print('Current layer refractive indices: ', self.refractive_indices)
        print('Current layer extinction coefficients: ', self.extinction_coefficients)

    def remove_sensor_layer(self, layer_index):
        """
        Removes a layer from sensor.
        :return:
        """
        # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
        self.layer_thicknesses = np.delete(self.layer_thicknesses, layer_index)
        self.refractive_indices = np.delete(self.refractive_indices, layer_index)
        self.extinction_coefficients = np.delete(self.extinction_coefficients, layer_index)
        print('Current layer thicknesses: ', self.layer_thicknesses)
        print('Current layer refractive indices: ', self.refractive_indices)
        print('Current layer extinction coefficients: ', self.extinction_coefficients)

    def plot_reflectivity_trace(self, index=0, xdata=None, ydata=None, rlines=None):
        """

        """
        pass

    def plot_sensorgram(self, xdata=None, ydata=None, rlines=None):
        """

        """
        pass


class ModelledReflectivityTrace:
    """
    This class defines how a modelled reflectivity trace behaves. Note that a new measurement object should be used for
    each new layer added to the sensor!
    """

    def __init__(self, spr_measurement_object):
        self.polarization = spr_measurement_object.polarization
        self.wavelength = spr_measurement_object.wavelength
        self.layer_thicknesses = spr_measurement_object.layer_thicknesses
        self.refractive_indices = spr_measurement_object.refractive_indices
        self.extinction_coefficients = spr_measurement_object.extinction_coefficients
        self.data_path = spr_measurement_object.data_path

    def calculate_trace(self):
        pass


class FittedReflectivityTrace(ModelledReflectivityTrace):

    """

    """

    def __init__(self, spr_measurement_object, ydata_type='R'):
        super().__init__(spr_measurement_object)  # Initializes the same way as parent objects, to shorten code
        self.ydata_type = ydata_type

    def calculate_fit(self):
        pass


def add_experiment(experiment_handle, session_handle):
    """

    :return:
    """
    # This should not add measurement id, it should be done with "add_measurement"
    pass


def add_analysis(analysis_handle, session_handle):
    """

    :return:
    """
    # This should not add measurement id, it should be done with "add_measurement"
    pass


def load_session(filename):
    """
    Loads a previously initiated session.
    :return:
    """
    pass


def save_session(session_handle):
    """
    Save a session to a binary pickle file.
    :return:
    """
    pass


def generate_id():
    """
    Each time this function is called within an object it can return a new ID.
    :yield:
    """
    new_id = 1
    while True:
        yield new_id
        new_id += 1


def set_new_data_path():
    print('Select the measurement datafile')
    new_data_path = askopenfilename(title='Select the data file (.csv)')
    return new_data_path


if __name__ == '__main__':

    # Create initial session
    active_session = Session()

    # Prompt user for initial measurement path
    print('Select initial measurement file (.csv)')
    data_path = askopenfilename(title='Select initial measurement file (.csv)')

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


