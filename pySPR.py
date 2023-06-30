# This is the main file where the webapp is initiated and further selections are made. It should ask for a datafile and
# either load .csv files directly or run conversion of .spr2 files using 'extract_SPR_spectra.py' in a separate thread
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

# TODO: It is time to test the object persistence for an "active_object" variable and when it is added to a dictionary
#  within a session object. Will the pointer inside the session object dictionary still point to the same object in
#  memory? (I think so, but good to confirm).
# TODO: Next, it is time to populate the Modelled/FittedReflectivityTrace classes with methods and attributes for
#  performing the fresnel_calculation() function


import extract_SPR_spectra
import fresnel_transfer_matrix as ftm
import numpy as np
import datetime
import os
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import pandas as pd
import plotly.express as px
import dash
import copy


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name='Session', directory=os.getcwd()):
        self.name = datetime.datetime.now().__str__()[0:16] + ' ' + name
        if not os.path.exists(directory + r'\\pySPR sessions'):
            os.mkdir(directory + r'\\pySPR sessions')
        self.location = directory + r'\\pySPR sessions'
        self.sensor_instances = {}
        self.sensor_ID = generate_id()
        self.analysis_instances = {}
        self.analysis_ID = generate_id()

    def remove_experiment(self):
        """
        Remove a measurement from the session.
        :return:
        """
        pass


class Sensor:

    """
    An SPR measurement typically have some things in common, such as the sensor layers, measured angles,
    measured reflectivity, measurement time, etc. This information can be shared between different analysis methods for
    one measurement. This class serves as a basis for describing the current sensor, containing information about its
    layers and their optical properties.
    """

    def __init__(self, data_path_, sensor_metal='Au', polarization=1):

        # Load sensor's default optical properties
        self.polarization = polarization
        self.wavelength = int(data_path_[-9:-6])
        self.sensor_metal = sensor_metal
        self.set_default_optical_properties(self.sensor_metal)

    def set_default_optical_properties(self, sensor_metal):

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

    def add_material_layer(self, thickness, n_re, n_im, layer_index_=-1):
        """
        Add additional layers on top of the sensor (before bulk medium).
        :return:
        """
        # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
        self.layer_thicknesses = np.insert(self.layer_thicknesses, layer_index_, thickness)
        self.refractive_indices = np.insert(self.refractive_indices, layer_index_, n_re)
        self.extinction_coefficients = np.insert(self.extinction_coefficients, layer_index_, n_im)
        print('Sensor thicknesses: ', self.layer_thicknesses)
        print('Sensor refractive indices: ', self.refractive_indices)
        print('Sensor extinction coefficients: ', self.extinction_coefficients)

    def remove_material_layer(self, layer_index_):

        """
        Removes a layer from a sensor.
        :param layer_index_: int, which layer index to remove (starting from 1)
        :return:
        """

        # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
        self.layer_thicknesses = np.delete(self.layer_thicknesses, layer_index_-1, axis=0)
        self.refractive_indices = np.delete(self.refractive_indices, layer_index_-1, axis=0)
        self.extinction_coefficients = np.delete(self.extinction_coefficients, layer_index_-1, axis=0)
        print('Sensor thicknesses: ', self.layer_thicknesses)
        print('Sensor refractive indices: ', self.refractive_indices)
        print('Sensor extinction coefficients: ', self.extinction_coefficients)


class ModelledReflectivityTrace:
    """
    This class defines how a modelled reflectivity trace behaves. Note that a different sensor object should be used for
    each layer added to the sensor!
    """

    def __init__(self, sensor_object, data_path_):
        self.polarization = sensor_object.polarization
        self.wavelength = sensor_object.wavelength
        self.layer_thicknesses = sensor_object.layer_thicknesses
        self.refractive_indices = sensor_object.refractive_indices
        self.extinction_coefficients = sensor_object.extinction_coefficients
        self.data_path = data_path_

    def calculate_trace(self):
        pass

    def plot_reflectivity_trace(self, time_index=0, xdata_=None, ydata_=None, rlines_=None):
        """

        """
        pass

    def plot_sensorgram(self, time_=None, ydata_=None, rlines_=None):
        """

        """
        pass


class FittedReflectivityTrace(ModelledReflectivityTrace):

    """

    """

    def __init__(self, sensor_object, data_path_, ydata_type='R'):
        super().__init__(sensor_object, data_path_)  # Initializes the same way as parent objects, to shorten code
        self.ydata_type = ydata_type

    def calculate_fit(self):
        pass


def add_sensor(session_handle, data_path_, sensor_metal='Au', polarization=1):
    """
    Adds sensor objects to a session object.
    :return: a sensor object
    """
    id_ = next(session_handle.sensor_ID)
    sensor_object = Sensor(data_path_, sensor_metal=sensor_metal, polarization=polarization)
    session_handle.sensor_instances[id_] = sensor_object

    return sensor_object


def add_modelled_reflectivity_trace(session_handle, sensor_object, data_path_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_handle.analysis_ID)
    analysis_object = ModelledReflectivityTrace(sensor_object, data_path_)
    session_handle.analysis_instances[id_] = analysis_object

    return analysis_object


def add_fitted_reflectivity_trace(session_handle, sensor_object, data_path_, ydata_type='R'):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_handle.analysis_ID)
    analysis_object = FittedReflectivityTrace(sensor_object, data_path_, ydata_type=ydata_type)
    session_handle.analysis_instances[id_] = analysis_object

    return analysis_object


def load_session(filename):
    """
    Loads a previous session.
    :return:
    """
    # TODO: It should load a session object and load its dictionaries containing measurement_instances and
    #  analysis_instances.

    pass


def save_to_session(session_handle, object_handle):
    """
    Save a session to a binary pickle file.
    :return:
    """
    # TODO: Save a measurement or analysis object to a session
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


def load_data():
    print('Select the measurement data file (.csv)')
    data_path_ = askopenfilename(title='Select the measurement data file', filetypes=[('CSV files', '*.csv')])

    # Load in the measurement data from a .csv file
    data_table = pd.read_csv(data_path_, delimiter=';', skiprows=1, header=None)
    time_ = data_table.iloc[:, 0]
    angles_ = data_table.iloc[0, 1:]
    ydata_ = data_table.iloc[1:, 1:]
    return data_path, time_, angles_, ydata_


if __name__ == '__main__':

    # Create initial session
    active_session = Session()

    # Prompt user for initial measurement data
    data_path, time, angles, ydata = load_data()

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


