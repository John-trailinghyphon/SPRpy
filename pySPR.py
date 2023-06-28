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
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory


# PLACEHOLDER ATTRIBUTES (should come from different sources later, such as .spr2 file, .csv files or from user input in
# dashapp
session_name_PH = '2306XX SPR experiment testing'
session_directory_PH = askdirectory(title='Select save folder')
PH_data_path = askopenfilename(title='Select a measurement file')
PH_polarization = 1
PH_wavelength = 670
PH_layer_thicknesses = np.array([np.NaN, 2, 50, np.NaN])
PH_n_re = np.array([1.5202, 3.3105, 0.2238, 1.0003])
PH_n_im = np.array([0, 3.4556, 3.9259, 0])


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name=session_name_PH, directory=session_directory_PH):
        self.name = name
        self.location = directory
        self.measurement_instances = []
        self.measurement_ID_generator = generate_measurement_id()

    def remove_measurement_instance(self):
        """
        Remove a measurement from the session.
        :return:
        """
        pass


class SPRExperiment:

    """
    An SPR experiment typically have some things in common, such as the sensor layers, measured angles, measured reflectivity, measurement
    time, etc. This information can be shared between different analysis methods for one measurement. This class serves
    as a base class for inheritance between different measurements within a session.
    """

    def __init__(self, data_path=PH_data_path, layer_thicknesses=PH_layer_thicknesses, n_re=PH_n_re, n_im=PH_n_im, wavelength=PH_wavelength, polarization=PH_polarization):

        # Note that these attributes can be changed for a measurement instance when it is created, but having them here
        # should make them default from initial startup if not specified
        self.data_path = data_path
        self.polarization = polarization
        self.wavelength = wavelength
        self.layer_thicknesses = layer_thicknesses
        self.refractive_indices = n_re
        self.extinction_coefficients = n_im

    # def add_new_measurement(self):
    #     """
    #     Should add new data file path to the session and initiate a new measurement of a particular type (when called
    #     from that measurement instance)
    #     :param new_path:
    #     :return:
    #     """


class PlottedData(SPRExperiment):
    """

    """

    def __init__(self, xdata, ydata):
        super().__init__()

        instance_id = next(super().measurement_ID_generator)
        self.measurement_instances = self.measurement_instances.append('test')

        self.angles = xdata
        self.reflectivity = ydata

    def show(self):
        pass


class ModelledReflectivityTrace(SPRExperiment):
    """

    """

    def __init__(self, xdata, variable):
        super().__init__()

    def calculate_trace(self):
        pass


class FittedReflectivityTrace(SPRExperiment):

    """

    """

    def __init__(self, xdata, ydata, ydata_type):
        super().__init__()

    def calculate_fit(self):
        pass


def load_session():
    """
    Loads a previously initiated session.
    :return:
    """
    pass


def save_session():
    """
    Save a session to a binary pickle file.
    :return:
    """
    pass


def generate_measurement_id():
    """
    Each time this function is called it will return a new session measurement ID.
    :yield:
    """
    measurement_id = 0
    while True:
        yield measurement_id
        measurement_id += 1


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


