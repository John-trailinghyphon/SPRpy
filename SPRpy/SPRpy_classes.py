# This file contains all SPRpy classes, and also methods and functions that depend on the class objects

import datetime
import os
import scipy
import pickle
import copy
import bottleneck
import numpy as np
from SPRpy_functions import *
from fresnel_transfer_matrix import fresnel_calculation, fit_spectrum, model_dip_angle, measured_dip_angle, DIP_ANGLE_WEIGHT, TIR_shift_overlay

EXCLUSION_HEIGHT_SEARCH_METHODS = {'adaptive': 'Adaptive: sweep until the crossing, then refine (recommended)',
                                   'grid': 'Fixed grid of height steps (SPRpy <= 1.2)'}

EXCLUSION_BULK_RI_METHODS = {'edge_overlay': 'TIR edge overlay (recommended)',
                             'TIR_angle': 'Separate TIR angles (v1.2)'}


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, version, SPR_TIR_fitting_parameters, name='Session', directory=None, current_data_path=None):
        self.version = version
        self.name = datetime.datetime.now().__str__()[0:16].replace(':', '_') + ' ' + name
        if not directory:
            directory = os.getcwd().replace('\\', '/')
            if os.path.exists(directory + '/SPRpy sessions'):
                self.location = directory + '/SPRpy sessions' + '/{name_}'.format(name_=self.name)
            else:
                os.mkdir(directory + '/SPRpy sessions')
                self.location = directory + '/SPRpy sessions' + '/{name_}'.format(name_=self.name)
        else:
            if os.path.exists(directory):
                self.location = directory + '/{name_}'.format(name_=self.name)
            else:
                directory = os.getcwd().replace('\\', '/')
                if os.path.exists(directory + '/SPRpy sessions'):
                    self.location = directory + '/SPRpy sessions' + '/{name_}'.format(name_=self.name)
                else:
                    os.mkdir(directory + '/SPRpy sessions')
                    self.location = directory + '/SPRpy sessions' + '/{name_}'.format(name_=self.name)
                print('Warning! Custom session folder in config.toml does not exist. Sessions placed in default SPRpy location instead.')

        if not os.path.exists(self.location):
            os.mkdir(self.location)
        if not os.path.exists(self.location + '/Sensors'):
            os.mkdir(self.location + '/Sensors')
        if not os.path.exists(self.location + '/Analysis instances'):
            os.mkdir(self.location + '/Analysis instances')
        self.sensor_instances = {}  # NOTE: The sessions in this list are also updated when modified as current sensor object
        self.sensor_ID_count = 0
        self.fresnel_analysis_instances = {}
        self.fresnel_analysis_ID_count = 0
        self.exclusion_height_analysis_instances = {}
        self.exclusion_height_analysis_ID_count = 0
        self.current_data_path = current_data_path
        self.SPR_TIR_fitting_parameters = SPR_TIR_fitting_parameters
        self.log = datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Welcome to SPRpy!' \
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Start your session by defining your SPR sensor layers.'

    def add_log(self, message, save=False):
        """
        Adds a time stamped message to the session log (shown in the GUI), and prints it in the console.
        :param message: string
        :param save: also save the session file
        """
        line = datetime.datetime.now().__str__()[0:16] + ' >> ' + message
        self.log = self.log + '\n' + line
        print(line)
        if save:
            self.save_session()

    def update_name_and_location(self, new_name):
        """
        Updates the name and location of the session.
        :param new_name: string
        :return:
        """

        self.name = new_name
        # old_location = self.location
        os.rename(self.location, self.location.replace(self.location.split('/')[-1], new_name))
        self.location = self.location.replace(self.location.split('/')[-1], new_name)
        return

    def remove_sensor(self, sensor_object_id):
        """
        Remove a sensor object from the session.
        :return:
        """
        removed = self.sensor_instances.pop(sensor_object_id)
        removed_file_path = self.location + '/Sensors' + '/S{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following sensor object: S{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def remove_fresnel_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.fresnel_analysis_instances.pop(analysis_object_id)
        removed_file_path = self.location + '/Analysis instances' + '/FM{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following analysis object: FM{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def remove_exclusion_height_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.exclusion_height_analysis_instances.pop(analysis_object_id)
        removed_file_path = self.location + '/Analysis instances' + '/EH{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following analysis object: EH{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def save_all(self):
        """
        Saves all objects stored in the session, and the session file itself.
        :return: None
        """

        # Save session object
        with open(self.location + '/Session file (v{version_}).pickle'.format(version_=self.version.replace('.', '_')), 'wb') as save_file:
            pickle.dump(self, save_file)

        # Save sensor instances
        for sensor_id in self.sensor_instances:
            with open(self.location + '/Sensors' + '/S{id} {name}.pickle'.format(id=sensor_id, name=self.sensor_instances[sensor_id].name), 'wb') as save_file:
                pickle.dump(self.sensor_instances[sensor_id], save_file)

        # Save fresnel analysis instances
        for analysis_id in self.fresnel_analysis_instances:
            with open(self.location + '/Analysis instances' + '/FM{id} {name}.pickle'.format(id=analysis_id, name=self.fresnel_analysis_instances[analysis_id].name), 'wb') as save_file:
                pickle.dump(self.fresnel_analysis_instances[analysis_id], save_file)

        # Save exclusion height analysis instances
        for analysis_id in self.exclusion_height_analysis_instances:
            with open(self.location + '/Analysis instances' + '/EH{id} {name}.pickle'.format(id=analysis_id, name=self.exclusion_height_analysis_instances[analysis_id].name), 'wb') as save_file:
                pickle.dump(self.exclusion_height_analysis_instances[analysis_id], save_file)

        return

    def save_session(self):

        # Save session object
        with open(self.location + '/Session file (v{version_}).pickle'.format(version_=self.version.replace('.', '_')), 'wb') as save_file:
            pickle.dump(self, save_file)

        return

    def save_sensor(self, sensor_id):
        """
        Saves a single sensor object to the session.
        :return: None
        """

        with open(self.location + '/Sensors' + '/S{id} {name}.pickle'.format(id=sensor_id, name=self.sensor_instances[
            sensor_id].name), 'wb') as save_file:
            pickle.dump(self.sensor_instances[sensor_id], save_file)

        return

    def save_fresnel_analysis(self, analysis_id):
        """
        Saves a single fresnel analysis object to the session.
        :return: None
        """

        with open(self.location + '/Analysis instances' + '/FM{id} {name}.pickle'.format(id=analysis_id, name=self.fresnel_analysis_instances[analysis_id].name), 'wb') as save_file:
            pickle.dump(self.fresnel_analysis_instances[analysis_id], save_file)

        return

    def save_exclusion_height_analysis(self, analysis_id):
        """
        Saves a single fresnel analysis object to the session.
        :return: None
        """

        with open(self.location + '/Analysis instances' + '/EH{id} {name}.pickle'.format(id=analysis_id, name=self.exclusion_height_analysis_instances[analysis_id].name), 'wb') as save_file:
            pickle.dump(self.exclusion_height_analysis_instances[analysis_id], save_file)

        return

    def import_sensor(self):

        file_path_ = select_file(prompt='Select the sensor object', prompt_folder=self.location + '/Sensors')
        self.sensor_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            sensor_object = pickle.load(import_file)

        sensor_object.object_id = self.sensor_ID_count
        self.sensor_instances[self.sensor_ID_count] = sensor_object

        return

    def import_fresnel_analysis(self):
        file_path_ = select_file(prompt='Select the analysis object', prompt_folder=self.location + '/Analysis instances')
        self.fresnel_analysis_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            analysis_object = pickle.load(import_file)

        analysis_object.object_id = self.fresnel_analysis_ID_count
        self.fresnel_analysis_instances[analysis_object.object_id] = analysis_object

        return

    def import_exclusion_height_analysis(self):
        file_path_ = select_file(prompt='Select the analysis object', prompt_folder=self.location + '/Analysis instances')
        self.exclusion_height_analysis_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            analysis_object = pickle.load(import_file)

        analysis_object.object_id = self.exclusion_height_analysis_ID_count
        self.exclusion_height_analysis_instances[analysis_object.object_id] = analysis_object

        return

class Sensor:

    """
      An SPR measurement typically have some things in common, such as the sensor layers, measured angles,
    measured reflectivity, measurement time, etc. This information can be shared between different analysis methods for
    one measurement. This class serves as a basis for describing the current sensor, containing information about its
    layers and their optical properties.
    """

    def __init__(self, data_path_, object_id_, default_sensor_values, object_name_='Gold sensor', sensor_metal='Au', data_type='R'):

        # Load sensor's default optical properties
        self.object_id = object_id_
        self.name = object_name_
        self.data_path = data_path_
        self.data_type = data_type
        self.wavelength = int(data_path_[-9:-6])
        self.channel = data_path_[-12:-4].replace('_', ' ')
        self.sensor_metal = sensor_metal
        self.default_sensor_values = default_sensor_values
        self.set_default_optical_properties(self.sensor_metal)
        self.fitted_var = self.optical_parameters.iloc[self.fitted_layer_index]

    def set_default_optical_properties(self, sensor_metal):

        # These default parameters should be set based on material layer and wavelength from loaded .csv file

        match sensor_metal:
            case 'Au' | 'gold' | 'Gold' | 'GOLD':
                self.layer_thicknesses = np.array([np.nan, float(self.default_sensor_values["d_nm"][1]), float(self.default_sensor_values["d_nm"][2]), np.nan])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_670"][0]), float(self.default_sensor_values["n_670"][1]), float(self.default_sensor_values["n_670"][2]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_670"][0]), float(self.default_sensor_values["k_670"][1]), float(self.default_sensor_values["k_670"][2]), 0.0])
                    case 785:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_785"][0]), float(self.default_sensor_values["n_785"][1]), float(self.default_sensor_values["n_785"][2]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_785"][0]), float(self.default_sensor_values["k_785"][1]), float(self.default_sensor_values["k_785"][2]), 0.0])
                    case 850:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_850"][0]), float(self.default_sensor_values["n_850"][1]), float(self.default_sensor_values["n_850"][2]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_850"][0]), float(self.default_sensor_values["k_850"][1]), float(self.default_sensor_values["k_850"][2]), 0.0])
                    case 980:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_980"][0]), float(self.default_sensor_values["n_980"][1]), float(self.default_sensor_values["n_980"][2]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_980"][0]), float(self.default_sensor_values["k_980"][1]), float(self.default_sensor_values["k_980"][2]), 0.0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Au', 'Bulk'],  # NOTE: sensor_object.optical_parameters will auto-update when changing sensor_object.layer_thicknesses, sensor_object.refractive_indices or sensor_object.extinction_coefficients
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'sio2' | 'SiO2' | 'SIO2' | 'Glass' | 'glass' | 'silica':
                # Fused silica values source: L. V. Rodríguez-de Marcos, J. I. Larruquert, J. A. Méndez, J. A. Aznárez. Self-consistent optical constants of SiO2 and Ta2O5 films. Opt. Mater. Express 6, 3622-3637 (2016)
                self.layer_thicknesses = np.array([np.nan, float(self.default_sensor_values["d_nm"][6]), float(self.default_sensor_values["d_nm"][2]),  float(self.default_sensor_values["d_nm"][3]), np.nan])
                self.fitted_layer_index = (3, 1)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_670"][0]), float(self.default_sensor_values["n_670"][6]), float(self.default_sensor_values["n_670"][2]), float(self.default_sensor_values["n_670"][3]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_670"][0]), float(self.default_sensor_values["k_670"][6]), float(self.default_sensor_values["k_670"][2]), float(self.default_sensor_values["k_670"][3]), 0.0])
                    case 785:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_785"][0]), float(self.default_sensor_values["n_785"][6]), float(self.default_sensor_values["n_785"][2]), float(self.default_sensor_values["n_785"][3]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_785"][0]), float(self.default_sensor_values["k_785"][6]), float(self.default_sensor_values["k_785"][2]), float(self.default_sensor_values["k_785"][3]), 0.0])
                    case 850:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_850"][0]), float(self.default_sensor_values["n_850"][6]), float(self.default_sensor_values["n_850"][2]), float(self.default_sensor_values["n_850"][3]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_850"][0]), float(self.default_sensor_values["k_850"][6]), float(self.default_sensor_values["k_850"][2]), float(self.default_sensor_values["k_850"][3]), 0.0])
                    case 980:
                        self.refractive_indices = np.array([float(self.default_sensor_values["n_980"][0]), float(self.default_sensor_values["n_980"][6]), float(self.default_sensor_values["n_980"][2]), float(self.default_sensor_values["n_980"][3]), 1.0003])
                        self.extinction_coefficients = np.array([float(self.default_sensor_values["k_980"][0]), float(self.default_sensor_values["k_980"][6]), float(self.default_sensor_values["k_980"][2]), float(self.default_sensor_values["k_980"][3]), 0.0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Ti', 'Au', 'SiO2', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'Pd' | 'palladium' | 'Palladium' | 'PALLADIUM':
                # Source Pd 670nm : Andersson, John, et al. "Surface plasmon resonance sensing with thin films of palladium and platinum–quantitative and real-time analysis." Physical Chemistry Chemical Physics 24.7 (2022): 4588-4594.
                # Source Pd 785nm, 980nm: Similar procedure as Andersson, John et al., but these were not reported.
                self.layer_thicknesses = np.array([np.nan, float(self.default_sensor_values["d_nm"][1]), float(self.default_sensor_values["d_nm"][4]), np.nan])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_670"][0]), float(self.default_sensor_values["n_670"][1]), float(self.default_sensor_values["n_670"][4]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_670"][0]), float(self.default_sensor_values["k_670"][1]), float(self.default_sensor_values["k_670"][4]), 0.0])
                    case 785:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_785"][0]), float(self.default_sensor_values["n_785"][1]), float(self.default_sensor_values["n_785"][4]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_785"][0]), float(self.default_sensor_values["k_785"][1]), float(self.default_sensor_values["k_785"][4]), 0.0])
                    case 850:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_850"][0]), float(self.default_sensor_values["n_850"][1]), float(self.default_sensor_values["n_850"][4]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_850"][0]), float(self.default_sensor_values["k_850"][1]), float(self.default_sensor_values["k_850"][4]), 0.0])
                    case 980:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_980"][0]), float(self.default_sensor_values["n_980"][1]), float(self.default_sensor_values["n_980"][4]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_980"][0]), float(self.default_sensor_values["k_980"][1]), float(self.default_sensor_values["k_980"][4]), 0.0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Pd', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'Pt' | 'platinum' | 'Platinum' | 'PLATINUM':
                # Source Pt 670 nm: Andersson, John, et al. "Surface plasmon resonance sensing with thin films of palladium and platinum–quantitative and real-time analysis." Physical Chemistry Chemical Physics 24.7 (2022): 4588-4594.
                # Source Pt 785, 850 and 980 nm: W. S. M. Werner, K. Glantschnig, C. Ambrosch-Draxl. Optical constants and inelastic electron-scattering data for 17 elemental metals, J. Phys Chem Ref. Data 38, 1013-1092 (2009)
                self.layer_thicknesses = np.array([np.nan, 2.00, 20.00, np.nan])
                self.fitted_layer_index = (2, 3)  # Tuple with index for df.iloc[fitted_layer_index]
                match self.wavelength:
                    case 670:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_670"][0]), float(self.default_sensor_values["n_670"][1]), float(self.default_sensor_values["n_670"][5]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_670"][0]), float(self.default_sensor_values["k_670"][1]), float(self.default_sensor_values["k_670"][5]), 0.0])
                    case 785:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_785"][0]), float(self.default_sensor_values["n_785"][1]), float(self.default_sensor_values["n_785"][5]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_785"][0]), float(self.default_sensor_values["k_785"][1]), float(self.default_sensor_values["k_785"][5]), 0.0])
                    case 850:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_850"][0]), float(self.default_sensor_values["n_850"][1]), float(self.default_sensor_values["n_850"][5]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_850"][0]), float(self.default_sensor_values["k_850"][1]), float(self.default_sensor_values["k_850"][5]), 0.0])
                    case 980:
                        self.refractive_indices = np.array(
                            [float(self.default_sensor_values["n_980"][0]), float(self.default_sensor_values["n_980"][1]), float(self.default_sensor_values["n_980"][5]), 1.0003])
                        self.extinction_coefficients = np.array(
                            [float(self.default_sensor_values["k_980"][0]), float(self.default_sensor_values["k_980"][1]), float(self.default_sensor_values["k_980"][5]), 0.0])
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Pt', 'Bulk'],
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})
        return

class FresnelModel:
    """
    This class defines how a modelled reflectivity trace behaves.

    TODO: Each object should also have a .csv export function.

    """

    def __init__(self, session_object, sensor_object_, data_path_, reflectivity_df_, object_id_, object_name_):
        self.name = object_name_
        self.object_id = object_id_
        self.sensor_object = sensor_object_
        self.sensor_object_label = ''
        self.initial_data_path = data_path_
        self.measurement_data = reflectivity_df_
        self.SPR_TIR_fitting_parameters = session_object.SPR_TIR_fitting_parameters  # We want the SPR and TIR fitting parameters to update when they change globally in the session
        self.angle_range = session_object.SPR_TIR_fitting_parameters['Fresnel_angle_range_points']
        self.polarization = 1.0  # 0-1, for degree of s (=0) and p(=1) polarization,
        self.ini_guess = np.array(4)
        self.bounds = [0, 50]  # or [(lb1, lb2), (ub1, ub2)] etc
        self.extinction_correction = 0
        self.y_offset = 0
        self.fit_offset = True
        self.fit_prism_k = True
        self.fitted_data = None
        self.fitted_result = np.array([0, 0, 0])
        self.update_bulk_RI_from_TIR = True  # If False, the bulk refractive index in the sensor table is used as is
        self.match_dip_angle = True  # Fit the model dip at the measured SPR angle (see fresnel_residuals); False for older sessions
        self.measured_dip_angle = np.nan  # SPR angle of the measured spectrum (deg)
        self.dip_angle_mismatch = np.nan  # Angle of the fitted model dip minus the measured SPR angle (deg)
        self.fitted_layer_index = copy.deepcopy(self.sensor_object.fitted_layer_index)  # Deepcopy to not overwrite if fitting multiple layers in the same sensor object
        self.fitted_layer = copy.deepcopy(self.sensor_object.optical_parameters.iloc[self.fitted_layer_index[0], 0])

    def calculate_fresnel_trace(self):

        fresnel_coefficients_ = fresnel_calculation(None,
                                                    angles=self.angle_range,
                                                    fitted_layer_index=self.sensor_object.fitted_layer_index,
                                                    wavelength=self.sensor_object.wavelength,
                                                    layer_thicknesses=self.sensor_object.layer_thicknesses,
                                                    n_re=self.sensor_object.refractive_indices,
                                                    n_im=self.sensor_object.extinction_coefficients,
                                                    ydata=None,
                                                    ydata_type=self.sensor_object.data_type,
                                                    polarization=self.polarization)
        return fresnel_coefficients_

    def model_reflectivity_trace(self):

        xdata_ = self.measurement_data['angles']
        ydata_ = self.measurement_data['ydata']

        # Calculate TIR angle and bulk refractive index (unless the bulk refractive index is set manually)
        if getattr(self, 'update_bulk_RI_from_TIR', True):
            TIR_angle, _, _, _, _ = TIR_determination(xdata_, ydata_, self.SPR_TIR_fitting_parameters)
            self.sensor_object.refractive_indices[-1] = self.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)

        # Add extinction correction to fitted surface layer extinction value
        self.sensor_object.extinction_coefficients[0] += self.extinction_correction  # Manually correct prism layer extinction

        # Selecting a range of measurement data to use for fitting, and including an offset in reflectivity (iterated 3 times)
        selection_xdata_ = xdata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])]
        selection_ydata_ = ydata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])]

        # Weighing options
        weights = None

        # weights_ = np.abs(np.diff(selection_ydata_))+1  # Highest derivative
        # weights = np.append(weights_, 0)

        # weights = 1/selection_ydata_

        # weights = np.ones(len(selection_ydata_))
        # weights[selection_ydata_.argmin():-1] = 2

        # SPR angle of the measured spectrum: the model dip is placed there with match_dip_angle, and the difference
        # between the model dip and this angle is reported in either case
        match_dip_angle = getattr(self, 'match_dip_angle', False)
        self.measured_dip_angle = measured_dip_angle(xdata_, ydata_, self.SPR_TIR_fitting_parameters)

        # Perform the first fitting
        result = fit_spectrum(self.ini_guess, self.bounds,
                              {'fitted_layer_index': self.sensor_object.fitted_layer_index,
                               'wavelength': self.sensor_object.wavelength,
                               'layer_thicknesses': self.sensor_object.layer_thicknesses,
                               'n_re': self.sensor_object.refractive_indices,
                               'n_im': self.sensor_object.extinction_coefficients,
                               'angles': selection_xdata_,
                               'ydata': selection_ydata_,
                               'weights': weights,
                               'ydata_type': self.sensor_object.data_type,
                               'polarization': self.polarization,
                               'ydata_offset': self.y_offset},
                              measured_dip_angle=self.measured_dip_angle,
                              dip_weight=DIP_ANGLE_WEIGHT if match_dip_angle else 0.0,
                              loss='huber', ftol=1e-12, xtol=1e-12, gtol=1e-12)

        # Collect the results from least_squares object and calculate corresponding fresnel coefficients
        self.fitted_result = np.array(result['x'])

        if self.fit_offset:
            self.y_offset = self.fitted_result[1]
        else:
            self.y_offset = 0

        fresnel_coefficients = fresnel_calculation(self.fitted_result,
                                                   fitted_layer_index=self.sensor_object.fitted_layer_index,
                                                   angles=selection_xdata_,
                                                   wavelength=self.sensor_object.wavelength,
                                                   layer_thicknesses=self.sensor_object.layer_thicknesses,
                                                   n_re=self.sensor_object.refractive_indices,
                                                   n_im=self.sensor_object.extinction_coefficients,
                                                   ydata=None,
                                                   weights=weights,
                                                   ydata_type='R',
                                                   polarization=self.polarization,
                                                   ydata_offset=self.y_offset
                                                   )

        # Compile into fresnel_coefficients data frame
        self.fitted_data = pd.DataFrame(data={'angles': selection_xdata_, 'ydata': fresnel_coefficients})

        # Angle of the fitted model dip relative to the measured SPR angle
        measured = self.measured_dip_angle if np.isfinite(self.measured_dip_angle) else float(selection_xdata_.iloc[int(np.argmin(selection_ydata_.to_numpy()))])
        self.dip_angle_mismatch = model_dip_angle(self.fitted_result, measured,
                                                  fitted_layer_index=self.sensor_object.fitted_layer_index, wavelength=self.sensor_object.wavelength,
                                                  layer_thicknesses=self.sensor_object.layer_thicknesses, n_re=self.sensor_object.refractive_indices,
                                                  n_im=self.sensor_object.extinction_coefficients, polarization=self.polarization,
                                                  ydata_offset=self.y_offset) - measured if np.isfinite(self.measured_dip_angle) else np.nan

        return self.fitted_data

    def export_fitted_results(self):
        """
        Exporting the result (including parameters) of a particular analysis as a .csv file
        :return:
        """
        pass

    def export_calculated_results(self):
        """
        Exporting the calculated fresnel traces (including parameters) of a particular analysis as a .csv file
        :return:
        """
        pass


class ExclusionHeight:

    """
        This class defines an analysis object for determining the exclusion height from measurement data with probe
        injections. The underlying method is described as the "non-interacting probe method" in the literature.
    """

    def __init__(self, session_object, fresnel_object_, sensorgram_df_, data_path_,  object_id_, object_name_):
        self.name = object_name_
        self.object_id = object_id_
        self.fresnel_object = fresnel_object_
        self.SPR_TIR_fitting_parameters = session_object.SPR_TIR_fitting_parameters  # We want the SPR and TIR fitting parameters to update when they change globally in the session
        self.fresnel_object_label = 'Fresnel background: FM{analysis_number} {analysis_name}'.format(
                                                                    analysis_number=fresnel_object_.object_id,
                                                                    analysis_name=fresnel_object_.name)
        self.sensor_object = fresnel_object_.sensor_object
        self.polarization = fresnel_object_.polarization
        self.initial_data_path = data_path_
        self.sensorgram_data = sensorgram_df_.iloc[:,:3]
        self.sensorgram_offset_ind = 0
        self.d_n_pair_resolution = 200
        self.height_bounds = [0, 200]
        self.RI_bounds = [1.2, 1.6]
        self.RI_initial_guess = 1.35
        self.height_steps = np.linspace(self.height_bounds[0], self.height_bounds[1], self.d_n_pair_resolution)
        self.points_below_SPR_min_ind = None  # This will be calculated based on the most updated fits of the Fresnel background
        self.points_above_SPR_min_ind = None  # This will be calculated based on the most updated fits of the Fresnel background
        self.injection_points = []
        self.buffer_points = []
        self.probe_points = []
        self.SPR_vs_TIR_dfs = []  # List of dataframes with labels 'SPR angles' and 'TIR angles' for indexing each step result
        self.buffer_reflectivity_dfs = []  # Use labels 'buffer reflectivity' and 'buffer angles' (and likewise probe) for indexing
        self.buffer_bulk_RIs = []  # Calculated from TIR angle of each reflectivity DF
        self.probe_reflectivity_dfs = []  # Use labels 'buffer reflectivity' and 'buffer angles' (and likewise probe) for indexing
        self.probe_bulk_RIs = []  # Calculated from TIR angle of each reflectivity DF
        self.d_n_pair_dfs = []  # Use labels 'height' and 'buffer RI' and 'probe RI' for indexing
        self.all_exclusion_results = np.array([[np.nan, np.nan],[np.nan, np.nan]])
        self.mean_exclusion_height_result = [np.nan]  # Tuple of mean value of exclusion height from all injection steps, and standard deviation
        self.mean_exclusion_RI_result = [np.nan]  # Tuple of mean value of exclusion RI from all injection steps, and standard deviation
        self.abort_flag = False
        self.fit_offset = False
        self.fit_prism = False
        self.bulk_RI_method = 'edge_overlay'  # How the bulk RI change between buffer and probe is determined, see EXCLUSION_BULK_RI_METHODS
        self.height_search = 'adaptive'  # How the crossing of the buffer and probe curves is searched, see EXCLUSION_HEIGHT_SEARCH_METHODS
        self.match_dip_angle = True  # Fit the model dip at the measured SPR angle of each buffer and probe spectrum (False for older sessions)
        self.TIR_shifts = []  # TIR shift (deg) between each buffer and its probe, used for the bulk RI change
        self.height_sensitivity = np.array([np.nan])  # Change of each exclusion height (nm) per +1 mdeg error in the TIR shift

    def initialize_model(self, ydata_df, bulk_RI_method=None):
        """
        Calculates the mean buffer and probe spectra, their bulk refractive indices and angle ranges for the fitting.

        :param ydata_df: reflectivity data frame of the measurement
        :param bulk_RI_method: 'edge_overlay' or 'TIR_angle' (see EXCLUSION_BULK_RI_METHODS). If None, the method stored
                               in the analysis object is used (sessions from SPRpy <= 1.2 use 'TIR_angle').
        """

        if bulk_RI_method is not None:
            self.bulk_RI_method = bulk_RI_method
        elif not hasattr(self, 'bulk_RI_method'):
            self.bulk_RI_method = 'TIR_angle'

        # Calculate number of points above and below minimum point based on fresnel model background range
        background_reflectivity = self.fresnel_object.measurement_data['ydata']
        background_angles = self.fresnel_object.measurement_data['angles']
        selection_criterion = (background_angles >= self.fresnel_object.angle_range[0]) & (background_angles <= self.fresnel_object.angle_range[1])
        selection_ydata_series = background_reflectivity[selection_criterion]
        smoothened_selection = bottleneck.move_mean(selection_ydata_series.to_numpy(), window=4, min_count=1)  # Ensures closer fit to minimum position
        smoothened_selection_series = pd.Series(smoothened_selection)
        self.points_below_SPR_min_ind = len(smoothened_selection_series[(smoothened_selection_series.index < smoothened_selection_series.idxmin())])
        self.points_above_SPR_min_ind = len(smoothened_selection_series[(smoothened_selection_series.index > smoothened_selection_series.idxmin())])

        # Overwrite previous probe and buffer reflectivity data frames and bulk RI lists and SPR vs TIR data frames
        self.buffer_reflectivity_dfs = []
        self.buffer_bulk_RIs = []
        self.probe_reflectivity_dfs = []
        self.probe_bulk_RIs = []
        self.SPR_vs_TIR_dfs = []

        # Mean reflectivity spectrum and TIR angle of every buffer and probe interval (pairs of selected points)
        background_angles_array = background_angles.to_numpy()

        def mean_spectra_and_TIR_angles(points, point_type):
            spectra, TIR_angles = [], []
            for point_index in range(0, 2 * (len(points) // 2), 2):
                mean_reflectivity = ydata_df.iloc[slice(*point_pair_indices(points, point_index)), :].mean(axis=0).squeeze()
                try:
                    TIR_angle, _, _, _, _ = TIR_determination(background_angles_array, mean_reflectivity.to_numpy(), self.SPR_TIR_fitting_parameters)
                except (TypeError, ValueError, IndexError):
                    raise TypeError('Something went wrong when selecting {kind} points. Please clear and reselect them and avoid clicking other point markers.'.format(kind=point_type))
                spectra.append(mean_reflectivity)
                TIR_angles.append(TIR_angle)
            return spectra, TIR_angles

        buffer_spectra, buffer_TIR_angles = mean_spectra_and_TIR_angles(self.buffer_points, 'buffer')
        probe_spectra, probe_TIR_angles = mean_spectra_and_TIR_angles(self.probe_points, 'probe')

        # Each probe has two buffer intervals (before and after), buffer interval k belongs to probe k // 2. The
        # exclusion height depends on the bulk RI change between them (not on the absolute bulk RI), i.e. on the TIR
        # shift. With 'edge_overlay' the shift is found by overlaying the TIR edges of buffer and probe, and the buffer
        # TIR angle is set to the probe TIR angle minus this shift.
        self.TIR_shifts = []
        for buffer_index in range(len(buffer_spectra)):
            probe_index = buffer_index // 2
            if probe_index >= len(probe_spectra):
                self.TIR_shifts.append(np.nan)
                continue
            if self.bulk_RI_method == 'edge_overlay':
                shift = TIR_shift_overlay(background_angles_array, buffer_spectra[buffer_index].to_numpy(), probe_spectra[probe_index].to_numpy(),
                                          buffer_TIR_angles[buffer_index], probe_TIR_angles[probe_index])
                buffer_TIR_angles[buffer_index] = probe_TIR_angles[probe_index] - shift
            else:
                shift = probe_TIR_angles[probe_index] - buffer_TIR_angles[buffer_index]
            self.TIR_shifts.append(shift)

        prism_RI = self.sensor_object.refractive_indices[0]
        self.buffer_bulk_RIs = [prism_RI * np.sin(np.pi / 180 * TIR_angle) for TIR_angle in buffer_TIR_angles]
        self.probe_bulk_RIs = [prism_RI * np.sin(np.pi / 180 * TIR_angle) for TIR_angle in probe_TIR_angles]

        # Angle range around the SPR minimum of each mean spectrum, with the same number of points as the Fresnel background
        def fitting_range(mean_reflectivity):
            minimum_ind = pd.Series(bottleneck.move_mean(mean_reflectivity.to_numpy(), window=4, min_count=1)).idxmin()
            selection = slice(minimum_ind - self.points_below_SPR_min_ind, minimum_ind + self.points_above_SPR_min_ind)
            return pd.DataFrame(data={'reflectivity': mean_reflectivity.iloc[selection], 'angles': background_angles.iloc[selection]})

        self.buffer_reflectivity_dfs = [fitting_range(spectrum) for spectrum in buffer_spectra]
        self.probe_reflectivity_dfs = [fitting_range(spectrum) for spectrum in probe_spectra]

        # Create SPR vs TIR data frames
        injectionpoint_index = 0
        offset_SPR_sensorgram = self.sensorgram_data['SPR angle'] - self.sensorgram_data.loc[self.sensorgram_offset_ind, 'SPR angle']
        offset_TIR_sensorgram = self.sensorgram_data['TIR angle'] - self.sensorgram_data.loc[self.sensorgram_offset_ind, 'TIR angle']

        for reflectivity_index in range(int(len(self.injection_points) / 2)):

            self.SPR_vs_TIR_dfs.append(pd.DataFrame(data={'SPR angles': offset_SPR_sensorgram.iloc[slice(*point_pair_indices(self.injection_points, injectionpoint_index))],
                                                          'TIR angles': offset_TIR_sensorgram.iloc[slice(*point_pair_indices(self.injection_points, injectionpoint_index))],
                                                          })
                                       )

            # Next pair of injection point indices
            injectionpoint_index += 2

        return


def point_pair_indices(points, first_index):
    """
    Start and end index of a pair of selected sensorgram points (points[first_index] and points[first_index + 1]),
    sorted so that the two points of a pair can be selected in any order.
    """
    start, end = sorted((int(points[first_index][0]), int(points[first_index + 1][0])))
    return start, end


def exclusion_fit_settings(exclusion_height_analysis_object):
    """
    Initial guess and bounds for fitting the surface layer RI (and optionally the intensity offset and prism k).
    """
    eh = exclusion_height_analysis_object
    if not eh.fit_offset:
        return [eh.RI_initial_guess], (eh.RI_bounds[0], eh.RI_bounds[1])
    offset_guess = float(np.ravel(eh.fresnel_object.ini_guess)[1]) if np.size(eh.fresnel_object.ini_guess) >= 2 else eh.fresnel_object.y_offset
    if not eh.fit_prism:
        return [eh.RI_initial_guess, offset_guess], [(eh.RI_bounds[0], -np.inf), (eh.RI_bounds[1], np.inf)]
    prism_k_guess = float(np.ravel(eh.fresnel_object.ini_guess)[2]) if np.size(eh.fresnel_object.ini_guess) >= 3 else 0.01
    return [eh.RI_initial_guess, offset_guess, prism_k_guess], [(eh.RI_bounds[0], -np.inf, 0), (eh.RI_bounds[1], np.inf, 0.1)]


def fit_surface_layer(exclusion_height_analysis_object, reflectivity_df, bulk_RI, height, initial_guess=None, bounds=None):
    """
    Fits the RI of the surface layer (and optionally the offset and prism k) at a given surface layer height and bulk RI.
    The sensor object of the analysis is not changed.

    :return: fitted values (RI, [offset, [prism k]])
    """
    eh = exclusion_height_analysis_object
    if initial_guess is None:
        initial_guess, bounds = exclusion_fit_settings(eh)

    refractive_indices = np.array(eh.sensor_object.refractive_indices, dtype=np.float64)
    refractive_indices[-1] = bulk_RI
    layer_thicknesses = np.array(eh.sensor_object.layer_thicknesses, dtype=np.float64)
    layer_thicknesses[-2] = height  # Surface layer height
    extinction_coefficients = np.array(eh.sensor_object.extinction_coefficients, dtype=np.float64)

    # SPR angle of the measured spectrum, at which the model dip is placed (with match_dip_angle)
    match_dip_angle = getattr(eh, 'match_dip_angle', False)
    dip_angle = measured_dip_angle(reflectivity_df['angles'], reflectivity_df['reflectivity'], eh.SPR_TIR_fitting_parameters) if match_dip_angle else None

    result = fit_spectrum(initial_guess, bounds,
                          {'fitted_layer_index': (-2, 2),  # Should always be the RI of the surface layer
                           'wavelength': eh.sensor_object.wavelength,
                           'layer_thicknesses': layer_thicknesses,
                           'n_re': refractive_indices,
                           'n_im': extinction_coefficients,
                           'angles': reflectivity_df['angles'].to_numpy(),
                           'ydata': reflectivity_df['reflectivity'].to_numpy(),
                           'ydata_type': eh.sensor_object.data_type,
                           'ydata_offset': eh.fresnel_object.y_offset,
                           'polarization': eh.polarization},
                          measured_dip_angle=dip_angle, dip_weight=DIP_ANGLE_WEIGHT if match_dip_angle else 0.0)
    return result['x']


def calculate_exclusion_height(exclusion_height_analysis_object_copy, buffer_or_probe_flag, data_frame_index):
    """
    This function calculates the exclusion height for a single injection step

    :param exclusion_height_analysis_object_copy: object containing all parameters and data
    :param buffer_or_probe_flag: either 'buffer' or 'probe'
    :param data_frame_index: index of dataframe and buffer or probe RI bulk value
    :return: [RI results, offset results, prism k results] for each height step
    """

    eh = exclusion_height_analysis_object_copy
    if buffer_or_probe_flag == 'buffer':
        reflectivity_df = eh.buffer_reflectivity_dfs[data_frame_index]
        bulk_RI = eh.buffer_bulk_RIs[data_frame_index]
    elif buffer_or_probe_flag == 'probe':
        reflectivity_df = eh.probe_reflectivity_dfs[data_frame_index]
        bulk_RI = eh.probe_bulk_RIs[data_frame_index]
    else:
        raise ValueError('Only buffer or probe allowed')

    initial_guess, bounds = exclusion_fit_settings(eh)
    RI_results, offset_results, prism_k_results = [], [], []
    for height in eh.height_steps:
        fitted = fit_surface_layer(eh, reflectivity_df, bulk_RI, height, initial_guess, bounds)
        RI_results.append(fitted[0])
        if eh.fit_offset:
            offset_results.append(fitted[1])
            if eh.fit_prism:
                prism_k_results.append(fitted[2])

    return [RI_results, offset_results, prism_k_results]


def exclusion_height_crossing(height_steps, buffer_RI_result, probe_RI_result):
    """
    Exclusion height and RI: the first height where the fitted surface layer RI of the buffer curve drops below that of
    the probe curve (located within the height step with a 200 point linear interpolation).

    :return: (height, RI), or (nan, nan) if the curves do not cross within the height steps
    """
    for height_ind in range(1, len(height_steps)):
        if buffer_RI_result[height_ind] < probe_RI_result[height_ind]:
            buffer_RI_zoom_range = np.linspace(buffer_RI_result[height_ind - 1], buffer_RI_result[height_ind], 200)
            probe_RI_zoom_range = np.linspace(probe_RI_result[height_ind - 1], probe_RI_result[height_ind], 200)
            height_zoom_range = np.linspace(height_steps[height_ind - 1], height_steps[height_ind], 200)
            index = int(np.argmax(buffer_RI_zoom_range < probe_RI_zoom_range))
            return float(height_zoom_range[index]), float(buffer_RI_zoom_range[index])
    return np.nan, np.nan


def exclusion_height_sensitivity(exclusion_height_analysis_object, probe_index, exclusion_height, buffer_RI_result=None, probe_RI_result=None, slopes=None):
    """
    Change of the exclusion height (nm) if the TIR shift between buffer and probe (i.e. the bulk RI change) is 1 mdeg
    larger. The exclusion height is where the buffer and probe curves of surface layer RI vs height cross, so it moves
    by dh = (dRI_probe/dn_bulk) / (dRI_buffer/dh - dRI_probe/dh) * dn_bulk.
    """
    eh = exclusion_height_analysis_object
    heights = np.asarray(eh.height_steps, dtype=np.float64)
    if not np.isfinite(exclusion_height) or (slopes is None and len(heights) < 3):
        return np.nan

    prism_RI = eh.sensor_object.refractive_indices[0]
    probe_bulk_RI = eh.probe_bulk_RIs[probe_index]
    TIR_angle = np.arcsin(probe_bulk_RI / prism_RI)
    delta_bulk_RI = prism_RI * np.cos(TIR_angle) * np.radians(1e-3)  # Bulk RI change for 1 mdeg TIR shift

    probe_df = eh.probe_reflectivity_dfs[probe_index]
    RI_plus = fit_surface_layer(eh, probe_df, probe_bulk_RI + delta_bulk_RI, exclusion_height)[0]
    RI_minus = fit_surface_layer(eh, probe_df, probe_bulk_RI - delta_bulk_RI, exclusion_height)[0]
    probe_RI_per_bulk_RI = (RI_plus - RI_minus) / (2 * delta_bulk_RI)

    if slopes is not None:
        buffer_slope, probe_slope = slopes
    else:
        buffer_slope = np.interp(exclusion_height, heights, np.gradient(np.asarray(buffer_RI_result, dtype=np.float64), heights))
        probe_slope = np.interp(exclusion_height, heights, np.gradient(np.asarray(probe_RI_result, dtype=np.float64), heights))
    if buffer_slope == probe_slope:
        return np.nan
    return float(probe_RI_per_bulk_RI / (buffer_slope - probe_slope) * delta_bulk_RI)


def adaptive_exclusion_height(exclusion_height_analysis_object, buffer_index, probe_index, tolerance=1e-3, slope_step=0.25,
                              minimum_height=0.5, coarse_steps=31):
    """
    Exclusion height of one injection step with a two-pass search (as in SPR_WetHeight.m by J. Andersson and
    J. Medin):
     1) Coarse pass: coarse_steps heights between the height bounds (the first and last of height_steps) are swept until
        the buffer and probe curves of surface layer RI vs height cross (the buffer RI drops below the probe RI).
     2) Fine pass: the curves are calculated with len(height_steps) points (the resolution) in a window around the
        coarse crossing (the coarse step containing it, plus half a coarse step on each side), and the crossing is
        located within the tolerance with a root finder (Brent's method) on the fitted curves.
    Heights where a fitted RI lies at the RI bounds, and heights below minimum_height (nm), are skipped since the curves
    are not defined there. If the curves have already crossed at the first usable height, the heights before it are
    swept with finer steps. A wide height range (e.g. 0-300 nm for thick brushes) therefore costs little extra time, and
    the curves are well resolved around the crossing.

    :return: dict with 'heights', 'buffer' and 'probe' (fitted values at each height of the fine pass, or of the coarse
             pass if no crossing was found, sorted by height), 'height' and 'RI' (crossing, nan if none) and
             'sensitivity' (nm per mdeg TIR shift)
    """
    eh = exclusion_height_analysis_object
    initial_guess, bounds = exclusion_fit_settings(eh)
    buffer_df, buffer_bulk_RI = eh.buffer_reflectivity_dfs[buffer_index], eh.buffer_bulk_RIs[buffer_index]
    probe_df, probe_bulk_RI = eh.probe_reflectivity_dfs[probe_index], eh.probe_bulk_RIs[probe_index]
    RI_low, RI_high = eh.RI_bounds[0], eh.RI_bounds[1]
    fits = {}

    def fitted(height):
        height = float(height)
        if height not in fits:
            fits[height] = (fit_surface_layer(eh, buffer_df, buffer_bulk_RI, height, initial_guess, bounds),
                            fit_surface_layer(eh, probe_df, probe_bulk_RI, height, initial_guess, bounds))
        return fits[height]

    def difference(height):
        buffer_fit, probe_fit = fitted(height)
        if min(buffer_fit[0], probe_fit[0]) <= RI_low + 1e-6 or max(buffer_fit[0], probe_fit[0]) >= RI_high - 1e-6:
            return np.nan  # At the RI bounds
        return buffer_fit[0] - probe_fit[0]

    def sweep(heights):
        """
        First bracket (h1, d1, h2, d2) where the difference changes from positive to negative, or the gap before the
        first usable height if the curves have already crossed there (h1, nan, h2, d2), or None.
        """
        previous = None
        for index, height in enumerate(heights):
            value = difference(height)
            if not np.isfinite(value):
                continue
            if previous is None and value <= 0:
                return (heights[index - 1], np.nan, height, value) if index > 0 else None
            if previous is not None and previous[1] > 0 and value <= 0:
                return previous[0], previous[1], height, value
            previous = (height, value)
        return None

    def find_bracket(heights):
        for _ in range(4):  # The sweep is refined if the curves have already crossed at the first usable height
            bracket = sweep(heights)
            if bracket is None or not np.isnan(bracket[1]):
                return bracket
            if bracket[2] - bracket[0] < 10 * tolerance:
                return None
            heights = np.linspace(bracket[0], bracket[2], 11)
        return None

    height_steps = np.asarray(eh.height_steps, dtype=np.float64)
    lower, upper = max(float(height_steps[0]), minimum_height), max(float(height_steps[-1]), minimum_height)
    fine_points = max(len(height_steps), 3)

    # 1) Coarse pass
    coarse_heights = np.unique(np.linspace(lower, upper, max(2, min(coarse_steps, fine_points))))
    bracket = find_bracket(coarse_heights)
    crossing_height, crossing_RI, sensitivity = np.nan, np.nan, np.nan
    curve_heights = coarse_heights

    # 2) Fine pass around the coarse crossing
    if bracket is not None:
        step = bracket[2] - bracket[0]
        window = (max(lower, bracket[0] - step / 2), min(upper, bracket[2] + step / 2))
        curve_heights = np.linspace(window[0], window[1], fine_points)
        fine_bracket = find_bracket(curve_heights)
        if fine_bracket is None:
            fine_bracket = bracket  # Not resolved on the fine grid (e.g. fits at the RI bounds), use the coarse step
        low, low_value, high, high_value = fine_bracket
        try:
            crossing_height = scipy.optimize.brentq(lambda h: np.nan_to_num(difference(h), nan=0.0), low, high, xtol=tolerance)
        except ValueError:
            crossing_height = low + (high - low) * low_value / (low_value - high_value) if np.isfinite(low_value) else high
        crossing_RI = float(fitted(crossing_height)[0][0])
        curve_heights = np.unique(np.concatenate([curve_heights, [crossing_height]]))

    if np.isfinite(crossing_height):
        try:
            low, high = max(crossing_height - slope_step, lower), crossing_height + slope_step
            slopes = ((fitted(high)[0][0] - fitted(low)[0][0]) / (high - low), (fitted(high)[1][0] - fitted(low)[1][0]) / (high - low))
            sensitivity = exclusion_height_sensitivity(eh, probe_index, crossing_height, slopes=slopes)
        except Exception as error:
            print('Could not estimate the sensitivity of the exclusion height: {error}'.format(error=error))

    heights = np.array(sorted(float(h) for h in curve_heights))
    return {'heights': heights,
            'buffer': np.array([fitted(h)[0] for h in heights]),
            'probe': np.array([fitted(h)[1] for h in heights]),
            'height': float(crossing_height), 'RI': float(crossing_RI), 'sensitivity': sensitivity}


def store_exclusion_height_step(exclusion_height_analysis_object, step, heights, buffer_fits, probe_fits, height, RI, sensitivity):
    """
    Stores the curves (fitted values vs height) and the result of one injection step in the analysis object.
    """
    eh = exclusion_height_analysis_object
    buffer_fits, probe_fits = np.atleast_2d(buffer_fits), np.atleast_2d(probe_fits)
    data = {'height': heights, 'buffer RI': buffer_fits[:, 0], 'probe RI': probe_fits[:, 0]}
    if eh.fit_offset and buffer_fits.shape[1] > 1:
        data.update({'buffer offsets': buffer_fits[:, 1], 'probe offsets': probe_fits[:, 1]})
        if eh.fit_prism and buffer_fits.shape[1] > 2:
            data.update({'buffer prism k': buffer_fits[:, 2], 'probe prism k': probe_fits[:, 2]})
    eh.d_n_pair_dfs[step] = pd.DataFrame(data=data)
    eh.all_exclusion_results[0, step] = height
    eh.all_exclusion_results[1, step] = RI
    eh.height_sensitivity[step] = sensitivity


def run_exclusion_height_calculation(exclusion_height_analysis_object, logical_cores=None):
    """
    Runs the exclusion height calculation of an initialized analysis object (height steps and settings set) for all
    injection steps (each buffer interval with its probe interval), and calculates the mean exclusion height and RI with
    standard deviations. The calculations run in the current process: with the vectorised fresnel model a full
    calculation takes seconds, less than starting parallel processes would (logical_cores is no longer used).
    The calculation stops after the current step if eh.abort_flag is set (e.g. by the abort button).
    """
    eh = exclusion_height_analysis_object
    eh.abort_flag = False
    steps = len(eh.buffer_reflectivity_dfs)
    eh.all_exclusion_results = np.full((2, steps), np.nan)
    eh.d_n_pair_dfs = [0] * steps
    eh.height_sensitivity = np.full(steps, np.nan)
    eh.mean_exclusion_height_result = [np.nan]
    eh.mean_exclusion_RI_result = [np.nan]

    probe_curves = {}
    for step in range(steps):
        if eh.abort_flag:
            print('Exclusion height calculation aborted.')
            eh.abort_flag = False
            break
        probe_index = step // 2
        if getattr(eh, 'height_search', 'grid') == 'adaptive':
            result = adaptive_exclusion_height(eh, step, probe_index)
            store_exclusion_height_step(eh, step, result['heights'], result['buffer'], result['probe'], result['height'], result['RI'], result['sensitivity'])
        else:
            buffer_curve = calculate_exclusion_height(eh, 'buffer', step)
            if probe_index not in probe_curves:
                probe_curves[probe_index] = calculate_exclusion_height(eh, 'probe', probe_index)
            probe_curve = probe_curves[probe_index]
            height, RI = exclusion_height_crossing(eh.height_steps, buffer_curve[0], probe_curve[0])
            sensitivity = np.nan
            if np.isfinite(height):
                try:
                    sensitivity = exclusion_height_sensitivity(eh, probe_index, height, buffer_curve[0], probe_curve[0])
                except Exception as error:
                    print('Could not estimate the sensitivity of the exclusion height: {error}'.format(error=error))
            store_exclusion_height_step(eh, step, np.asarray(eh.height_steps, dtype=np.float64),
                                        np.column_stack([np.asarray(curve, dtype=np.float64) for curve in buffer_curve if len(curve)]),
                                        np.column_stack([np.asarray(curve, dtype=np.float64) for curve in probe_curve if len(curve)]),
                                        height, RI, sensitivity)

    if np.all(np.isnan(eh.all_exclusion_results[0, :])):
        eh.mean_exclusion_height_result = (np.nan, np.nan)
        eh.mean_exclusion_RI_result = (np.nan, np.nan)
    else:
        eh.mean_exclusion_height_result = (np.nanmean(eh.all_exclusion_results[0, :]), np.nanstd(eh.all_exclusion_results[0, :]))
        eh.mean_exclusion_RI_result = (np.nanmean(eh.all_exclusion_results[1, :]), np.nanstd(eh.all_exclusion_results[1, :]))
    return eh


def add_sensor_backend(session_object, data_path_, default_sensor_values, sensor_metal='Au'):

    """
    Adds sensor objects to a session object.
    :return: a sensor object
    """
    session_object.sensor_ID_count += 1
    sensor_object = Sensor(data_path_, session_object.sensor_ID_count, default_sensor_values, sensor_metal=sensor_metal)
    session_object.sensor_instances[session_object.sensor_ID_count] = sensor_object

    return sensor_object


def copy_sensor_backend(session_object, sensor_object):

    """
    Copies sensor object to a session object.
    :return: a sensor object
    """
    session_object.sensor_ID_count += 1
    copied_sensor_object = copy.deepcopy(sensor_object)
    copied_sensor_object.object_id = session_object.sensor_ID_count
    session_object.sensor_instances[session_object.sensor_ID_count] = copied_sensor_object

    return copied_sensor_object


def add_fresnel_model_object(session_object, sensor_object, data_path_, reflectivity_df_, object_name_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """
    session_object.fresnel_analysis_ID_count += 1
    analysis_object = FresnelModel(session_object, sensor_object, data_path_, reflectivity_df_, session_object.fresnel_analysis_ID_count, object_name_)
    session_object.fresnel_analysis_instances[session_object.fresnel_analysis_ID_count] = analysis_object

    return analysis_object


def add_exclusion_height_object(session_object, fresnel_object, sensorgram_df_, data_path_, object_name_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """
    session_object.exclusion_height_analysis_ID_count += 1
    analysis_object = ExclusionHeight(session_object, fresnel_object, sensorgram_df_, data_path_, session_object.exclusion_height_analysis_ID_count, object_name_)
    session_object.exclusion_height_analysis_instances[session_object.exclusion_height_analysis_ID_count] = analysis_object

    return analysis_object
