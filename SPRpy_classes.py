import datetime
import os
import scipy
import pickle
import copy
import bottleneck
import multiprocessing
from SPRpy_functions import *
from fresnel_transfer_matrix import fresnel_calculation


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name='Session', directory=os.getcwd(), current_data_path=None):
        self.name = datetime.datetime.now().__str__()[0:16].replace(':', '_') + ' ' + name
        if not os.path.exists(directory + r'\SPRpy sessions'):
            os.mkdir(directory + r'\SPRpy sessions')
        self.location = directory + r'\SPRpy sessions' + r'\{name_}'.format(name_=self.name)
        if not os.path.exists(self.location):
            os.mkdir(self.location)
        if not os.path.exists(self.location + r'\Sensors'):
            os.mkdir(self.location + r'\Sensors')
        if not os.path.exists(self.location + r'\Analysis instances'):
            os.mkdir(self.location + r'\Analysis instances')
        self.sensor_instances = {}  # NOTE: The sessions in this list are also updated when modified as current sensor object
        self.sensor_ID_count = 0
        self.fresnel_analysis_instances = {}
        self.fresnel_analysis_ID_count = 0
        self.exclusion_height_analysis_instances = {}
        self.exclusion_height_analysis_ID_count = 0
        self.current_data_path = current_data_path
        self.log = datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Welcome to SPRpy!' \
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Start your session by defining your SPR sensor layers.' \
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'You can import previous results under "File and session controls".'

    def remove_sensor(self, sensor_object_id):
        """
        Remove a sensor object from the session.
        :return:
        """
        removed = self.sensor_instances.pop(sensor_object_id)
        removed_file_path = self.location + r'\Sensors' + r'\S{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following sensor object: S{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def remove_fresnel_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.fresnel_analysis_instances.pop(analysis_object_id)
        removed_file_path = self.location + r'\Analysis instances' + r'\FM{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following analysis object: FM{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def remove_exclusion_height_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.exclusion_height_analysis_instances.pop(analysis_object_id)
        removed_file_path = self.location + r'\Analysis instances' + r'\EH{id} {name}.pickle'.format(id=removed.object_id, name=removed.name)
        os.remove(removed_file_path)
        print('Removed the following analysis object: FM{id} {name}'.format(id=removed.object_id, name=removed.name))

        return

    def save_all(self):
        """
        Saves all objects stored in the session, and the session file itself.
        :return: None
        """

        # Save session object
        with open(self.location + r'\Session file.pickle', 'wb') as save_file:
            pickle.dump(self, save_file)

        # Save sensor instances
        for sensor_id in self.sensor_instances:
            with open(self.location + r'\Sensors' + r'\S{id} {name}.pickle'.format(id=sensor_id, name=self.sensor_instances[sensor_id].name), 'wb') as save_file:
                pickle.dump(self.sensor_instances[sensor_id], save_file)

        # Save fresnel analysis instances
        for analysis_id in self.fresnel_analysis_instances:
            with open(self.location + r'\Analysis instances' + r'\FM{id} {name}.pickle'.format(id=analysis_id, name=self.fresnel_analysis_instances[analysis_id].name), 'wb') as save_file:
                pickle.dump(self.fresnel_analysis_instances[analysis_id], save_file)

        # Save exclusion height analysis instances
        for analysis_id in self.exclusion_height_analysis_instances:
            with open(self.location + r'\Analysis instances' + r'\FM{id} {name}.pickle'.format(id=analysis_id, name=self.exclusion_height_analysis_instances[analysis_id].name), 'wb') as save_file:
                pickle.dump(self.exclusion_height_analysis_instances[analysis_id], save_file)

        return

    def save_session(self):

        # Save session object
        with open(self.location + r'\Session file.pickle', 'wb') as save_file:
            pickle.dump(self, save_file)

        return

    def save_sensor(self, sensor_id):
        """
        Saves a single sensor object to the session.
        :return: None
        """

        with open(self.location + r'\Sensors' + r'\S{id} {name}.pickle'.format(id=sensor_id, name=self.sensor_instances[
            sensor_id].name), 'wb') as save_file:
            pickle.dump(self.sensor_instances[sensor_id], save_file)

        return

    def save_fresnel_analysis(self, analysis_id):
        """
        Saves a single fresnel analysis object to the session.
        :return: None
        """

        with open(self.location + r'\Analysis instances' + r'\FM{id} {name}.pickle'.format(id=analysis_id, name=self.fresnel_analysis_instances[analysis_id].name), 'wb') as save_file:
            pickle.dump(self.fresnel_analysis_instances[analysis_id], save_file)

        return

    def save_exclusion_height_analysis(self, analysis_id):
        """
        Saves a single fresnel analysis object to the session.
        :return: None
        """

        with open(self.location + r'\Analysis instances' + r'\EH{id} {name}.pickle'.format(id=analysis_id, name=self.exclusion_height_analysis_instances[analysis_id].name), 'wb') as save_file:
            pickle.dump(self.exclusion_height_analysis_instances[analysis_id], save_file)

        return

    def import_sensor(self):

        file_path_ = select_file(prompt='Select the sensor object', prompt_folder=self.location + r'\Sensors')
        self.sensor_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            sensor_object = pickle.load(import_file)

        sensor_object.object_id = self.sensor_ID_count
        self.sensor_instances[self.sensor_ID_count] = sensor_object

        return

    def import_fresnel_analysis(self):
        file_path_ = select_file(prompt='Select the analysis object', prompt_folder=self.location + r'\Analysis instances')
        self.fresnel_analysis_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            analysis_object = pickle.load(import_file)

        analysis_object.object_id = self.fresnel_analysis_ID_count
        self.fresnel_analysis_instances[analysis_object.object_id] = analysis_object

        return

    def import_exclusion_height_analysis(self):
        file_path_ = select_file(prompt='Select the analysis object', prompt_folder=self.location + r'\Analysis instances')
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

    def __init__(self, data_path_, object_id_, object_name_='Gold sensor', sensor_metal='Au', data_type='R', polarization=1):
        """
        :param data_path_: string
        :param sensor_metal: string, see options in method "set_default_optical_properties"
        :param polarization: int, default 1 or "p-polarization"
        """
        # Load sensor's default optical properties
        self.object_id = object_id_
        self.name = object_name_
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
                self.layer_thicknesses = np.array([np.NaN, 2.00, 50.00, np.NaN])
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
                self.optical_parameters = pd.DataFrame(data={'Layers': ['Prism', 'Cr', 'Au', 'Bulk'],  # NOTE: sensor_object.optical_parameters will auto-update when changing sensor_object.layer_thicknesses, sensor_object.refractive_indices or sensor_object.extinction_coefficients
                                                             'd [nm]': self.layer_thicknesses,
                                                             'n': self.refractive_indices,
                                                             'k': self.extinction_coefficients})

            case 'sio2' | 'SiO2' | 'SIO2' | 'Glass' | 'glass' | 'silica':
                # Fused silica values source: L. V. Rodríguez-de Marcos, J. I. Larruquert, J. A. Méndez, J. A. Aznárez.
                # Self-consistent optical constants of SiO2 and Ta2O5 films
                # Opt. Mater. Express 6, 3622-3637 (2016) (Numerical data kindly provided by Juan Larruquert)
                self.layer_thicknesses = np.array([np.NaN, 2.00, 50.00, 14.00, np.NaN])
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
                self.layer_thicknesses = np.array([np.NaN, 2.00, 20.00, np.NaN])
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
                self.layer_thicknesses = np.array([np.NaN, 2.00, 20.00, np.NaN])
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
        return


class FresnelModel:
    """
    This class defines how a modelled reflectivity trace behaves. Note that a different sensor object should be used for
    each layer added to the sensor!

    TODO: Each object should also have a .csv export function.

    """

    def __init__(self, sensor_object_, data_path_, reflectivity_df_, TIR_range_, scanspeed_, object_id_, object_name_):
        self.name = object_name_
        self.object_id = object_id_
        self.sensor_object = sensor_object_
        self.sensor_object_label = ''
        self.initial_data_path = data_path_
        self.measurement_data = reflectivity_df_
        self.TIR_range = TIR_range_  # Default for air. For water: (60.8, 63)
        self.scanspeed = scanspeed_  # Scan speed from .spr2 file, 1 (slow), 5 (medium) or 10 (fast)
        self.angle_range = [40, 80]
        self.ini_guess = 4
        self.bounds = [0, 50]
        self.extinction_correction = 0
        self.y_offset = 0
        self.fitted_data = None
        self.fitted_result = None

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
                                                    polarization=self.sensor_object.polarization)
        return fresnel_coefficients_

    def model_reflectivity_trace(self):
        """

        :param ini_guess:
        :param bounds:
        :param TIR_range:
        :param scanspeed:
        :return:
        """

        xdata_ = self.measurement_data['angles']
        ydata_ = self.measurement_data['ydata']

        # Calculate TIR angle and bulk refractive index
        TIR_angle, TIR_fitted_angles, TIR_fitted_ydata = TIR_determination(xdata_, ydata_, self.TIR_range, self.scanspeed)
        self.sensor_object.refractive_indices[-1] = self.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)

        # Add extinction correction to fitted surface layer extinction value
        extinction_corrected = self.sensor_object.extinction_coefficients
        extinction_corrected[0] += self.extinction_correction  # Correct Prism layer extinction

        # Selecting a range of measurement data to use for fitting, and including an offset in reflectivity (iterated 3 times)
        selection_xdata_ = xdata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])]

        for offset_ind in range(3):
            selection_ydata_ = ydata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])] - self.y_offset

            # Perform the fitting
            result = scipy.optimize.least_squares(fresnel_calculation,
                                                  self.ini_guess,
                                                  bounds=self.bounds,
                                                  kwargs={'fitted_layer_index': self.sensor_object.fitted_layer_index,
                                                          'wavelength': self.sensor_object.wavelength,
                                                          'layer_thicknesses': self.sensor_object.layer_thicknesses,
                                                          'n_re': self.sensor_object.refractive_indices,
                                                          'n_im': extinction_corrected,
                                                          'angles': selection_xdata_,
                                                          'ydata': selection_ydata_,
                                                          'ydata_type': self.sensor_object.data_type,
                                                          'polarization': self.sensor_object.polarization}
                                                  )
            # Collect the results from least_squares object and calculate corresponding fresnel coefficients
            self.fitted_result = result['x'][0]
            fresnel_coefficients = fresnel_calculation(self.fitted_result,
                                                       fitted_layer_index=self.sensor_object.fitted_layer_index,
                                                       angles=selection_xdata_,
                                                       wavelength=self.sensor_object.wavelength,
                                                       layer_thicknesses=self.sensor_object.layer_thicknesses,
                                                       n_re=self.sensor_object.refractive_indices,
                                                       n_im=extinction_corrected,
                                                       ydata=None,
                                                       ydata_type='R',
                                                       polarization=1
                                                       )
            if offset_ind < 2:
                # Calculate new y_offset
                self.y_offset = self.y_offset + np.min(selection_ydata_) - np.min(fresnel_coefficients)

        # Shift fresnel coefficients back to measurement data  level
        fresnel_ydata = fresnel_coefficients + self.y_offset

        # Compile into fresnel_coefficients data frame
        self.fitted_data = pd.DataFrame(data={'angles': selection_xdata_, 'ydata': fresnel_ydata})

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

    def __init__(self, fresnel_object_, sensorgram_df_, data_path_,  object_id_, object_name_):
        self.name = object_name_
        self.object_id = object_id_
        self.fresnel_object = fresnel_object_
        self.fresnel_object_label = 'Fresnel background: FM{analysis_number} {analysis_name}'.format(
                                                                    analysis_number=fresnel_object_.object_id,
                                                                    analysis_name=fresnel_object_.name)
        self.sensor_object = fresnel_object_.sensor_object
        self.initial_data_path = data_path_
        self.sensorgram_data = sensorgram_df_
        self.sensorgram_offset_ind = 0
        self.d_n_pair_resolution = 200
        self.height_steps = np.linspace(0, 200, self.d_n_pair_resolution)
        self.points_below_SPR_min_ind = None
        self.points_above_SPR_min_ind = None
        self.injection_points = []
        self.buffer_points = []
        self.probe_points = []
        self.SPR_vs_TIR_dfs = []  # List of dataframes with labels 'SPR angles' and 'TIR angles' for indexing each step result
        self.buffer_reflectivity_dfs = []  # Use labels 'buffer reflectivity' and 'buffer angles' (and likewise probe) for indexing
        self.buffer_bulk_RIs = []  # Calculated from TIR angle of each reflectivity DF
        self.probe_reflectivity_dfs = []  # Use labels 'buffer reflectivity' and 'buffer angles' (and likewise probe) for indexing
        self.probe_bulk_RIs = []  # Calculated from TIR angle of each reflectivity DF
        self.buffer_d_n_pair_dfs = []  # Use labels 'buffer thickness' and 'buffer refractive index' (and likewise probe) for indexing
        self.probe_d_n_pair_dfs = []  # Use labels 'buffer thickness' and 'buffer refractive index' (and likewise probe) for indexing
        self.all_exclusion_RI_steps = []
        self.all_exclusion_height_results = []
        self.mean_exclusion_height_result = None

    # TODO: The methods running calculations here need to use mutliprocessing and whould be run inside background callbacks in the dash app to prevent timeout after 30s of calculations.
    # TODO: Make sure the quality of fit for each d,n pair can be viewed with a pagination passing over each injection.
    #  Should probably do a list or dictionary backend-wise containing each injection.
    # TODO: Backend-wise, do stepping in height and perform fitting of the refractive index (opposite to matlab script).
    #  It means user enters a plausible range for the heights instead, which is easier to relate to and weird
    #  bulk effects are automatically detected (like PEG+MCH demonstrated).
    # TODO: Pull most data from current sensor object and current fresnel analysis object (which should have
    #  been performed on an angular trace containing buffer+layer, the resulting height corresponds to a 0 %
    #  swollen version of the layer), but needs input for a range of plausible heights. This range can be calculated as
    #  a default based on the dry height of the layer in air, but maybe easier to just add a tooltip stating that fact than implementing it.
    # TODO: Make it so the progress bar updates with each finished injection, and that the results page is
    #  updated with each injection so it can be aborted if necessary
    # TODO: There is no need to include different offsets for buffer and probe, they can simply use the offset from fresnel background object (if it is modeled from the liquid!)

    def initialize_model(self, ydata_df):

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

        # Calculate average reflectivity traces based on selected points
        bufferpoint_index = 0
        for reflectivity_index in range(int(len(self.buffer_points) / 2)):

            sliced_buffer_reflectivity_spectras = ydata_df.iloc[self.buffer_points[bufferpoint_index][0]:self.buffer_points[bufferpoint_index + 1][0], :]  # Selecting all spectras between the pairwise selected buffer points
            mean_buffer_reflectivity = sliced_buffer_reflectivity_spectras.mean(axis=0).squeeze()

            # Calculate TIR and bulk RI for each mean spectra
            buffer_TIR_angle, _, _ = TIR_determination(background_angles.to_numpy(), mean_buffer_reflectivity.to_numpy(), self.fresnel_object.TIR_range, self.fresnel_object.scanspeed)
            self.buffer_bulk_RIs.append(self.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * buffer_TIR_angle))

            # Calculate appropriate range selection
            buffer_reflectivity_minimum_ind = pd.Series(bottleneck.move_mean(mean_buffer_reflectivity.to_numpy(), window=4, min_count=1)).idxmin()
            self.buffer_reflectivity_dfs.append(pd.DataFrame(data={'reflectivity': mean_buffer_reflectivity.iloc[buffer_reflectivity_minimum_ind - self.points_below_SPR_min_ind:buffer_reflectivity_minimum_ind + self.points_above_SPR_min_ind],
                                                                   'angles': background_angles.iloc[buffer_reflectivity_minimum_ind - self.points_below_SPR_min_ind:buffer_reflectivity_minimum_ind + self.points_above_SPR_min_ind]
                                                                   })
                                                )

            # Next pair of buffer point indices
            bufferpoint_index += 2

        probepoint_index = 0
        for reflectivity_index in range(int(len(self.probe_points) / 2)):
            sliced_probe_reflectivity_spectras = ydata_df.iloc[self.probe_points[probepoint_index][0]:self.probe_points[probepoint_index + 1][0], :]  # Selecting all spectras between the pairwise selected probe point
            mean_probe_reflectivity = sliced_probe_reflectivity_spectras.mean(axis=0).squeeze()

            # Calculate TIR and bulk RI for each mean spectra
            probe_TIR_angle, _, _ = TIR_determination(background_angles.to_numpy(), mean_probe_reflectivity.to_numpy(), self.fresnel_object.TIR_range, self.fresnel_object.scanspeed)
            self.probe_bulk_RIs.append(self.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * probe_TIR_angle))

            # Calculate appropriate range selection
            probe_reflectivity_minimum_ind = pd.Series(bottleneck.move_mean(mean_probe_reflectivity.to_numpy(), window=4, min_count=1)).idxmin()
            self.probe_reflectivity_dfs.append(pd.DataFrame(data={'reflectivity': mean_probe_reflectivity.iloc[probe_reflectivity_minimum_ind - self.points_below_SPR_min_ind:probe_reflectivity_minimum_ind + self.points_above_SPR_min_ind],
                                                                  'angles': background_angles.iloc[probe_reflectivity_minimum_ind - self.points_below_SPR_min_ind:probe_reflectivity_minimum_ind + self.points_above_SPR_min_ind]
                                                                  })
                                               )

            # Next pair of probe point indices
            probepoint_index += 2


        # Create SPR vs TIR data frames
        injectionpoint_index = 0
        offset_SPR_sensorgram = self.sensorgram_data['SPR angle'] - self.sensorgram_data.loc[self.sensorgram_offset_ind, 'SPR angle']
        offset_TIR_sensorgram = self.sensorgram_data['TIR angle'] - self.sensorgram_data.loc[self.sensorgram_offset_ind, 'TIR angle']

        for reflectivity_index in range(int(len(self.injection_points) / 2)):

            self.SPR_vs_TIR_dfs.append(pd.DataFrame(data={'SPR angles': offset_SPR_sensorgram.iloc[self.injection_points[injectionpoint_index][0]:self.injection_points[injectionpoint_index + 1][0]],
                                                          'TIR angles': offset_TIR_sensorgram.iloc[self.injection_points[injectionpoint_index][0]:self.injection_points[injectionpoint_index + 1][0]],
                                                          })
                                       )

            # Next pair of injection point indices
            injectionpoint_index += 2

        return


def model_buffer_reflectivity_trace(exclusion_height_analysis_object, step_index_, height):
    """

    :param ini_guess:
    :param bounds:
    :param TIR_range:
    :param scanspeed:
    :return:
    """

    # TODO: Adapt this code to model layer RI instead of height

    selection_xdata_ = exclusion_height_analysis_object.buffer_reflectivity_dfs[step_index_]['angles']  # This should already be the selected range
    selection_ydata_ = exclusion_height_analysis_object.buffer_reflectivity_dfs[step_index_]['reflectivity']

    # Calculate TIR angle and bulk refractive index
    TIR_angle, TIR_fitted_angles, TIR_fitted_ydata = TIR_determination(selection_xdata_, selection_ydata_, exclusion_height_analysis_object.fresnel_object.TIR_range, exclusion_height_analysis_object.fresnel_object.scanspeed)
    exclusion_height_analysis_object.sensor_object.refractive_indices[-1] = exclusion_height_analysis_object.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)

    # Selecting a range of measurement data to use for fitting, and including an offset in reflectivity (iterated 3 times)
    # for offset_ind in range(3): # TODO: There is no need to include different offsets for buffer and probe
    offset_ydata_ = selection_ydata_ - exclusion_height_analysis_object.fresnel_object.y_offset

    # Perform the fitting
    result = scipy.optimize.least_squares(fresnel_calculation,
                                          exclusion_height_analysis_object.fresnel_object.ini_guess,
                                          bounds=exclusion_height_analysis_object.fresnel_object.bounds,
                                          kwargs={'fitted_layer_index': exclusion_height_analysis_object.sensor_object.fitted_layer_index,
                                                  'wavelength': exclusion_height_analysis_object.sensor_object.wavelength,
                                                  'layer_thicknesses': exclusion_height_analysis_object.sensor_object.layer_thicknesses,
                                                  'n_re': exclusion_height_analysis_object.sensor_object.refractive_indices,
                                                  'n_im': exclusion_height_analysis_object.sensor_object.extinction_coefficients,
                                                  'angles': selection_xdata_,
                                                  'ydata': offset_ydata_,
                                                  'ydata_type': exclusion_height_analysis_object.sensor_object.data_type,
                                                  'polarization': exclusion_height_analysis_object.sensor_object.polarization}
                                          )
    # Collect the results from least_squares object and calculate corresponding fresnel coefficients
    fitted_result = result['x'][0]
    fresnel_coefficients = fresnel_calculation(fitted_result,
                                               fitted_layer_index=exclusion_height_analysis_object.sensor_object.fitted_layer_index,
                                               angles=selection_xdata_,
                                               wavelength=exclusion_height_analysis_object.sensor_object.wavelength,
                                               layer_thicknesses=exclusion_height_analysis_object.sensor_object.layer_thicknesses,
                                               n_re=exclusion_height_analysis_object.sensor_object.refractive_indices,
                                               n_im=exclusion_height_analysis_object.sensor_object.extinction_coefficients,
                                               ydata=None,
                                               ydata_type='R',
                                               polarization=1
                                               )
    # if offset_ind < 2: # TODO: There is no need to include different offsets for buffer and probe
    #     # Calculate new y_offset
    #     exclusion_height_analysis_object.fresnel_object.y_offset = exclusion_height_analysis_object.fresnel_object.y_offset + np.min(offset_ydata_) - np.min(fresnel_coefficients)

    # Shift fresnel coefficients back to measurement data  level
    fresnel_ydata = fresnel_coefficients + exclusion_height_analysis_object.fresnel_object.y_offset

    # Compile into fresnel_coefficients data frame
    exclusion_height_analysis_object.fitted_data = pd.DataFrame(data={'angles': selection_xdata_, 'ydata': fresnel_ydata})

    return exclusion_height_analysis_object.fitted_data


def calculate_exclusion_height(exclusion_height_analysis_object, buffer_or_probe_flag, data_frame_index):
    """
    TODO: this should be the main function that is called to calculate exclusion heights. It should output a list with the RI results for each step.

    :param exclusion_height_analysis_object: object containing all parameters and data
    :param data_frame_index: index of dataframe and buffer or probe RI bulk value
    :return RI_results: list of calculated RI results for each height step
    """
    # Check if calculations should be performed on buffer or probe data
    if buffer_or_probe_flag == 'buffer':
        # TODO: Don't forget to swap the buffer RI bulk value of exclusion_height_analysis_object.fresnel_object.refractive_indeces to the exclusion_height_analysis_object.buffer_bulk_RIs[data_frame_index]
        RI_results = []

        return RI_results
    elif buffer_or_probe_flag == 'probe':
        # TODO: Don't forget to swap the probe RI bulk value of exclusion_height_analysis_object.fresnel_object.refractive_indeces to the exclusion_height_analysis_object.probe_bulk_RIs[data_frame_index]
        RI_results = []

        return RI_results
    else:
        raise ValueError('only buffer or probe allowed')


def exclusion_height_process(exclusion_height_analysis_object, buffer_or_probe_flag, data_frame_index, connection):
    '''
    This function initiates the calculations and send back the result to the spawning process

    :param exclusion_height_analysis_object: object containing all parameters and data
    :param buffer_or_probe_flag: either 'buffer' or 'probe'
    :param data_frame_index: index of dataframe
    :param connection: child pipe connection object from multiprocessing.Pipe()
    :return: None
    '''

    result = calculate_exclusion_height(exclusion_height_analysis_object, buffer_or_probe_flag, data_frame_index)
    connection.send(result)
    connection.close()


def process_all_exclusion_heights(exclusion_height_analysis_object):
    """
    TODO: This function should start multiprocessing of several calls to calculate_exclusion_height() for different injection steps
     I think the best way to get the results are using the pipe method, where each process writes to a pipe and the main process reads from the pipe.
     The main process should then update the progress bar and the results page with each finished injection.
    TODO: See chatGPT history for multiprocessing example using Pipe connection objects to recieve the result from each process.

    """

    # TODO: Need to check for the number of cores available and limit the number of processes to that number.
    #  Also need to check if the number of processes is larger than the number of steps and limit it to the number of steps.
    #  The number of steps is the number of injections, which is len(exclusion_height_analysis_object.injection_points) / 2
    #  The number of cores can be found using multiprocessing.cpu_count()-2 (to save cores for other applications)
    #  The number of processes should be limited to the number of cores, or the number of steps, whichever is smallest.

    buffer_connections = []
    probe_connections = []
    buffer_processes = []
    probe_processes = []

    buffer_index = 0
    probe_index = 0
    process_step = 0
    while process_step < min(len(exclusion_height_analysis_object.injection_points) / 2, multiprocessing.cpu_count()-2):
        # Setup buffer process
        buffer_parent_conn, buffer_child_conn = multiprocessing.Pipe()
        buffer_connections.append(buffer_parent_conn)
        buffer_process = multiprocessing.Process(target=exclusion_height_process, args=(exclusion_height_analysis_object, 'buffer', buffer_index, buffer_child_conn))
        buffer_processes.append(buffer_process)
        buffer_index += 1

        # Setup probe process
        probe_parent_conn, probe_child_conn = multiprocessing.Pipe()
        probe_connections.append(probe_parent_conn)
        probe_process = multiprocessing.Process(target=exclusion_height_process, args=(exclusion_height_analysis_object, 'probe', probe_index, probe_child_conn))
        probe_processes.append(probe_process)
        probe_index += 1

        # Start processes
        buffer_process.start()
        probe_process.start()

        # Update step
        process_step += 2

    # If there are more steps than cpu cores, wait for the next two processes to finish and start new processes for the remaining steps
    if process_step < len(exclusion_height_analysis_object.injection_points) / 2:
        for process_index in range(len(exclusion_height_analysis_object.injection_points) / 2 - (multiprocessing.cpu_count()-2)):
            buffer_processes[process_index].join()
            probe_processes[process_index].join()

            # Setup buffer process
            buffer_parent_conn, buffer_child_conn = multiprocessing.Pipe()
            buffer_connections.append(buffer_parent_conn)
            buffer_process = multiprocessing.Process(target=exclusion_height_process, args=(exclusion_height_analysis_object, 'buffer', buffer_index, buffer_child_conn))
            buffer_processes.append(buffer_process)
            buffer_index += 1

            # Setup probe process
            probe_parent_conn, probe_child_conn = multiprocessing.Pipe()
            probe_connections.append(probe_parent_conn)
            probe_process = multiprocessing.Process(target=exclusion_height_process, args=(exclusion_height_analysis_object, 'probe', probe_index, probe_child_conn))
            probe_processes.append(probe_process)
            probe_index += 1

            # Start processes
            buffer_process.start()
            probe_process.start()


    # Wait for each process to finish and collect and record results
    for buffer_process, probe_process, buffer_conn, probe_conn in zip(buffer_processes, probe_processes, buffer_connections, probe_connections):
        buffer_process.join()
        probe_process.join()

        buffer_RI_result = buffer_conn.recv()
        probe_RI_result = probe_conn.recv()

        exclusion_height_analysis_object.all_exclusion_RI_steps.append(pd.DataFrame(data={'buffer RI': buffer_RI_result, 'probe RI': probe_RI_result}))

        # Calculate exclusion height from buffer and probe height steps and RI result intersection and add to exclusion_height_analysis_object.all_exclusion_height_results
        # TODO: Base this on matlab intersection code for finding the exclusion height



def add_sensor_backend(session_object, data_path_, sensor_metal='Au', polarization=1):

    """
    Adds sensor objects to a session object.
    :return: a sensor object
    """
    session_object.sensor_ID_count += 1
    sensor_object = Sensor(data_path_, session_object.sensor_ID_count, sensor_metal=sensor_metal, polarization=polarization)
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


def add_fresnel_model_object(session_object, sensor_object, data_path_, reflectivity_df_, TIR_range_, scanspeed_, object_name_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """
    session_object.fresnel_analysis_ID_count += 1
    analysis_object = FresnelModel(sensor_object, data_path_, reflectivity_df_, TIR_range_, scanspeed_, session_object.fresnel_analysis_ID_count, object_name_)
    session_object.fresnel_analysis_instances[session_object.fresnel_analysis_ID_count] = analysis_object

    return analysis_object


def add_exclusion_height_object(session_object, fresnel_object, sensorgram_df_, data_path_, object_name_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """
    session_object.exclusion_height_analysis_ID_count += 1
    analysis_object = ExclusionHeight(fresnel_object, sensorgram_df_, data_path_, session_object.exclusion_height_analysis_ID_count, object_name_)
    session_object.exclusion_height_analysis_instances[session_object.exclusion_height_analysis_ID_count] = analysis_object

    return analysis_object
