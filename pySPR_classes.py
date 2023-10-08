import datetime
import os
import scipy
import pickle
import copy
from pySPR_functions import *
from fresnel_transfer_matrix import fresnel_calculation


class Session:

    """
    A base class for storing and loading a group of measurements for a particular session in one file. This should be
    the first thing that a user is prompted for before they start their analysis.
    """

    def __init__(self, name='Experiments', directory=os.getcwd(), current_data_path=None):
        self.name = datetime.datetime.now().__str__()[0:16].replace(':', '_') + ' ' + name
        if not os.path.exists(directory + r'\pySPR sessions'):
            os.mkdir(directory + r'\pySPR sessions')
        self.location = directory + r'\pySPR sessions' + r'\{name_}'.format(name_=self.name)
        if not os.path.exists(self.location):
            os.mkdir(self.location)
        if not os.path.exists(self.location + r'\Sensors'):
            os.mkdir(self.location + r'\Sensors')
        if not os.path.exists(self.location + r'\Analysis instances'):
            os.mkdir(self.location + r'\Analysis instances')
        self.sensor_instances = {}  # NOTE: The sessions in this list are also updated when modified as current sensor object
        self.sensor_ID_count = 0
        self.analysis_instances = {}
        self.analysis_ID_count = 0
        self.current_data_path = current_data_path
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

        return

    def remove_analysis(self, analysis_object_id):
        """
        Remove an analysis object from the session.
        :return:
        """
        removed = self.analysis_instances.pop(analysis_object_id)
        print('Removed the following analysis object: ' + str(removed))

        return

    def save_all(self):
        """
        Saves all objects stored in the session, and the session file itself.
        :return: None
        """

        # Save session object
        with open(self.location + r'\Session_file.pickle', 'wb') as save_file:
            pickle.dump(self, save_file)

        # Save sensor instances
        for sensor_id in self.sensor_instances:
            with open(self.location + r'\Sensors' + r'\Sensor_{id}.pickle'.format(id=sensor_id), 'wb') as save_file:
                pickle.dump(self.sensor_instances[sensor_id], save_file)

        # Save analysis instances
        for analysis_name in self.analysis_instances:
            with open(self.location + r'\Analysis instances' + r'\{name}.pickle'.format(name=analysis_name), 'wb') as save_file:
                pickle.dump(self.analysis_instances[analysis_name], save_file)

        return

    def save_session(self):

        # Save session object
        with open(self.location + r'\Session_file.pickle', 'wb') as save_file:
            pickle.dump(self, save_file)

        return

    def save_sensor(self, sensor_id):
        """
        Saves a single sensor object to the session.
        :return: None
        """

        with open(self.location + r'\Sensors' + r'\Sensor_{id}.pickle'.format(id=sensor_id), 'wb') as save_file:
            pickle.dump(self.sensor_instances[sensor_id], save_file)

        return

    def save_analysis(self, analysis_name):
        """
        Saves a single analysis object to the session.
        :return: None
        """

        with open(self.location + r'\Analysis instances' + r'\{name}.pickle'.format(name=analysis_name), 'wb') as save_file:
            pickle.dump(self.analysis_instances[analysis_name], save_file)

        return

    def import_sensor(self):

        file_path_ = select_file(prompt='Select the sensor object', prompt_folder=self.location + r'\Sensors')
        self.sensor_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            sensor_object = pickle.load(import_file)

        sensor_object.object_id = self.sensor_ID_count
        self.sensor_instances[self.sensor_ID_count] = sensor_object

        return

    def import_analysis(self):
        file_path_ = select_file(prompt='Select the analysis object', prompt_folder=self.location + r'\Analysis instances')
        self.analysis_ID_count += 1

        with open(file_path_, 'rb') as import_file:
            analysis_object = pickle.load(import_file)

        analysis_object.object_id = self.analysis_ID_count
        self.analysis_instances[analysis_object.object_id] = analysis_object

        return


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
        return

    # TODO: These functions are not necessary
    # def add_material_layer(self, thickness, n_re, n_im, layer_index_=-1):
    #     """
    #     Add additional layers on top of the sensor (before bulk medium). not used by dash app
    #     :return:
    #     """
    #     # Use negative indexing for this, so it always add the layer on the top no matter what was there previously
    #     self.layer_thicknesses = np.insert(self.layer_thicknesses, layer_index_, thickness)
    #     self.refractive_indices = np.insert(self.refractive_indices, layer_index_, n_re)
    #     self.extinction_coefficients = np.insert(self.extinction_coefficients, layer_index_, n_im)
    #     self.fitted_layer = 'h_surf'
    #     print('Sensor thicknesses: ', self.layer_thicknesses)
    #     print('Sensor refractive indices: ', self.refractive_indices)
    #     print('Sensor extinction coefficients: ', self.extinction_coefficients)
    #
    # def remove_material_layer(self, layer_index_):
    #
    #     """
    #     Removes a layer from a sensor. (Not used by dash app UI.)
    #
    #     :param layer_index_: int, which layer to remove (starting from 1)
    #     :return:
    #     """
    #
    #     self.layer_thicknesses = np.delete(self.layer_thicknesses, layer_index_-1, axis=0)
    #     if len(self.layer_thicknesses) == 4:
    #         self.fitted_layer = 'n_im_metal'
    #     self.refractive_indices = np.delete(self.refractive_indices, layer_index_-1, axis=0)
    #     self.extinction_coefficients = np.delete(self.extinction_coefficients, layer_index_-1, axis=0)
    #
    #     print('Sensor thicknesses: ', self.layer_thicknesses)
    #     print('Sensor refractive indices: ', self.refractive_indices)
    #     print('Sensor extinction coefficients: ', self.extinction_coefficients)


class ModelledReflectivityTrace:
    """
    This class defines how a modelled reflectivity trace behaves. Note that a different sensor object should be used for
    each layer added to the sensor!

    TODO: Each object should also have a .csv export function.
    TODO: IN dash app, add functionality to update current sensor object with the model object optical parameters
    TODO: Add button and graph object to dashapp for the calculate_fresnel_trace method
    TODO: Add attributes for storing calculated fresnel traces as part of results
    """

    def __init__(self, sensor_object_, data_path_, TIR_range_, angle_range_, scanspeed_, name_):
        self.name = name_
        self.sensor_object = sensor_object_  # TODO: Check if this updates the sensor object and if an output callback can update the sensor data table with the result
        self.data_path = data_path_
        self.fit_result = None
        self.y_offset = 0
        self.TIR_range = TIR_range_  # Default for air. For water: (60.8, 63)
        self.angle_range = angle_range_
        self.scanspeed = scanspeed_  # Scan speed from .spr2 file, 1 (slow), 5 (medium) or 10 (fast)

    def calculate_fresnel_trace(self, angles_=np.linspace(39, 50, 1567)):

        fresnel_coefficients_ = fresnel_calculation(None,
                                                    angles=angles_,
                                                    fitted_layer_index=self.sensor_object.fitted_layer_index,
                                                    wavelength=self.sensor_object.wavelength,
                                                    layer_thicknesses=self.sensor_object.layer_thicknesses,
                                                    n_re=self.sensor_object.refractive_indices,
                                                    n_im=self.sensor_object.extinction_coefficients,
                                                    ydata=None,
                                                    ydata_type=self.sensor_object.data_type,
                                                    polarization=self.sensor_object.polarization)
        return fresnel_coefficients_

    def model_reflectivity_trace(self, current_data_path, reflectivity_df, ini_guess, bounds):
        """

        :param ini_guess:
        :param bounds:
        :param TIR_range:
        :param scanspeed:
        :return:
        """

        # Check if current data path matches data_path when object was first initialized, otherwise load previous data
        if current_data_path == self.data_path:
            xdata_ = reflectivity_df['angles']
            ydata_ = reflectivity_df['ydata']

        else:
            _, _, _, _, _, reflectivity_df_ = load_csv_data(path=self.data_path)
            xdata_ = reflectivity_df_['angles']
            ydata_ = reflectivity_df_['ydata']

        # Calculate TIR angle and bulk refractive index
        TIR_angle, TIR_fitted_angles, TIR_fitted_ydata = TIR_determination(xdata_, ydata_, self.TIR_range, self.scanspeed)
        self.sensor_object.refractive_indices[-1] = self.sensor_object.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)  # TODO: Currently, the sensor bulk RI in optical parameters do not update according to the TIR angle, or?

        # Selecting a range of measurement data to use for fitting, and including an offset in reflectivity
        selection_xdata_ = xdata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])]
        selection_ydata_ = ydata_[(xdata_ >= self.angle_range[0]) & (xdata_ <= self.angle_range[1])] - self.y_offset

        # Perform the fitting
        result = scipy.optimize.least_squares(fresnel_calculation,
                                              ini_guess,
                                              bounds=bounds,
                                              kwargs={'fitted_layer_index': self.sensor_object.fitted_layer_index,
                                                      'wavelength': self.sensor_object.wavelength,
                                                      'layer_thicknesses': self.sensor_object.layer_thicknesses,
                                                      'n_re': self.sensor_object.refractive_indices,
                                                      'n_im': self.sensor_object.extinction_coefficients,
                                                      'angles': selection_xdata_,
                                                      'ydata': selection_ydata_,
                                                      'ydata_type': self.sensor_object.data_type,
                                                      'polarization': self.sensor_object.polarization}
                                              )
        # Collect the results from least_squares object and calculate corresponding fresnel coefficients
        self.fit_result = result['x']
        fresnel_coefficients = fresnel_calculation(self.fit_result,
                                                   fitted_layer_index=self.sensor_object.fitted_layer_index,
                                                   angles=selection_xdata_,
                                                   wavelength=self.sensor_object.wavelength,
                                                   layer_thicknesses=self.sensor_object.layer_thicknesses,
                                                   n_re=self.sensor_object.refractive_indices,
                                                   n_im=self.sensor_object.extinction_coefficients,
                                                   ydata=None,
                                                   ydata_type='R',
                                                   polarization=1
                                                   )

        return self.fit_result, selection_xdata_, fresnel_coefficients

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


class NonInteractingProbe(ModelledReflectivityTrace):

    """

    """

    def __init__(self, sensor_object, data_path_, TIR_range_, angle_range_, scanspeed_, object_id_, ydata_type='R'):
        super().__init__(sensor_object, data_path_, TIR_range_, angle_range_, scanspeed_, object_id_)  # Initializes the same way as parent objects, to shorten code
        self.object_id = object_id_
        self.ydata_type = ydata_type

    # TODO: The methods running calculations here need to use background callbacks (https://dash.plotly.com/background-callbacks)
    def calculate_fit(self):
        pass


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


def add_modelled_reflectivity_trace(session_object, sensor_object, data_path_, TIR_range_, angle_range_, scanspeed_, object_name_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    analysis_object = ModelledReflectivityTrace(sensor_object, data_path_, TIR_range_, angle_range_, scanspeed_, object_name_)
    session_object.analysis_instances[object_name_] = analysis_object

    return analysis_object
