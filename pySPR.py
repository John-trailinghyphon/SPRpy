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

# TODO: It is probably smartest to build the dash webapp last, since all functionality will be tied to it
#  and piped through it. The callbacks need to be made through function calls, targeting the active analysis objects.
#  Best to flush out the underlying logic first, but also think a bit about how the dash app layout will be?

# TODO: Next, it is time to populate the Modelled/FittedReflectivityTrace classes with methods and attributes for
#  performing the fresnel_calculation() function. Note that the calculation results should be saved in the object, but
#  it is better to avoid saving the raw data, unless that is all that was used (simply plotting reflectivity trace for
#  instance).

# TODO: For the non-interacting probe method where a previously calculated background is used in the calculations, there
#  is a need of some way to select a previous analysis. The easiest way is probably to make sure that each analysis is
#  named something? I imagine that the user selects from a list of names of each analysis (and there should be a default
#  name that auto-increments.

# TODO: There should be a way to load objects from one previous session into the active session. So that the background
#  does not have to be remade every time for instance.

# Regarding the dash app below

# TODO: My overall design vision is to have DIV elements with buttons, dropdown menus and graphs to control and visualize
#  each different type of analysis (like fresnel modelling of reflectivity traces, Non-interacting probe method etc.).
#  Maybe toggle buttons could be used to control which analysis method is used? This would also change the layout of the
#  displayed divs under the "main control DIV".


# TODO: "Main control DIV". Need buttons for controlling a lot of things:
#  - Which measurement file to load data from (.csv)
#  - Choosing analysis method (default should be fitting reflectivity traces, like background)
#  - Session control
#     * loading previous session
#     * importing previous session into current one
#     * importing sensor and analysis objects from previous session (like background)
#     * saving session
#     * removing sensor objects and analysis objects from the active session
#     *
#  - Adding sensor object
#     * adding material layers to sensor (maybe using a table interface)
#     * modifying existing layers
#     * removing layers
#  - Adding new analysis object (starting an analysis DIV)
#     * Selecting between different available types, which will change the analysis interface DIV and its options
#     * run calculations, fitting, plotting, selecting values, etc.
#  - Exporting the finished analysis as a HTML file with retained interactivity (omitting control DIV elements). This
#  can then be added to obsidian notes for instance, allowing straight forward result documentation.

#  TODO: Dash concepts to implement:
#   * DataTable to display current sensor parameters
#   * DataTable to display fitting parameters

import fresnel_transfer_matrix as ftm
import numpy as np
import datetime
import os
import sys
from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
import pandas as pd
import dash
import dash_bootstrap_components as dbc


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
            + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'You can load previous sessions using the "File" tab.'

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
    # TODO: Is this class actually necessary? If the dash_table.DataTable() class is used, then updating the layer
    #  properties from there is actually easier? Maybe default values should be determined from that instead?
    """
      An SPR measurement typically have some things in common, such as the sensor layers, measured angles,
    measured reflectivity, measurement time, etc. This information can be shared between different analysis methods for
    one measurement. This class serves as a basis for describing the current sensor, containing information about its
    layers and their optical properties.
    """

    def __init__(self, data_path_, object_id_, sensor_metal='Au', polarization=1):
        """
        :param data_path_: string
        :param sensor_metal: string, see options in method "set_default_optical_properties"
        :param polarization: int, default 1 or "p-polarization"
        """
        # Load sensor's default optical properties
        self.object_id = object_id_
        self.polarization = polarization
        self.wavelength = int(data_path_[-9:-6])
        self.sensor_metal = sensor_metal
        self.set_default_optical_properties(self.sensor_metal)

    def set_default_optical_properties(self, sensor_metal):

        # These default parameters should be set based on material layer and wavelength from loaded .csv file
        match sensor_metal:
            case 'Au' | 'gold' | 'Gold' | 'GOLD':
                self.layer_thicknesses = np.array([np.NaN, 2, 50, np.NaN])
                self.fitted_layer = 'n_im_metal'
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
                self.fitted_layer = 'h_SiO2'
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
                self.fitted_layer = 'n_im_metal'
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
                self.fitted_layer = 'n_im_metal'
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
        self.fitted_layer = 'h_surf'
        print('Sensor thicknesses: ', self.layer_thicknesses)
        print('Sensor refractive indices: ', self.refractive_indices)
        print('Sensor extinction coefficients: ', self.extinction_coefficients)

    def remove_material_layer(self, layer_index_):

        """
        Removes a layer from a sensor.
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

    Each object should also have a .csv export function.
    """

    def __init__(self, sensor_object, data_path_, object_id_):
        self.object_id = object_id_
        self.sensor_id = sensor_object.object_id
        self.polarization = sensor_object.polarization
        self.wavelength = sensor_object.wavelength
        self.layer_thicknesses = sensor_object.layer_thicknesses
        self.fitted_layer = sensor_object.fitted_layer
        self.refractive_indices = sensor_object.refractive_indices
        self.extinction_coefficients = sensor_object.extinction_coefficients
        self.data_path = data_path_

    def calculate_trace(self,  xdata_, ydata_, ini_guess, lower_bond, upper_bound):
        # TODO: This must handle data correctly, in that the data path should be checked to match the current global
        #  data path so that the calculation is performed with the right data. If the self.data_path_ attribute doesn't
        #  match the global current data path then load in the correct one.
        ftm.fresnel_calculation(ini_guess,
                                angles=xdata_,
                                layer=self.fitted_layer,
                                wavelength=self.wavelength,
                                layer_thicknesses=self.layer_thicknesses,
                                n_re=self.refractive_indices,
                                n_im=self.extinction_coefficients,
                                ydata=None,
                                ydata_type='R',
                                polarization=1
                                )

    def plot_reflectivity_trace(self, xdata_, ydata_, time_index=0, rlines_=None):
        """

        :param xdata_:
        :param ydata_:
        :param time_index:
        :param rlines_:
        :return:
        """
        pass

    def plot_sensorgram(self, time_=None, ydata_=None, rlines_=None):
        """

        """
        pass

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


def add_sensor(session_handle, data_path_, sensor_metal='Au', polarization=1):
    """
    Adds sensor objects to a session object.
    :return: a sensor object
    """
    id_ = next(session_handle.sensor_ID)
    sensor_object = Sensor(data_path_, id_, sensor_metal=sensor_metal, polarization=polarization)
    session_handle.sensor_instances[id_] = sensor_object

    return sensor_object


def add_modelled_reflectivity_trace(session_handle, sensor_object, data_path_):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_handle.analysis_ID)
    analysis_object = ModelledReflectivityTrace(sensor_object, data_path_, id_)
    session_handle.analysis_instances[id_] = analysis_object

    return analysis_object


def add_fitted_reflectivity_trace(session_handle, sensor_object, data_path_, ydata_type='R'):
    """
    Adds analysis objects to a session object.
    :return: an analysis object
    """

    id_ = next(session_handle.analysis_ID)
    analysis_object = FittedReflectivityTrace(sensor_object, data_path_, id_, ydata_type=ydata_type)
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


def load_csv_data():
    print('Select the measurement data file (.csv)')
    data_path_ = askopenfilename(title='Select the measurement data file', filetypes=[('CSV files', '*.csv')])

    # Load in the measurement data from a .csv file
    data_table = pd.read_csv(data_path_, delimiter=';', skiprows=1, header=None)
    time_ = data_table.iloc[:, 0]
    angles_ = data_table.iloc[0, 1:]
    ydata_ = data_table.iloc[1:, 1:]
    return data_path_, time_, angles_, ydata_


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

    # TODO: Add date and time to log messages
    # Create initial session
    current_session = Session()

    # # Prompt user for initial measurement data
    # data_path, time, angles, ydata = load_csv_data()

    # Launch Dash app
    app = dash.Dash(external_stylesheets=[dbc.themes.SPACELAB])

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
                style={'width': '99%', 'height': '150px'}
            )
        ], style={'margin-top': '40px', 'margin-left': '10px', 'text-align': 'left'}),

        # PLACEHOLDER: Test button for session log
        dash.html.Div([
            dbc.InputGroup(
                [
                    dbc.Button('Add note to log', id='submit-button', n_clicks=0, color='info'),
                    dbc.Input(id='test-input', value='', type='text')
                ]
            )

        ]),

        # File and session control
        dash.html.H3("File and session controls", className='dash-bootstrap', style={'margin-top': '20px', 'text-align': 'center'}),
        dbc.Container([
            dbc.ButtonGroup([
                dbc.Button('Load session', id='load-session', n_clicks=0, title='Load a previous session in its entirety'),
                dbc.Button('Import from session', id='import-from-session', n_clicks=0, title='Use this to import previous sensors or analysis from another session'),
                dbc.Button('Load data', id='load-data', n_clicks=0, title='Load data from another measurement'),
                dbc.Button('New sensor', id='new-sensor', n_clicks=0),
                dbc.DropdownMenu(
                    label='Choose sensor',
                    color='secondary',
                    children=[
                        dbc.DropdownMenuItem('Sensor 1'),
                        dbc.DropdownMenuItem('Sensor 2')
                    ])
            ])
        ], style={'display': 'flex', 'justify-content': 'center'}),

        # Sensor datatable
        # dash.dash_table.DataTable

    ])

    # PLACEHOLDER: Function for test button
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


    # Connect textarea

    # sys.stdout = output_redirector
    # sys.stderr = output_redirector

    app.run_server(debug=True)
