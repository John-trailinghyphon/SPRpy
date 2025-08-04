# This is the main file where SPRpy is initiated. It is run by executing the file in a Python interpreter.
# The webapp is built using Dash (https://dash.plotly.com/), which is a Python framework for building webapps.

import time
import tomllib
import types
import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from SPRpy_classes import *
from __about__ import version

# Dash app theme (you have to delete browser cookies before changing theme takes effect)
dash_app_theme = dbc.themes.SPACELAB  # Options: CERULEAN, COSMO, CYBORG, DARKLY, FLATLY, JOURNAL, LITERA, LUMEN, LUX,
# MATERIA, MINTY, MORPH, PULSE, QUARTZ, SANDSTONE, SIMPLEX, SKETCHY, SLATE, SOLAR, SPACELAB, SUPERHERO, UNITED, VAPOR, YETI, ZEPHYR.

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

if __name__ == '__main__':

    # Read configuration parameters
    with open('config.toml', 'r') as f:
        config = tomllib.loads(f.read())

    # Access individual parameters as variables
    TIR_range_water_or_long_measurement = config["TIR_range_water_or_long_measurement"]  # TIR range for water --> Automatically used for 50 or more scans per file
    TIR_range_air_or_few_scans = config["TIR_range_air_or_few_scans"]  # TIR range for dry scans --> Automatically used for less than 50 scans per file
    auto_angle_range_points = config["auto_angle_range_points"]  # Number of data points below and above the SPR minimum for auto-detection of the SPR peak
    ask_for_previous_session = config["ask_for_previous_session"]
    default_data_folder = config["default_data_folder"]
    session_host = config["session_host"]
    default_sensor_values = config["default_sensor_values"]
    max_logical_cores = config["max_logical_cores"]

    # Determine how many processes can be used for calculations at a time
    if max_logical_cores == 0:
        logical_cores = multiprocessing.cpu_count()
    elif max_logical_cores > multiprocessing.cpu_count():
        print('Warning: max_logical_cores exceeding system specifications. Using all available cores.')
        logical_cores = multiprocessing.cpu_count()
    else:
        logical_cores = max_logical_cores

    load_session_flag = False
    if ask_for_previous_session is True:

        session_prompt = str(input(
            r'Would you like to load a previous session? Type "y" for yes, or simply skip by pressing enter.'))

        if session_prompt == 'y' or session_prompt == 'Y' or session_prompt == '"y"' or session_prompt == '\'y\'':

            print('Loading previous session, please wait...')
            load_session_flag = True
            session_file = select_file(r'Choose a previous session file', prompt_folder=os.getcwd())

            with open(session_file, 'rb') as file:
                current_session = pickle.load(file)

            # Assert that the session file is of the same version as the current SPRpy version
            if current_session.version != version:
                print('WARNING: The session file was created with a different version of SPRpy than the current one. '
                      'This WILL cause errors for differing X and MAY cause errors for differing Y in version X.Y.Z'
                      'The session file version was ' + current_session.version + ' and the '
                      'currently used SPRpy version is ' + version + '.')
                print('In case of errors, consider pip installing the version of SPRpy (python -m pip install SPRpy==' + current_session.version + ') that was used to create the session file in a separate virtual environment and run SPRpy from there instead.')

                # Compatibility fix for loading older sessions. NOTE: Will still cause erroneous results in result summary tab if multiple layers were fitted for one sensor object
                if current_session.version < '0.3.0' and current_session.fresnel_analysis_instances != {}:
                    for analysis_instance in current_session.fresnel_analysis_instances:
                        current_session.fresnel_analysis_instances[analysis_instance].fitted_layer_index = copy.deepcopy(current_session.fresnel_analysis_instances[analysis_instance].sensor_object.fitted_layer_index)
                        current_session.fresnel_analysis_instances[analysis_instance].fitted_layer = copy.deepcopy(current_session.fresnel_analysis_instances[analysis_instance].sensor_object.optical_parameters.iloc[current_session.fresnel_analysis_instances[analysis_instance].fitted_layer_index[0], 0])

            # Make sure the location and name of the session file is updated
            current_session.location = os.path.dirname(session_file)
            current_session.name = current_session.location.split('/')[-1]

            # Load measurement data
            try:
                current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(
                    path=current_session.current_data_path, default_data_folder=default_data_folder)
            except FileNotFoundError:
                current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(prompt='Select the original data file matching '+current_session.current_data_path)

            # Calculate sensorgram (assume air or liquid medium for TIR calculation based on number of scans)
            if ydata_df.shape[0] > 50:
                TIR_range = TIR_range_water_or_long_measurement
            else:
                TIR_range = TIR_range_air_or_few_scans

            sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, TIR_range, scanspeed)

            # Offset to start at 0 degrees at 0 minutes
            sensorgram_df_selection = sensorgram_df
            sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                                   sensorgram_df_selection['SPR angle'][0]
            sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                                   sensorgram_df_selection['TIR angle'][0]

            # Set current sensor and analysis objects to be the latest one of the session (highest index value)
            current_sensor = current_session.sensor_instances[max(current_session.sensor_instances.keys())]

            try:
                current_fresnel_analysis = current_session.fresnel_analysis_instances[max(current_session.fresnel_analysis_instances.keys())]
            except ValueError:
                current_fresnel_analysis = None

            try:
                current_exclusion_height_analysis = current_session.exclusion_height_analysis_instances[max(current_session.exclusion_height_analysis_instances.keys())]
            except ValueError:
                current_exclusion_height_analysis = None

            # Add note to log
            current_session.log = current_session.log + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + 'Reopened session'

    # If no previous session data was loaded
    if (ask_for_previous_session is False) or (load_session_flag is False):

        # Prompt user for initial measurement data
        print('Please wait...')
        current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(default_data_folder=default_data_folder)

        # Create initial session
        current_session = Session(version, current_data_path=current_data_path)

        # Calculate sensorgram (assume air or liquid medium for TIR calculation based on number of scans)
        if ydata_df.shape[0] > 50:
            TIR_range = TIR_range_water_or_long_measurement
        else:
            TIR_range = TIR_range_air_or_few_scans

        sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, TIR_range, scanspeed)

        # Offset to start at 0 degrees at 0 minutes
        sensorgram_df_selection = sensorgram_df
        sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                               sensorgram_df_selection['SPR angle'][0]
        sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                               sensorgram_df_selection['TIR angle'][0]

        # Add sensor object based on chosen measurement data
        current_sensor = add_sensor_backend(current_session, current_data_path)

        # Calculate TIR angle and update current_sensor.refractive_indices accordingly
        TIR_angle, TIR_fitted_angles, TIR_fitted_ydata = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], TIR_range, scanspeed)
        current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)
        current_sensor.optical_parameters.replace(current_sensor.optical_parameters['n'].iloc[-1], current_sensor.refractive_indices[-1], inplace=True)
        current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
            sensor_number=current_sensor.object_id,
            sensor_name=current_sensor.name,
            channel=current_sensor.channel,
            fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
            fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

        current_session.save_session()
        current_session.save_sensor(current_sensor.object_id)

        # Add empty analysis objects
        current_fresnel_analysis = None
        current_exclusion_height_analysis = None

    # Dash app
    app = dash.Dash(name='SPRpy', title='SPRpy', external_stylesheets=[dash_app_theme])
    app._favicon = 'icon.ico'
    # app.config.suppress_callback_exceptions = True  # NOTE: Comment out this line for debugging purposes if callbacks do not fire when supposed to

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

    d_n_pair_fig = go.Figure(go.Scatter(
                x=[0],
                y=[0],
                mode='lines',
                name='Buffer',
                showlegend=False,
                line_color='#636EFA'
            ))
    d_n_pair_fig.update_layout(
        xaxis_title=r'$\large{\text{Refractive index}}$',
        yaxis_title=r'$\large{\text{Height [nm]}}$',
        font_family='Balto',
        font_size=19,
        margin_r=25,
        margin_l=60,
        margin_t=40,
        template='simple_white',
        uirevision=True)
    d_n_pair_fig.update_xaxes(mirror=True, showline=True)
    d_n_pair_fig.update_yaxes(mirror=True, showline=True)

    try:
        x_barplot = [[current_session.fresnel_analysis_instances[
                          fresnel_inst].fitted_layer for
                      fresnel_inst in current_session.fresnel_analysis_instances],
                     [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in
                      current_session.fresnel_analysis_instances]]
        y_barplot = [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for
                     fresnel_inst in current_session.fresnel_analysis_instances]
        result_barplot_fig = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
        result_barplot_fig.update_layout(
            yaxis_title='Fitted value',
            font_family='Balto',
            font_size=19,
            margin_r=25,
            margin_l=60,
            margin_t=40,
            template='simple_white',
            uirevision=True,
            height=600)
        result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
        result_barplot_fig.update_yaxes(mirror=True, showline=True)
    except:
        result_barplot_fig = go.Figure(go.Bar(x=[0], y=[0]))
        result_barplot_fig.update_layout(
            yaxis_title='Fitted value',
            font_family='Balto',
            font_size=19,
            margin_r=25,
            margin_l=60,
            margin_t=40,
            template='simple_white',
            uirevision=True,
            height=600)
        result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
        result_barplot_fig.update_yaxes(mirror=True, showline=True)

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
                dash.html.Div([
                    dash.dcc.Markdown('''
                # **#SPRpy#**
                ''', className='dash-bootstrap'),
                    dash.dcc.Markdown('''
                    #### **v ''' + version + '''**
                    ''', className='dash-bootstrap', style={'display': 'flex', 'justify-content': 'center', 'margin-right': '10px'}),
                ], style={'margin-top': '6rem', 'margin-left': '5rem', 'margin-right': '5rem'}),
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
        # TODO: Add an Interval component that updates the session log once per minute (when/if starting to add automatic log messages)
        # Session log div
        dash.html.Div([
            dash.html.H3('{name_} - Session log'.format(name_=current_session.name),
                         className='dash-bootstrap',
                         id='session-title'),
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
                    dbc.Button('Rename session', id='rename-session-button', n_clicks=0, color='warning'),
                    dbc.Input(id='test-input', value='', type='text', style={'margin-right': '2%'})
                ]
            )

        ], style={'margin-left': '2%'}),

        # File and session control
        dash.html.H3("File and sensor controls", className='dash-bootstrap', style={'margin-top': '20px', 'text-align': 'center'}),
        dash.html.Div(['Current measurement file:    ', current_data_path.split('/')[-1]],
                      id='datapath-textfield',
                      style={'margin-right': '10px', 'textAlign': 'center'}),
        dbc.Container([
            dbc.ButtonGroup([
                dbc.Button('Load new data',
                           id='load-data',
                           n_clicks=0,
                           color='primary',
                           title='Load data from another measurement. Analysis is always performed on this active measurement'),
                dash.dcc.Store(id='loaded-new-measurement', storage_type='memory'),
                # TODO: Add functionality for this button
                # dbc.Button('Import result',
                #            id='import-from-session',
                #            n_clicks=0,
                #            color='primary',
                #            title='Use this to import previous results from another session'),
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
                dbc.Button('Remove current sensor',
                           id='remove-sensor-button',
                           n_clicks=0,
                           color='primary',
                           title='Removes the currently selected sensor from the session.'),
                dbc.Modal([
                    dbc.ModalHeader(dbc.ModalTitle('Removing sensor object')),
                    dbc.ModalBody('Are you sure you want to delete the currently selected sensor?\n(at least one remaining required)'),
                    dbc.ModalFooter(
                        dbc.ButtonGroup([
                            dbc.Button('Confirm', id='remove-sensor-confirm',
                                       color='success',
                                       n_clicks=0),
                            dbc.Button('Cancel', id='remove-sensor-cancel',
                                       color='danger',
                                       n_clicks=0)
                        ])
                    )
                ],
                    id='remove-sensor-modal',
                    size='sm',
                    is_open=False,
                    backdrop='static',
                    keyboard=False),
                dbc.DropdownMenu(
                    id='chosen-sensor-dropdown',
                    label='Sensors',
                    color='primary',
                    children=[
                        dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name, id={'type': 'sensor-list', 'index': sensor_id},
                                             n_clicks=0) for sensor_id in current_session.sensor_instances], style={'margin-left': '-5px'})
            ])
        ], style={'margin-bottom': '20px', 'display': 'flex', 'justify-content': 'center'}),

        # Sensor datatable
        dash.html.Div([
            dash.html.Div([
                dash.html.H4(['S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                    sensor_number=current_sensor.object_id,
                    sensor_name=current_sensor.name,
                    channel=current_sensor.channel,
                    fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                    fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])
                ], id='sensor-table-title', style={'text-align': 'center'}),
                dash.html.Div([
                    dash.dash_table.DataTable(data=current_sensor.optical_parameters.to_dict('records'),
                                              columns=[{'name': 'Layers', 'id': 'Layers', 'type': 'text'},
                                                       {'name': 'd [nm]', 'id': 'd [nm]', 'type': 'numeric'},
                                                       {'name': 'n', 'id': 'n', 'type': 'numeric'},
                                                       {'name': 'k', 'id': 'k', 'type': 'numeric'}],
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
                    dbc.Button('Save edited values',
                               id='table-update-values',
                               n_clicks=0,
                               color='danger',
                               title='Save the displayed values to sensor after editing'),
                    dbc.Button('Select variable to fit',
                               id='table-select-fitted',
                               n_clicks=0,
                               color='success',
                               title='Click this button after selecting a different parameter to fit by clicking it such'
                                     ' that it is marked in red. NOTE: First click "Save edited values" if new layers were added.'),
                    dbc.Button('Rename sensor',
                               id='rename-sensor-button',
                               n_clicks=0,
                               color='warning',
                               title='Rename the current sensor'),
                    dbc.Modal([
                        dbc.ModalHeader(dbc.ModalTitle('Rename sensor')),
                        dbc.ModalBody(dbc.Input(id='rename-sensor-input', placeholder='Give a name...', type='text')),
                        dbc.ModalFooter(dbc.Button('Confirm', id='rename-sensor-confirm', color='success', n_clicks=0))
                    ],
                        id='rename-sensor-modal',
                        size='sm',
                        is_open=False,
                        backdrop='static',
                        keyboard=False),
                    dbc.Button(
                        "Show default values",
                        id="show-default-param-button",
                        color="secondary",
                        n_clicks=0,
                        title='CTRL+Z not supported for table. Check default values here if needed.'
                    ),
                ], style={'width': '672px', 'margin-left': '4px', 'margin-top': '5px', 'margin-bottom': '20px'}),
            ], style={'width': '675px'}),
            dash.html.Div([
                dbc.Collapse(
                    dbc.Card(
                        dbc.CardBody(
                            dbc.Table.from_dataframe(pd.DataFrame(
                                default_sensor_values,
                            ), size='sm', striped=True, bordered=True, hover=True)
                        ), style={'width': '800px'}),
                    id='default-values-collapse',
                    is_open=False)
            ], style={'margin-top': '40px', 'margin-left': '10px'}),
        ], style={'display': 'flex', 'justify-content': 'center'}),

        # Analysis tabs
        dash.html.Div([
            dash.html.H1(['Analysis options']),
            dbc.Tabs([

                # Response quantification tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.dcc.Graph(id='quantification-reflectivity-graph',
                                           figure=reflectivity_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.Button('Add data trace',
                                           id='quantification-reflectivity-add-data-trace',
                                           n_clicks=0,
                                           color='danger',
                                           title='Add a measurement trace to the figure from an external dry scan .csv file. The most recent scan in the file is used.'),
                                dbc.Button('Add fresnel trace',
                                           id='quantification-reflectivity-add-fresnel-trace',
                                           n_clicks=0,
                                           color='success',
                                           title='Add a fresnel calculation trace to the figure based on current sensor values.'),
                                dbc.Button('Clear traces',
                                           id='quantification-reflectivity-clear-traces',
                                           n_clicks=0,
                                           color='warning',
                                           title='Clear added traces (required to regain sensorgram hover data selection).'),
                                dbc.DropdownMenu(
                                    id='reflectivity-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[
                                        dbc.DropdownMenuItem('.PNG', id='quantification-reflectivity-save-png', n_clicks=0),
                                        dbc.DropdownMenuItem('.SVG', id='quantification-reflectivity-save-svg', n_clicks=0),
                                        dbc.DropdownMenuItem('.HTML', id='quantification-reflectivity-save-html', n_clicks=0),
                                        dbc.DropdownMenuItem('.csv', id='quantification-reflectivity-save-csv', n_clicks=0)],
                                )
                            ], style={'margin-left': '13%'}),
                        ], style={'width': '35%'}),
                        # TODO: Add option to plot our PureKinetics stuff in this graph window
                        dash.html.Div([
                            dash.dcc.Graph(id='quantification-sensorgram-graph',
                                           figure=sensorgram_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.Switch(
                                    id='hover-selection-switch',
                                    label='Lock hover selection',
                                    value=False),
                                dbc.DropdownMenu(
                                    id='sensorgram-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[dbc.DropdownMenuItem('.PNG', id='quantification-sensorgram-save-png', n_clicks=0),
                                              dbc.DropdownMenuItem('.SVG', id='quantification-sensorgram-save-svg', n_clicks=0),
                                              dbc.DropdownMenuItem('.HTML', id='quantification-sensorgram-save-html', n_clicks=0),
                                              dbc.DropdownMenuItem('.csv', id='quantification-sensorgram-save-csv', n_clicks=0)],
                                    style={'margin-left': '100%'}),
                            ], style={'margin-left': '27.5%'}),
                        ], style={'width': '60%'})
                    ], id='quantification-tab-content', style={'display': 'flex', 'justify-content': 'center'})
                ], label='Response quantification', tab_id='quantification-tab', style={'margin-top': '10px'}),

                # Fresnel modelling tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.html.H3(['Settings']),
                            dbc.Form([
                                dash.html.Div([
                                    dbc.ButtonGroup([
                                        dbc.Button('Add new fresnel analysis',
                                                   id='add-fresnel-analysis-button',
                                                   n_clicks=0,
                                                   color='primary',
                                                   title='Add a new fresnel analysis object for the current sensor.'),
                                        dash.dcc.Store(id='add-fresnel-analysis-signal', storage_type='session'),
                                        dbc.Modal([
                                            dbc.ModalHeader(dbc.ModalTitle('New fresnel analysis object')),
                                            dbc.ModalBody(
                                                dbc.Input(id='fresnel-analysis-name-input', placeholder='Give a name...', type='text')),
                                            dbc.ModalFooter(
                                                dbc.Button('Confirm', id='add-fresnel-analysis-confirm', color='success',
                                                           n_clicks=0))
                                        ],
                                            id='add-fresnel-analysis-modal',
                                            size='sm',
                                            is_open=False,
                                            backdrop='static',
                                            keyboard=False),
                                        dbc.DropdownMenu(id='fresnel-analysis-dropdown',
                                                         label='Choose analysis',
                                                         color='primary',
                                                         children=[dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                                                                        id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                                                                                        n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]),
                                        dbc.Button('Rename analysis',
                                                   id='rename-fresnel-analysis-button',
                                                   n_clicks=0,
                                                   color='warning',
                                                   title='Rename the current fresnel analysis object.'),
                                        dbc.Modal([
                                            dbc.ModalHeader(dbc.ModalTitle('Rename fresnel analysis')),
                                            dbc.ModalBody(
                                                dbc.Input(id='rename-fresnel-analysis-input', placeholder='Give a name...',
                                                          type='text')),
                                            dbc.ModalFooter(
                                                dbc.Button('Confirm', id='rename-fresnel-analysis-confirm', color='success',
                                                           n_clicks=0))
                                        ],
                                            id='rename-fresnel-analysis-modal',
                                            size='sm',
                                            is_open=False,
                                            backdrop='static',
                                            keyboard=False),
                                        dbc.Button('Remove analysis',
                                                   id='remove-fresnel-analysis-button',
                                                   n_clicks=0,
                                                   color='primary',
                                                   title='Remove the currently selected analysis.'),
                                        dbc.Modal([
                                            dbc.ModalHeader(dbc.ModalTitle('Removing fresnel analysis object')),
                                            dbc.ModalBody('Are you sure you want to delete the currently selected analysis?'),
                                            dbc.ModalFooter(
                                                dbc.ButtonGroup([
                                                    dbc.Button('Confirm', id='remove-fresnel-analysis-confirm',
                                                               color='success',
                                                               n_clicks=0),
                                                    dbc.Button('Cancel', id='remove-fresnel-analysis-cancel',
                                                               color='danger',
                                                               n_clicks=0)
                                                ])
                                            )
                                        ],
                                            id='remove-fresnel-analysis-modal',
                                            size='sm',
                                            is_open=False,
                                            backdrop='static',
                                            keyboard=False),
                                        dbc.Button('Batch analysis',
                                                   id='batch-fresnel-analysis-button',
                                                   n_clicks=0,
                                                   color='primary',
                                                   title='Perform automatic batch fresnel modelling on several similar measurement files based on a selected example sensor and example analysis.'),
                                        dbc.Modal([
                                            dbc.ModalHeader(dbc.ModalTitle('Start automatic batch fresnel modelling')),
                                            dbc.ModalBody([dash.html.Div(['Prerequisites:']),
                                                           dash.html.Div([' - All files must be in the same folder.']),
                                                           dash.html.Div([' - All files must use the same sensor layer structure and wavelength.']),
                                                           dbc.Button('Choose measurement files',
                                                                      id='batch-fresnel-analysis-choose-files',
                                                                      n_clicks=0,
                                                                      style={'margin-bottom': '20px'}),
                                                           dash.dcc.Store(id='batch-fresnel-analysis-files', storage_type='session'),
                                                           dbc.Row([
                                                               dash.dcc.Dropdown(
                                                                   id='batch-fresnel-analysis-example-sensor-dropdown',
                                                                   placeholder='Select example sensor',
                                                                   options=[{'label': 'S' + str(sensor_id) + ' ' +
                                                                                      current_session.sensor_instances[
                                                                                          sensor_id].name,
                                                                             'value': sensor_id} for sensor_id in
                                                                            current_session.sensor_instances]),
                                                               dash.dcc.Dropdown(
                                                                   id='batch-fresnel-analysis-example-analysis-dropdown',
                                                                   placeholder='Select example analysis',
                                                                   options=[{'label': 'FM' + str(
                                                                       fresnel_analysis_id) + ' ' +
                                                                                      current_session.fresnel_analysis_instances[
                                                                                          fresnel_analysis_id].name,
                                                                             'value': fresnel_analysis_id} for
                                                                            fresnel_analysis_id in
                                                                            current_session.fresnel_analysis_instances])
                                                           ], id='batch-fresnel-analysis-example-row'),
                                                           dash.dcc.RadioItems(options=[{'label': 'Use copies of the example background sensor', 'value': 0},
                                                                                        {'label': 'Match several background sensors individually', 'value': 1}],
                                                                               value=0,
                                                                               id='batch-fresnel-analysis-radio-selection',
                                                                               style={'margin-bottom': '30px'}),
                                                           dash.dcc.RadioItems(options=[{'label': 'Add layer to sensor backgrounds directly', 'value': 0},
                                                                                        {'label': 'Use new copy of sensor backgrounds, then add layer', 'value': 1}],
                                                                               value=0,
                                                                               id='batch-fresnel-analysis-newlayer-radio-selection',
                                                                               style={'visibility': 'hidden'}),
                                                           dbc.Row([
                                                               dash.dcc.Dropdown(
                                                                   id='batch-fresnel-analysis-background-sensors-dropdown',
                                                                   placeholder='Select individual background sensors',
                                                                   multi=True,
                                                                   clearable=True,
                                                                   options=[{'label': 'S' + str(sensor_id) + ' ' +
                                                                                      current_session.sensor_instances[
                                                                                          sensor_id].name,
                                                                             'value': sensor_id} for sensor_id in
                                                                            current_session.sensor_instances])
                                                           ], id='batch-fresnel-analysis-background-sensors-dropdown-row', style={'visibility': 'hidden'}),
                                                           dbc.Button('Submit backgrounds',
                                                                      id='batch-fresnel-analysis-background-sensors-button-submit',
                                                                      n_clicks=0,
                                                                      style={'visibility': 'hidden'}),
                                                           dbc.Table.from_dataframe(pd.DataFrame({'Measurement file': [''], 'Matching background sensor': ['']}), bordered=True, id='batch-fresnel-analysis-table', style={'visibility': 'hidden', 'margin-top': '20px'}),
                                                           dash.dcc.Store(
                                                               id='batch-fresnel-analysis-background-sensors',
                                                               storage_type='session'),
                                                           dash.dcc.Store(id='batch-fresnel-analysis-start',
                                                                          storage_type='memory'),
                                                           dash.dcc.Store(id='batch-fresnel-analysis-done',
                                                                          storage_type='memory'),
                                                           dash.dcc.Store(id='batch-fresnel-analysis-finish',
                                                                          storage_type='memory'),
                                                           ]),
                                            dbc.ModalFooter([
                                                dbc.ButtonGroup([
                                                    dbc.Spinner(color='success', type='border',
                                                                id='batch-fresnel-spinner',
                                                                spinner_style={'visibility': 'hidden',
                                                                               'margin-top': '10px',
                                                                               'margin-right': '10px',
                                                                               'width': '2rem',
                                                                               'height': '2rem'}),
                                                    dbc.Button('Start', id='batch-fresnel-analysis-confirm',
                                                               color='success',
                                                               n_clicks=0),
                                                    dbc.Button('Cancel', id='batch-fresnel-analysis-cancel',
                                                               color='danger',
                                                               n_clicks=0)
                                                ])
                                            ])
                                        ],
                                            id='batch-fresnel-analysis-modal',
                                            size='xl',
                                            is_open=False,
                                            backdrop='static',
                                            keyboard=False),
                                    ])
                                ]),
                                dash.html.Div([
                                    dbc.Collapse(
                                        dbc.Card(
                                            dbc.CardBody(
                                                dbc.Form([
                                                    dbc.Row([
                                                        dbc.Label(
                                                            'Data path: \n' + current_data_path.split('/')[-1],
                                                            id='fresnel-fit-datapath')
                                                    ], style={'margin-bottom': '10px'}),
                                                    dbc.Row([
                                                        dbc.Label(
                                                            'Sensor: S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                                                                sensor_number=current_sensor.object_id,
                                                                sensor_name=current_sensor.name,
                                                                channel=current_sensor.channel,
                                                                fitted_layer=current_sensor.optical_parameters.iloc[
                                                                    current_sensor.fitted_layer_index[0], 0],
                                                                fitted_param=current_sensor.optical_parameters.columns[
                                                                    current_sensor.fitted_layer_index[1]]),
                                                            id='fresnel-fit-sensor')
                                                    ], style={'margin-bottom': '10px'}),
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
                                                                          value=float(current_sensor.fitted_var) - float(current_sensor.fitted_var) / 2,
                                                                          type='number'),
                                                                dbc.Input(id='fresnel-fit-option-upperbound',
                                                                          value=float(current_sensor.fitted_var) + float(current_sensor.fitted_var) / 2,
                                                                          type='number')
                                                            ])
                                                        ], width=4),
                                                        dbc.Label('p-factor', width='auto'),
                                                        dbc.Col([
                                                            dbc.Input(id='fresnel-fit-option-pfactor',
                                                                      value=1.0, type='number')
                                                        ], width=2),
                                                    ], style={'margin-bottom': '10px'}),
                                                    dbc.Row([
                                                        dbc.Label('Angle range', width='auto'),
                                                        dbc.Col([
                                                            dash.dcc.RangeSlider(value=[reflectivity_df['angles'].iloc[reflectivity_df['ydata'].idxmin()-auto_angle_range_points[0]], reflectivity_df['angles'].iloc[reflectivity_df['ydata'].idxmin()+auto_angle_range_points[1]]],
                                                                                 min=reflectivity_df['angles'].iloc[0],
                                                                                 max=reflectivity_df['angles'].iloc[-1],
                                                                                 marks={mark_ind: str(mark_ind) for mark_ind in range(reflectivity_df['angles'].iloc[0].astype('int'), reflectivity_df['angles'].iloc[-1].astype('int')+1, 1)},
                                                                                 step=0.005,
                                                                                 allowCross=False,
                                                                                 tooltip={"placement": "top",
                                                                                          "always_visible": True},
                                                                                 id='fresnel-fit-option-rangeslider')
                                                        ])
                                                    ], style={'margin-bottom': '10px'}),
                                                    dbc.Row([
                                                        dbc.Col([
                                                            dbc.Checkbox(
                                                                id='fresnel-analysis-offset-fit',
                                                                label="Fit reflectivity offset?",
                                                                value=True)
                                                        ], width='auto', style={'padding-top': '20px'}),
                                                        dbc.Col([
                                                            dbc.Checkbox(
                                                                id='fresnel-analysis-elastomer-fit',
                                                                label="Fit prism k?",
                                                                disabled=False,
                                                                value=True)
                                                        ], width='auto', style={'padding-top': '20px'}),
                                                    ]),
                                                    dbc.Row([
                                                        dbc.Label('Manual prism k correction (iterative) [1e-3]', width='auto'),
                                                        dbc.Col([
                                                            dash.dcc.Slider(min=-0.0005, max=0.0005,
                                                                            step=0.00005,
                                                                            marks={-0.0005: '-5', -0.0004: '-4',
                                                                                   -0.0003: '-3', -0.0002: '-2',
                                                                                   -0.0001: '-1', 0: '0',
                                                                                   0.0001: '1', 0.0002: '2',
                                                                                   0.0003: '3', 0.0004: '4',
                                                                                   0.0005: '5'},
                                                                            tooltip={"placement": "top"},
                                                                            id='fresnel-fit-option-extinctionslider')
                                                        ], style={'margin-top': '15px'})
                                                    ], id='fresnel-analysis-manual-prism-k-row', style={'margin-bottom': '10px', 'visibility': 'hidden'}),
                                                    dbc.Row([
                                                        dbc.Label('Fit result: ', id='fresnel-fit-result')
                                                    ], style={'margin-bottom': '10px'})
                                                ])
                                            )
                                        ), id='fresnel-analysis-option-collapse', is_open=False)
                                ])
                            ], id='fresnel-fit-options-form')
                        ], style={'margin-top': '1.9rem', 'width': '65%'}),
                        dash.html.Div([
                            dash.dcc.Graph(id='fresnel-reflectivity-graph',
                                           figure=reflectivity_fig,
                                           mathjax=True),
                            dbc.ButtonGroup([
                                dbc.Button('Run modelling',
                                           id='fresnel-reflectivity-run-model',
                                           n_clicks=0,
                                           color='success',
                                           title='Run the fresnel model',
                                           disabled=False),
                                dash.dcc.Store(id='fresnel-reflectivity-run-finished', storage_type='session'),
                                dbc.DropdownMenu(
                                    id='fresnel-save-dropdown',
                                    label='Save as...',
                                    color='info',
                                    children=[
                                        dbc.DropdownMenuItem('.PNG', id='fresnel-reflectivity-save-png', n_clicks=0),
                                        dbc.DropdownMenuItem('.SVG', id='fresnel-reflectivity-save-svg', n_clicks=0),
                                        dbc.DropdownMenuItem('.HTML', id='fresnel-reflectivity-save-html', n_clicks=0),
                                        dbc.DropdownMenuItem('.csv', id='fresnel-reflectivity-save-csv', n_clicks=0)],
                                    style={'margin-left': '-5px'})
                            ], style={'margin-left': '30%'}),
                        ], style={'width': '35%', 'margin-top': '1.9rem', 'margin-left': '5%'}),
                    ], id='fresnel-tab-content', style={'display': 'flex', 'justify-content': 'center'})
                ], label='Fresnel modelling', tab_id='fresnel-tab', style={'margin-top': '10px'}),

                # Exclusion height determination tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.html.Div([
                                dash.html.H3(['Settings']),
                                dbc.Form([
                                    dash.html.Div([
                                        dbc.ButtonGroup([
                                            dbc.Button('Add new exclusion height analysis',
                                                       id='add-exclusion-height-analysis-button',
                                                       n_clicks=0,
                                                       color='primary',
                                                       title='Add a new exclusion analysis object for the current sensor.'
                                                       ),
                                            dash.dcc.Store(id='add-exclusion-height-analysis-signal', storage_type='session'),
                                            dbc.Modal([
                                                dbc.ModalHeader(dbc.ModalTitle('New exclusion height analysis object')),
                                                dbc.ModalBody([
                                                    dash.dcc.Dropdown(id='exclusion-choose-background-dropdown',
                                                                      placeholder='Choose background...',
                                                                      options=[{'label': 'FM' + str(fresnel_id) + ' ' +
                                                                                         current_session.fresnel_analysis_instances[
                                                                                             fresnel_id].name,
                                                                                'value': fresnel_id} for fresnel_id in
                                                                               current_session.fresnel_analysis_instances]),
                                                    dbc.Input(id='exclusion-height-analysis-name-input',
                                                              placeholder='Give a name...', type='text')
                                                ]),
                                                dbc.ModalFooter(
                                                    dbc.Button('Confirm', id='add-exclusion-height-analysis-confirm',
                                                               color='success',
                                                               n_clicks=0)
                                                )
                                            ],
                                                id='add-exclusion-height-analysis-modal',
                                                size='sm',
                                                is_open=False,
                                                backdrop='static',
                                                keyboard=False),
                                            dbc.DropdownMenu(id='exclusion-height-analysis-dropdown',
                                                             label='Choose analysis',
                                                             color='primary',
                                                             children=[dbc.DropdownMenuItem(
                                                                 'EH' + str(exclusion_id) + ' ' +
                                                                 current_session.exclusion_height_analysis_instances[
                                                                     exclusion_id].name,
                                                                 id={'type': 'exclusion-analysis-list',
                                                                     'index': exclusion_id},
                                                                 n_clicks=0) for exclusion_id in
                                                                 current_session.exclusion_height_analysis_instances]),
                                            dbc.Button('Remove analysis',
                                                       id='remove-exclusion-height-analysis-button',
                                                       n_clicks=0,
                                                       color='primary',
                                                       title='Remove the currently selected analysis.'),
                                            dbc.Modal([
                                                dbc.ModalHeader(dbc.ModalTitle('Removing exclusion-height analysis object')),
                                                dbc.ModalBody(
                                                    'Are you sure you want to delete the currently selected analysis?'),
                                                dbc.ModalFooter(
                                                    dbc.ButtonGroup([
                                                        dbc.Button('Confirm', id='remove-exclusion-height-analysis-confirm',
                                                                   color='success',
                                                                   n_clicks=0),
                                                        dbc.Button('Cancel', id='remove-exclusion-height-analysis-cancel',
                                                                   color='danger',
                                                                   n_clicks=0)
                                                    ])
                                                )
                                            ],
                                                id='remove-exclusion-height-analysis-modal',
                                                size='sm',
                                                is_open=False,
                                                backdrop='static',
                                                keyboard=False),
                                        ])
                                    ]),
                                    dash.html.Div([
                                        dbc.Collapse(
                                            dbc.Card(
                                                dbc.CardBody(
                                                    dbc.Form([
                                                        dbc.Row([
                                                            dbc.Label('Selected analysis: ',
                                                                      id='exclusion-height-analysis-label')
                                                        ], style={'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Label(
                                                                'Sensor: S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                                                                    sensor_number=current_sensor.object_id,
                                                                    sensor_name=current_sensor.name,
                                                                    channel=current_sensor.channel,
                                                                    fitted_layer=current_sensor.optical_parameters.iloc[
                                                                        current_sensor.fitted_layer_index[0], 0],
                                                                    fitted_param=
                                                                    current_sensor.optical_parameters.columns[
                                                                        current_sensor.fitted_layer_index[1]]),
                                                                id='exclusion-height-sensor-label')
                                                        ], style={'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Label(
                                                                'Fresnel analysis: FM{analysis_number} {analysis_name}'.format(
                                                                    analysis_number=1,
                                                                    analysis_name='Placeholder'),
                                                                id='exclusion-height-fresnel-analysis-label')
                                                        ], style={'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Label('Height bounds (min, max)', width='auto'),
                                                            dbc.Col([
                                                                dbc.InputGroup([
                                                                    dbc.Input(id='exclusion-height-option-lowerbound',
                                                                              value=float(0),
                                                                              type='number'),
                                                                    dbc.Input(id='exclusion-height-option-upperbound',
                                                                              value=float(200),
                                                                              type='number')
                                                                ])
                                                            ], width=7)
                                                        ], style={'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Label('Resolution', width='auto'),
                                                            dbc.Col([
                                                                dbc.InputGroup([
                                                                    dbc.Input(id='exclusion-height-option-resolution',
                                                                              value=int(100),
                                                                              type='number'),
                                                                ])
                                                            ], width=3),
                                                        ], style={'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Col([
                                                                dbc.Checkbox(
                                                                    id='exclusion-height-analysis-offset-refit',
                                                                    label="Refitting offset?",
                                                                    value=False)
                                                            ], width='auto'),
                                                            dbc.Col([
                                                                dbc.Checkbox(
                                                                    id='exclusion-height-analysis-prism-refit',
                                                                    label="Refitting prism k?",
                                                                    value=False)
                                                            ], width='auto'),
                                                        ]),
                                                        dbc.Row([
                                                            dbc.Label('Injection points: ', width='auto', id='exclusion-height-settings-injection-points')
                                                        ]),
                                                        dbc.Row([
                                                            dbc.Label('Buffer points: ', width='auto', id='exclusion-height-settings-buffer-points')
                                                        ]),
                                                        dbc.Row([
                                                            dbc.Label('Probe points: ', width='auto', id='exclusion-height-settings-probe-points')
                                                        ]),
                                                        dbc.Row([
                                                            dbc.Col([
                                                                dbc.Button('Initialize model',
                                                                           id='exclusion-height-initialize-model',
                                                                           color='primary',
                                                                           n_clicks=0,
                                                                           size='lg',
                                                                           title='Prepare model after selecting all points.')
                                                            ], width=6)
                                                        ])
                                                    ])
                                                )
                                            ), id='exclusion-height-analysis-option-collapse', is_open=False)
                                    ])
                                ], id='exclusion-height-fit-options-form'),
                                dbc.Collapse([
                                    dash.html.Div([
                                        dbc.ButtonGroup([
                                            dbc.Spinner(color='success', type='border', id='exclusion-height-spinner', spinner_style={'visibility': 'hidden', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}),
                                            dbc.Button('Run full calculation', id='exclusion-height-run-button',
                                                       color='success',
                                                       n_clicks=0,
                                                       size='lg',
                                                       disabled=False),
                                            dbc.Button('Abort', id='exclusion-height-abort-button',
                                                       color='danger',
                                                       n_clicks=0,
                                                       size='lg',
                                                       disabled=True,
                                                       title='Cancelling a running calculation. NOTE THAT PREVIOUS PROGRESS IS STILL OVERWRITTEN.'),
                                        ]),
                                        dash.dcc.Store(id='exclusion-run-finished', storage_type='session')
                                    ])
                                ], id='exclusion-height-progress-collapse', is_open=False, style={'margin-top': '120px'})
                            ], style={'margin-top': '80px'}),
                            dbc.Collapse([
                                dash.html.Div([
                                    dash.dcc.Graph(id='exclusion-height-sensorgram-graph',
                                                   figure=sensorgram_fig,
                                                   mathjax=True),
                                    dash.html.Div([
                                        dbc.Label('Click-action selector', style={'margin-left': '5%', 'margin-top': '35px'}),
                                        dbc.RadioItems(
                                            options=[
                                                {"label": "Offset data", "value": 1},
                                                {"label": "Choose injection points", "value": 2},
                                                {"label": "Choose buffer points", "value": 3},
                                                {"label": "Choose probe points", "value": 4}],
                                            value=1,
                                            id='exclusion-height-click-action-selector',
                                            style={'margin-left': '20px'}),
                                        dbc.Button('Clear selected points', id='exclusion-height-click-action-clear',
                                                   color='warning',
                                                   n_clicks=0,
                                                   style={'margin-left': '20px', 'margin-top': '35px', 'margin-bot': '35px', 'margin-right': '18%', 'line-height': '1.5'}),
                                        dbc.DropdownMenu(
                                            id='exclusion-height-sensorgram-save-dropdown',
                                            label='Save as...',
                                            color='info',
                                            children=[dbc.DropdownMenuItem('.PNG', id='exclusion-height-sensorgram-save-png',
                                                                           n_clicks=0),
                                                      dbc.DropdownMenuItem('.SVG', id='exclusion-height-sensorgram-save-svg',
                                                                           n_clicks=0),
                                                      dbc.DropdownMenuItem('.HTML', id='exclusion-height-sensorgram-save-html',
                                                                           n_clicks=0),
                                                      dbc.DropdownMenuItem('.csv', id='exclusion-height-sensorgram-save-csv', n_clicks=0)])
                                    ], style={'display': 'flex', 'justify-content': 'left'}),
                                ])
                            ], id='exclusion-height-sensorgram-collapse', is_open=False, style={'width': '60%', 'margin-left': '3%'})
                        ], style={'display': 'flex', 'justify-content': 'center'}),

                        # Results
                        dbc.Collapse([
                            dash.html.Div([
                                dash.html.H3(['Exclusion height results'], style={'display': 'flex', 'justify-content': 'left'}),
                                dash.html.Div([
                                    dbc.Label('Mean exclusion height: None',
                                              id='exclusion-height-result-mean-height')
                                ], style={'margin-top': '30px', 'display': 'flex', 'justify-content': 'left'}),
                                dash.html.Div([
                                    dbc.Label('Mean exclusion RI: None',
                                              id='exclusion-height-result-mean-RI')
                                ], style={'margin-top': '10px', 'display': 'flex', 'justify-content': 'left'}),
                                dash.html.Div([
                                    dbc.Label('All exclusion heights: None',
                                              id='exclusion-height-result-all-heights')
                                ], style={'margin-top': '10px', 'display': 'flex', 'justify-content': 'left'}),
                                dash.html.Div([
                                    dbc.Label('All exclusion RIs: None',
                                              id='exclusion-height-result-all-RI')
                                ], style={'margin-top': '10px', 'display': 'flex', 'justify-content': 'left'}),
                                dbc.Label('Injection step', id='exclusion-height-result-pagination-label',
                                          style={'display': 'flex', 'justify-content': 'center'}),
                                dash.html.Div([
                                    dbc.Pagination(max_value=2, previous_next=True, id='exclusion-height-result-pagination')
                                ], style={'display': 'flex', 'justify-content': 'center'}),
                                dash.html.Div([
                                    dash.html.Div([
                                        dash.dcc.Graph(id='exclusion-height-SPRvsTIR-graph',
                                                       figure=reflectivity_fig,
                                                       mathjax=True),
                                        dbc.ButtonGroup([
                                            dbc.DropdownMenu(
                                                id='exclusion-height-SPRvsTIR-save-dropdown',
                                                label='Save as...',
                                                color='info',
                                                children=[
                                                    dbc.DropdownMenuItem('.PNG',
                                                                         id='exclusion-height-SPRvsTIR-save-png',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.SVG',
                                                                         id='exclusion-height-SPRvsTIR-save-svg',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.HTML',
                                                                         id='exclusion-height-SPRvsTIR-save-html',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.csv',
                                                                         id='exclusion-height-SPRvsTIR-save-csv',
                                                                         n_clicks=0)])
                                        ], style={'margin-left': '13%'}),
                                    ], style={'width': '33%'}),
                                    dash.html.Div([
                                        dash.dcc.Graph(id='exclusion-height-reflectivity-graph',
                                                       figure=reflectivity_fig,
                                                       mathjax=True),
                                        dbc.ButtonGroup([
                                            dbc.DropdownMenu(
                                                id='exclusion-height-reflectivity-save-dropdown',
                                                label='Save as...',
                                                color='info',
                                                children=[
                                                    dbc.DropdownMenuItem('.PNG', id='exclusion-height-reflectivity-save-png',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.SVG', id='exclusion-height-reflectivity-save-svg',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.HTML', id='exclusion-height-reflectivity-save-html',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.csv',
                                                                         id='exclusion-height-reflectivity-save-csv',
                                                                         n_clicks=0)])
                                        ], style={'margin-left': '13%'}),
                                    ], style={'width': '33%'}),
                                    dash.html.Div([
                                        dash.dcc.Graph(id='exclusion-height-d-n-pair-graph',
                                                       figure=d_n_pair_fig,
                                                       mathjax=True),
                                        dbc.ButtonGroup([
                                            dbc.DropdownMenu(
                                                id='exclusion-height-d-n-pair-save-dropdown',
                                                label='Save as...',
                                                color='info',
                                                children=[
                                                    dbc.DropdownMenuItem('.PNG',
                                                                         id='exclusion-height-d-n-pair-save-png',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.SVG',
                                                                         id='exclusion-height-d-n-pair-save-svg',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.HTML',
                                                                         id='exclusion-height-d-n-pair-save-html',
                                                                         n_clicks=0),
                                                    dbc.DropdownMenuItem('.csv',
                                                                         id='exclusion-height-d-n-pair-save-csv',
                                                                         n_clicks=0)],
                                            )
                                        ], style={'margin-left': '13%'}),
                                    ], style={'width': '33%'})
                                ], style={'display': 'flex', 'justify-content': 'center'})
                            ], style={'margin-top': '40px'}),
                        ], id='exclusion-height-result-collapse', is_open=False)
                    ], id='exclusion-height-tab-content')
                ], label='Exclusion height determination', tab_id='exclusion-height-tab', style={'margin-top': '10px'}),

                # Result summary tab
                dbc.Tab([
                    dash.html.Div([
                        dash.html.Div([
                            dash.html.Div([
                                dash.html.Div([
                                    dash.html.H4(['Fresnel modelling results']),
                                        dbc.Button('Export into .csv file', id='export-single-file-button',
                                                   color='primary',
                                                   n_clicks=0,
                                                   disabled=False),
                                    dbc.Table.from_dataframe(pd.DataFrame(
                                        {'Analysis': [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Variable': ['{layer}|{parameter}-{channel}'.format(
                                             layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                                             parameter=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                                             channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel) for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Value': [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for fresnel_inst in current_session.fresnel_analysis_instances]}),
                                        bordered=True, id='result-summary-fresnel-table'),
                                    dash.html.H4(['Exclusion height results']),
                                    dbc.Table.from_dataframe(pd.DataFrame(
                                        {'Analysis': [current_session.exclusion_height_analysis_instances[exclusion_inst].name
                                                           for exclusion_inst in
                                                           current_session.exclusion_height_analysis_instances],
                                         'Exclusion height mean [all]': ['{mean_} {all_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[0], 2), all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2))) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion RI mean [all]': ['{mean_} {all_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0], 4), all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :], decimals=4))) for exclusion_inst in current_session.exclusion_height_analysis_instances]}),
                                        bordered=True, id='result-summary-exclusion-table')
                                ], style={'margin-right': '100px', 'flex': '0 0 400px'}),
                                dash.html.Div([
                                    dash.html.Div([
                                        dash.dcc.Graph(id='summary-fresnel-barplot-graph',
                                                       figure=result_barplot_fig,
                                                       mathjax=True),
                                        dbc.DropdownMenu(
                                            id='barplot-save-dropdown',
                                            label='Save as...',
                                            color='info',
                                            children=[
                                                dbc.DropdownMenuItem('.PNG', id='barplot-save-png', n_clicks=0),
                                                dbc.DropdownMenuItem('.SVG', id='barplot-save-svg', n_clicks=0),
                                                dbc.DropdownMenuItem('.HTML', id='barplot-save-html', n_clicks=0),
                                                dbc.DropdownMenuItem('.csv', id='barplot-save-csv', n_clicks=0)],
                                            style={'margin-left': '299px'}
                                        )
                                    ]),
                                ]),
                            ], style={'margin-top': '20px', 'display': 'flex', 'justify-content': 'center'}),
                        ], style={'width': '95%' , 'margin-top': '1.9rem', 'margin-left': '5%'}),
                    ], id='summary-tab-content')
                ], label='Result summary and export', tab_id='summary-tab', style={'margin-top': '10px'}),
            ], id='analysis-tabs', active_tab='quantification-tab'),
        ], style={'margin-left': '2%', 'margin-right': '2%'})
    ])

    # Adding note to session log
    @dash.callback(
        dash.Output('console', 'value'),
        dash.Output('session-title', 'children'),
        dash.Output('test-input', 'value'),
        dash.Input('submit-button', 'n_clicks'),
        dash.Input('rename-session-button', 'n_clicks'),
        dash.State('test-input', 'value'),
        prevent_initial_call=True)
    def update_session_log(input1, input2, state1):

        global current_session
        if 'submit-button' == dash.ctx.triggered_id:

            new_message = current_session.log + '\n' + datetime.datetime.now().__str__()[0:16] + ' >> ' + state1
            current_session.log = new_message
            current_session.save_session()

            return new_message, dash.no_update, ''

        elif 'rename-session-button' == dash.ctx.triggered_id:

            current_session.update_name_and_location(state1)
            current_session.save_session()
            new_name = '{name_} - Session log'.format(name_=current_session.name)

            return dash.no_update, new_name, ''

    # Load in new measurement data and send a Store signal to other callbacks to update appropriately
    @dash.callback(
        dash.Output('loaded-new-measurement', 'data'),
        dash.Output('datapath-textfield', 'children'),
        dash.Output('batch-fresnel-analysis-files', 'data'),
        dash.Input('load-data', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-choose-files', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-button', 'n_clicks'),
        prevent_initial_call=True)
    def update_measurement_data(load_data, load_data_batch, button):

        global current_data_path
        global current_session
        global scanspeed
        global time_df
        global angles_df
        global ydata_df
        global reflectivity_df
        global sensorgram_df
        global sensorgram_df_selection
        global TIR_range_water_or_long_measurement
        global TIR_range_air_or_few_scans
        global TIR_range
        global default_data_folder

        if 'load-data' == dash.ctx.triggered_id:
            # Load measurement data and update session current data path
            current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(default_data_folder=default_data_folder)
            current_session.current_data_path = current_data_path
            current_session.save_session()

            # Calculate sensorgram (assume air or liquid medium for TIR calculation based on number of scans)
            if ydata_df.shape[0] > 50:
                TIR_range = TIR_range_water_or_long_measurement
            else:
                TIR_range = TIR_range_air_or_few_scans

            sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, TIR_range, scanspeed)

            # Offset to start at 0 degrees at 0 minutes
            sensorgram_df_selection = sensorgram_df
            sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                                   sensorgram_df_selection['SPR angle'][0]
            sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                                   sensorgram_df_selection['TIR angle'][0]

            return 'signal', ['Current measurement file:    ', current_data_path.split('/')[-1]], dash.no_update

        elif 'batch-fresnel-analysis-choose-files' == dash.ctx.triggered_id:
            print('Select the measurement data files (.csv)')
            batch_data_paths_ = select_files('Select the measurement data files', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')])
            return dash.no_update, dash.no_update, batch_data_paths_

        elif 'batch-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, None

    # Updating the sensor table with new values and properties
    @dash.callback(
        dash.Output('sensor-table', 'data'),  # Update sensor table data
        dash.Output('sensor-table-title', 'children'),  # Update sensor table title
        dash.Output('chosen-sensor-dropdown', 'children'),  # Update chosen sensor dropdown
        dash.Output('rename-sensor-modal', 'is_open'),
        dash.Output('remove-sensor-modal', 'is_open'),
        dash.Input({'type': 'sensor-list', 'index': dash.ALL}, 'n_clicks'),
        dash.Input('new-sensor-gold', 'n_clicks'),
        dash.Input('new-sensor-glass', 'n_clicks'),
        dash.Input('new-sensor-palladium', 'n_clicks'),
        dash.Input('new-sensor-platinum', 'n_clicks'),
        dash.Input('rename-sensor-button', 'n_clicks'),
        dash.Input('rename-sensor-confirm', 'n_clicks'),
        dash.Input('remove-sensor-button', 'n_clicks'),
        dash.Input('remove-sensor-confirm', 'n_clicks'),
        dash.Input('copy-sensor', 'n_clicks'),
        dash.Input('add-table-layer', 'n_clicks'),
        dash.Input('table-update-values', 'n_clicks'),
        dash.Input('table-select-fitted', 'n_clicks'),
        dash.Input('fresnel-reflectivity-run-finished', 'data'),
        dash.Input('batch-fresnel-analysis-finish', 'data'),
        dash.State('sensor-table', 'data'),
        dash.State('sensor-table', 'columns'),
        dash.State('sensor-table', 'active_cell'),
        dash.State('rename-sensor-input', 'value'),
        prevent_initial_call=True)
    def update_sensor_table(n_clicks_sensor_list, add_gold, add_sio2, add_palladium, add_platinum, rename_button,
                            rename_confirm, remove_button,
                            remove_confirm, click_copy, n_clicks_add_row, n_clicks_update, n_clicks_fitted,
                            fitted_result_update, batch_result_update, table_rows, table_columns, active_cell, sensor_name_):
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
        global current_session
        global current_data_path
        global reflectivity_df
        global scanspeed
        global TIR_range

        if 'new-sensor-gold' == dash.ctx.triggered_id:
            current_sensor = add_sensor_backend(current_session, current_data_path, sensor_metal='Au')
            current_sensor.name = 'Gold sensor'

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], TIR_range,
                                                scanspeed)
            current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(
                np.pi / 180 * TIR_angle)
            current_sensor.optical_parameters['n'] = current_sensor.refractive_indices

            # Save sensor and session
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'new-sensor-glass' == dash.ctx.triggered_id:
            current_sensor = add_sensor_backend(current_session, current_data_path, sensor_metal='SiO2')
            current_sensor.name = 'Glass sensor'

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], TIR_range,
                                                scanspeed)
            current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(
                np.pi / 180 * TIR_angle)
            current_sensor.optical_parameters['n'] = current_sensor.refractive_indices

            # Save sensor and session
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'new-sensor-palladium' == dash.ctx.triggered_id:
            current_sensor = add_sensor_backend(current_session, current_data_path, sensor_metal='Pd')
            current_sensor.name = 'Palladium sensor'

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], TIR_range,
                                                scanspeed)
            current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(
                np.pi / 180 * TIR_angle)
            current_sensor.optical_parameters['n'] = current_sensor.refractive_indices

            # Save sensor and session
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'new-sensor-platinum' == dash.ctx.triggered_id:
            current_sensor = add_sensor_backend(current_session, current_data_path, sensor_metal='Pt')
            current_sensor.name = 'Platinum sensor'

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], TIR_range, scanspeed)
            current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)
            current_sensor.optical_parameters['n'] = current_sensor.refractive_indices

            # Save sensor and session
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'copy-sensor' == dash.ctx.triggered_id:
            new_sensor = copy_sensor_backend(current_session, current_sensor)
            new_sensor.name = current_sensor.name + ' copy'
            current_sensor = new_sensor
            current_sensor.channel = current_data_path[-12:-4].replace('_', ' ')
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'rename-sensor-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, True, dash.no_update

        elif 'rename-sensor-confirm' == dash.ctx.triggered_id:

            # First remove previous sensor pickle object file
            old_path = current_session.location + r'\Sensors\S{id} {name}.pickle'.format(id=current_sensor.object_id,
                                                                                         name=current_sensor.name)
            os.remove(old_path)

            # Change sensor name and save new sensor pickle file and session
            current_sensor.name = sensor_name_
            current_session.save_sensor(current_sensor.object_id)
            current_session.save_session()

            data_rows = current_sensor.optical_parameters.to_dict('records')
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, False, dash.no_update

        elif 'remove-sensor-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, True

        elif 'remove-sensor-confirm' == dash.ctx.triggered_id:

            # Only allow removing sensors if there are at least 1 sensor in the list, otherwise do nothing
            if len(current_session.sensor_instances) > 1:

                removed = current_sensor
                try:
                    current_sensor = current_session.sensor_instances[1]
                except KeyError:  # In case the first few instances have already been removed, try the next one
                    failed = True
                    attempted = 1
                    max_attempts = len(current_session.sensor_instances) + 1
                    while failed and not attempted == max_attempts:
                        attempted += 1
                        try:
                            current_sensor = current_session.sensor_instances[attempted+1]
                        except KeyError:
                            continue
                        failed = False

                current_session.remove_sensor(removed.object_id)
                current_session.save_session()

                sensor_options = [
                    dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                         id={'type': 'sensor-list', 'index': sensor_id},
                                         n_clicks=0) for sensor_id in current_session.sensor_instances]

                data_rows = current_sensor.optical_parameters.to_dict('records')

                return data_rows, current_sensor.sensor_table_title, sensor_options, dash.no_update, False

            else:
                raise dash.exceptions.PreventUpdate

        elif 'remove-sensor-cancel' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, False

        elif 'table-update-values' == dash.ctx.triggered_id:

            # Update background sensor object
            current_sensor.optical_parameters = pd.DataFrame.from_records(table_rows)
            current_sensor.layer_thicknesses = current_sensor.optical_parameters['d [nm]'].to_numpy()
            current_sensor.refractive_indices = current_sensor.optical_parameters['n'].to_numpy()
            current_sensor.extinction_coefficients = current_sensor.optical_parameters['k'].to_numpy()
            current_sensor.fitted_var = current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index]

            # Save new sensor to session and Sensor folder
            current_session.save_session()
            current_session.save_sensor(current_sensor.object_id)

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-table-layer' == dash.ctx.triggered_id:
            table_rows.insert(-1, {c['id']: '' for c in table_columns})

            return table_rows, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'table-select-fitted' == dash.ctx.triggered_id:

            current_sensor.fitted_layer_index = (active_cell['row'], active_cell['column'])
            current_sensor.fitted_var = current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index]
            current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                sensor_number=current_sensor.object_id,
                sensor_name=current_sensor.name,
                channel=current_sensor.channel,
                fitted_layer=current_sensor.optical_parameters.iloc[active_cell['row'], 0],
                fitted_param=current_sensor.optical_parameters.columns[active_cell['column']])

            # Save new sensor to session and Sensor folder
            current_session.save_session()
            current_session.save_sensor(current_sensor.object_id)

            return table_rows, current_sensor.sensor_table_title, dash.no_update, dash.no_update, dash.no_update

        elif 'fresnel-reflectivity-run-finished' == dash.ctx.triggered_id:

            data_rows = current_sensor.optical_parameters.to_dict('records')

            return data_rows, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'batch-fresnel-analysis-finish' == dash.ctx.triggered_id:

            data_rows = current_sensor.optical_parameters.to_dict('records')

            sensor_options = [
                dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                     id={'type': 'sensor-list', 'index': sensor_id},
                                     n_clicks=0) for sensor_id in current_session.sensor_instances]

            return data_rows, current_sensor.sensor_table_title, sensor_options, dash.no_update, dash.no_update

        else:
            current_sensor = current_session.sensor_instances[dash.callback_context.triggered_id.index]

            data_rows = current_sensor.optical_parameters.to_dict('records')

            return data_rows, current_sensor.sensor_table_title, dash.no_update, dash.no_update, dash.no_update


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

    # Update the reflectivity plot in the Response quantification tab
    @dash.callback(
        dash.Output('quantification-reflectivity-graph', 'figure', allow_duplicate=True),
        dash.Input('quantification-reflectivity-add-data-trace', 'n_clicks'),
        dash.Input('quantification-reflectivity-add-fresnel-trace', 'n_clicks'),
        dash.Input('quantification-reflectivity-clear-traces', 'n_clicks'),
        dash.Input('quantification-reflectivity-save-png', 'n_clicks'),
        dash.Input('quantification-reflectivity-save-svg', 'n_clicks'),
        dash.Input('quantification-reflectivity-save-html', 'n_clicks'),
        dash.Input('quantification-reflectivity-save-csv', 'n_clicks'),
        dash.Input('quantification-sensorgram-graph', 'hoverData'),
        dash.Input('loaded-new-measurement', 'data'),
        dash.State('quantification-reflectivity-graph', 'figure'),
        dash.State('hover-selection-switch', 'value'),
        prevent_initial_call=True)
    def update_reflectivity_quantification_graph(add_data_trace, add_fresnel_trace, clear_traces, save_png, save_svg, save_html, save_csv, hoverData, loaded_new, figure_JSON, lock_hover):

        global ydata_df
        global angles_df
        global current_sensor
        global reflectivity_df

        figure_object = go.Figure(figure_JSON)

        # Update based on hover over sensorgram figure
        if 'quantification-sensorgram-graph' == dash.ctx.triggered_id:

            # First make sure no other traces has been added and the very first value is ignored
            if figure_object.data.__len__() == 1:

                # Then also make sure lock hover switch is set to inactive
                if lock_hover is False:
                    time_index = hoverData['points'][0]['pointIndex']

                    reflectivity_df['ydata'] = ydata_df.loc[time_index+1]

                    new_figure = go.Figure(go.Scatter(x=reflectivity_df['angles'],
                                                      y=reflectivity_df['ydata'],
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
            else:
                return dash.no_update

        # This adds a trace to the reflectivity plot from a separate measurement file. The trace data is not stored.
        elif 'quantification-reflectivity-add-data-trace' == dash.ctx.triggered_id:
            _, _, _, _, _, trace_reflectivity_df = load_csv_data(default_data_folder=default_data_folder)
            figure_object.add_trace(go.Scatter(x=trace_reflectivity_df['angles'],
                                               y=trace_reflectivity_df['ydata'],
                                               mode='lines',
                                               showlegend=False))

        # This adds a fresnel calculation trace to the reflectivity plot
        elif 'quantification-reflectivity-add-fresnel-trace' == dash.ctx.triggered_id:
            fresnel_coefficients = fresnel_calculation(None,
                                                       angles=reflectivity_df['angles'],
                                                       fitted_layer_index=current_sensor.fitted_layer_index,
                                                       wavelength=current_sensor.wavelength,
                                                       layer_thicknesses=current_sensor.layer_thicknesses,
                                                       n_re=current_sensor.refractive_indices,
                                                       n_im=current_sensor.extinction_coefficients,
                                                       ydata=None,
                                                       ydata_type=current_sensor.data_type,
                                                       polarization=1.0)
            figure_object.add_trace(go.Scatter(x=reflectivity_df['angles'],
                                               y=fresnel_coefficients,
                                               mode='lines',
                                               showlegend=False))

        # Clear added traces from the reflectivity plot
        elif 'quantification-reflectivity-clear-traces' == dash.ctx.triggered_id:
            new_figure = go.Figure(go.Scatter(x=reflectivity_df['angles'],
                                              y=reflectivity_df['ydata'],
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

        elif 'loaded-new-measurement' == dash.ctx.triggered_id:
            new_figure = go.Figure(go.Scatter(x=reflectivity_df['angles'],
                                              y=reflectivity_df['ydata'],
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

        elif 'quantification-reflectivity-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            fig_keys_x = ['x' + str(i) for i in range(len(figure_object.data))]
            fig_keys_y = ['y' + str(i) for i in range(len(figure_object.data))]
            fig_keys = [key for sublist in zip(fig_keys_x, fig_keys_y) for key in sublist]
            fig_values_x = []
            for i in range(len(figure_object.data)):
                fig_values_x.append(list(figure_object.data[i].x["_inputArray"].values())[:-3])
            fig_values_y = []
            for i in range(len(figure_object.data)):
                fig_values_y.append(list(figure_object.data[i].y["_inputArray"].values())[:-3])
            fig_values = [value for sublist in zip(fig_values_x, fig_values_y) for value in sublist]
            fig_df = pd.DataFrame(data={key:value for (key, value) in zip(fig_keys, fig_values)})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

    # Update the sensorgram plot in the Response quantification tab
    #TODO: Add our PureKinetics plot somewhere here as well
    @dash.callback(
        dash.Output('quantification-sensorgram-graph', 'figure'),
        dash.Input('quantification-sensorgram-save-png', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-svg', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-html', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-csv', 'n_clicks'),
        dash.Input('quantification-sensorgram-graph', 'clickData'),
        dash.Input('loaded-new-measurement', 'data'),
        dash.State('quantification-sensorgram-graph', 'figure'),
        prevent_initial_call=True)  # Adding this fixed a weird bug with graph not updating after firing clickData callbacks
    def update_sensorgram_quantification_tab(save_png, save_svg, save_html, save_csv, clickData, data_update, figure_JSON):

        figure_object = go.Figure(figure_JSON)
        global sensorgram_df_selection

        if 'quantification-sensorgram-graph' == dash.ctx.triggered_id:

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

        elif 'loaded-new-measurement' == dash.ctx.triggered_id:

            new_sensorgram_fig = go.Figure(go.Scatter(x=sensorgram_df_selection['time'],
                                                      y=sensorgram_df_selection['SPR angle'],
                                                      name='SPR angle',
                                                      line_color='#636efa'))

            new_sensorgram_fig.add_trace(go.Scatter(x=sensorgram_df_selection['time'],
                                                    y=sensorgram_df_selection['TIR angle'],
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

        elif 'quantification-sensorgram-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            fig_df = pd.DataFrame(data={'Time': list(figure_object.data[0].x["_inputArray"].values())[:-3], 'SPR': list(figure_object.data[0].y["_inputArray"].values())[:-3], 'TIR': list(figure_object.data[1].y["_inputArray"].values())[:-3]})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate


    # Update the reflectivity plot in the Fresnel fitting tab
    @dash.callback(
        dash.Output('fresnel-reflectivity-graph', 'figure'),
        dash.Output('result-summary-fresnel-table', 'children'),
        dash.Output('summary-fresnel-barplot-graph', 'figure'),
        dash.Output('fresnel-reflectivity-run-finished', 'data'),
        dash.Output('fresnel-analysis-dropdown', 'children'),
        dash.Output('fresnel-analysis-option-collapse', 'is_open'),
        dash.Output('add-fresnel-analysis-modal', 'is_open'),
        dash.Output('remove-fresnel-analysis-modal', 'is_open'),
        dash.Output('fresnel-fit-option-rangeslider', 'value'),
        dash.Output('fresnel-fit-option-iniguess', 'value'),
        dash.Output('fresnel-fit-option-lowerbound', 'value'),
        dash.Output('fresnel-fit-option-upperbound', 'value'),
        dash.Output('fresnel-fit-option-extinctionslider', 'value', allow_duplicate=True),
        dash.Output('fresnel-fit-result', 'children'),
        dash.Output('fresnel-fit-sensor', 'children'),
        dash.Output('exclusion-choose-background-dropdown', 'options'),
        dash.Output('fresnel-fit-option-rangeslider', 'marks'),
        dash.Output('fresnel-fit-option-rangeslider', 'min'),
        dash.Output('fresnel-fit-option-rangeslider', 'max'),
        dash.Output('fresnel-fit-datapath', 'children'),
        dash.Output('rename-fresnel-analysis-modal', 'is_open'),
        dash.Output('batch-fresnel-analysis-done', 'data'),
        dash.Output('fresnel-analysis-offset-fit', 'value'),
        dash.Output('fresnel-analysis-elastomer-fit', 'value', allow_duplicate=True),
        dash.Input('fresnel-reflectivity-run-model', 'n_clicks'),
        dash.Input('add-fresnel-analysis-button', 'n_clicks'),
        dash.Input('add-fresnel-analysis-confirm', 'n_clicks'),
        dash.Input('remove-fresnel-analysis-button', 'n_clicks'),
        dash.Input('remove-fresnel-analysis-confirm', 'n_clicks'),
        dash.Input('remove-fresnel-analysis-cancel', 'n_clicks'),
        dash.Input('fresnel-fit-option-rangeslider', 'value'),
        dash.Input({'type': 'fresnel-analysis-list', 'index': dash.ALL}, 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-png', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-svg', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-html', 'n_clicks'),
        dash.Input('fresnel-reflectivity-save-csv', 'n_clicks'),
        dash.Input('rename-fresnel-analysis-button', 'n_clicks'),
        dash.Input('rename-fresnel-analysis-confirm', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-start', 'data'),
        dash.State('fresnel-analysis-name-input', 'value'),
        dash.State('fresnel-reflectivity-graph', 'figure'),
        dash.State('fresnel-fit-option-rangeslider', 'value'),
        dash.State('fresnel-fit-option-iniguess', 'value'),
        dash.State('fresnel-fit-option-lowerbound', 'value'),
        dash.State('fresnel-fit-option-upperbound', 'value'),
        dash.State('fresnel-fit-option-extinctionslider', 'value'),
        dash.State('rename-fresnel-analysis-input', 'value'),
        dash.State('batch-fresnel-analysis-files', 'data'),
        dash.State('batch-fresnel-analysis-background-sensors', 'data'),
        dash.State('batch-fresnel-analysis-radio-selection', 'value'),
        dash.State('batch-fresnel-analysis-newlayer-radio-selection', 'value'),
        dash.State('batch-fresnel-analysis-example-sensor-dropdown', 'value'),
        dash.State('batch-fresnel-analysis-example-analysis-dropdown', 'value'),
        dash.State('fresnel-analysis-offset-fit', 'value'),
        dash.State('fresnel-analysis-elastomer-fit', 'value'),
        dash.State('fresnel-fit-option-pfactor', 'value'),
        prevent_initial_call=True)
    def update_reflectivity_fresnel_graph(run_model, add_button, add_confirm_button, remove_button, remove_confirm, remove_cancel, rangeslider_inp,
                                          selected_fresnel_object, save_png, save_svg, save_html, save_csv, rename_button, rename_confirm, batch_start_signal, analysis_name, figure_JSON, rangeslider_state, ini_guess,
                                          lower_bound, upper_bound,
                                          extinction_correction, analysis_name_, batch_files, background_sensors, batch_radio_selection, batch_newlayer_radio_selection, batch_sensor_index, batch_analysis_index, offset_fit_flag, elastomer_fit_flag, polarization_factor):

        global current_fresnel_analysis
        global current_data_path
        global current_session
        global reflectivity_df
        global current_sensor
        global TIR_range
        global scanspeed

        figure_object = go.Figure(figure_JSON)

        if 'fresnel-fit-option-rangeslider' == dash.ctx.triggered_id:

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

            return new_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'fresnel-reflectivity-run-model' == dash.ctx.triggered_id:

            # Set analysis options from dash app
            current_fresnel_analysis.angle_range = rangeslider_state
            current_fresnel_analysis.polarization = polarization_factor
            current_fresnel_analysis.sensor_object_label = 'Sensor: ' + current_sensor.sensor_table_title
            current_fresnel_analysis.fit_offset = offset_fit_flag
            current_fresnel_analysis.fit_prism_k = elastomer_fit_flag

            # Determine number of simultaneous fitting variables
            if current_fresnel_analysis.fit_offset and current_fresnel_analysis.fit_prism_k:
                current_fresnel_analysis.ini_guess = np.array([ini_guess, current_fresnel_analysis.y_offset, 0.001])
                current_fresnel_analysis.bounds = [(lower_bound, -np.inf, 0), (upper_bound, np.inf, 0.1)]
                current_fresnel_analysis.extinction_correction = 0

            elif current_fresnel_analysis.fit_offset and not current_fresnel_analysis.fit_prism_k:
                current_fresnel_analysis.ini_guess = np.array([ini_guess, current_fresnel_analysis.y_offset])
                current_fresnel_analysis.bounds = [(lower_bound, -np.inf), (upper_bound, np.inf)]
                current_fresnel_analysis.extinction_correction = extinction_correction

            elif not current_fresnel_analysis.fit_offset:
                current_fresnel_analysis.ini_guess = np.array([ini_guess])
                current_fresnel_analysis.bounds = [lower_bound, upper_bound]
                current_fresnel_analysis.extinction_correction = 0

            # Run calculations and modelling
            fresnel_df = current_fresnel_analysis.model_reflectivity_trace()

            # Update current sensor object with the fit result and prism extinction value
            current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index] = round(current_fresnel_analysis.fitted_result[0], 4)

            if not current_fresnel_analysis.fit_prism_k:
                current_sensor.optical_parameters.iloc[(0, 3)] = current_sensor.extinction_coefficients[0]
            else:
                current_sensor.optical_parameters.iloc[(0, 3)] = round(current_fresnel_analysis.fitted_result[2], 5)

            # Save session and analysis object
            current_session.save_session()
            current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)

            # Fit result text
            result = 'Fit result: {res}'.format(res=round(current_fresnel_analysis.fitted_result[0], 4))

            # Plot fitted trace
            new_figure = go.Figure(go.Scatter(x=current_fresnel_analysis.measurement_data['angles'],
                                              y=current_fresnel_analysis.measurement_data['ydata'],
                                              mode='lines',
                                              showlegend=False,
                                              line_color='#636efa'
                                              ))
            new_figure.add_trace(go.Scatter(x=fresnel_df['angles'],
                                            y=fresnel_df['ydata'],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='#ef553b'
                                            ))
            new_figure.add_trace(go.Scatter(x=[current_fresnel_analysis.angle_range[0], current_fresnel_analysis.angle_range[0]],
                                            y=[min(current_fresnel_analysis.measurement_data['ydata']), max(current_fresnel_analysis.measurement_data['ydata'])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            new_figure.add_trace(go.Scatter(x=[current_fresnel_analysis.angle_range[1], current_fresnel_analysis.angle_range[1]],
                                            y=[min(current_fresnel_analysis.measurement_data['ydata']), max(current_fresnel_analysis.measurement_data['ydata'])],
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

            table_header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td(current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td('{layer}|{parameter}-{channel}'.format(
                                             layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                                             parameter=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                                             channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)), dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            x_barplot = [[current_session.fresnel_analysis_instances[
                          fresnel_inst].fitted_layer for
                      fresnel_inst in current_session.fresnel_analysis_instances],[current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in current_session.fresnel_analysis_instances]]
            y_barplot = [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for fresnel_inst in current_session.fresnel_analysis_instances]
            result_barplot_fig = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
            result_barplot_fig.update_layout(
                yaxis_title='Fitted value',
                font_family='Balto',
                font_size=19,
                margin_r=25,
                margin_l=60,
                margin_t=40,
                template='simple_white',
                uirevision=True,
                height=600)
            result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
            result_barplot_fig.update_yaxes(mirror=True, showline=True)

            return new_figure, fresnel_result_summary_dataframe, result_barplot_fig, 'finished', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, result, current_fresnel_analysis.sensor_object_label, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            current_fresnel_analysis = add_fresnel_model_object(current_session, current_sensor, current_data_path, reflectivity_df, TIR_range, scanspeed, analysis_name)

            # Calculate initial intensity offset from data
            FR_y = fresnel_calculation(
                angles=reflectivity_df['angles'].iloc[
                       reflectivity_df['ydata'].idxmin()-1:reflectivity_df['ydata'].idxmin()+1],
                wavelength=current_fresnel_analysis.sensor_object.wavelength,
                layer_thicknesses=current_fresnel_analysis.sensor_object.layer_thicknesses,
                n_re=current_fresnel_analysis.sensor_object.refractive_indices,
                n_im=current_fresnel_analysis.sensor_object.extinction_coefficients,
                ydata_type='R',
                polarization=current_fresnel_analysis.polarization)
            current_fresnel_analysis.y_offset = reflectivity_df['ydata'].min() - FR_y[1]  # Can't calculate only 1 angle, so use middle of 3 around minimum

            current_fresnel_analysis.ini_guess = np.array([float(current_sensor.fitted_var), current_fresnel_analysis.y_offset])
            current_fresnel_analysis.bounds = [(current_fresnel_analysis.ini_guess[0] / 4, -np.inf), (current_fresnel_analysis.ini_guess[0] * 2, np.inf)]
            current_fresnel_analysis.angle_range = [reflectivity_df['angles'].iloc[reflectivity_df['ydata'].idxmin()-auto_angle_range_points[0]], reflectivity_df['angles'].iloc[reflectivity_df['ydata'].idxmin()+auto_angle_range_points[1]]]

            current_session.save_session()
            current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)

            analysis_options = [
                dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                     id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                                     n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]

            exclusion_analysis_dropdown = [{'label': 'FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name, 'value': fresnel_id} for fresnel_id in current_session.fresnel_analysis_instances]

            current_fresnel_analysis.sensor_object_label = 'Sensor: ' + current_sensor.sensor_table_title

            # Update fresnel plot with current measurement data
            new_figure = go.Figure(go.Scatter(x=current_fresnel_analysis.measurement_data['angles'],
                                              y=current_fresnel_analysis.measurement_data['ydata'],
                                              mode='lines',
                                              showlegend=False,
                                              line_color='#636efa'
                                              ))

            # Update angle range markers
            angle_range_marks = {mark_ind: str(mark_ind) for mark_ind in range(current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 1)}

            # Add lines for angle range
            new_figure.add_trace(
                go.Scatter(x=[current_fresnel_analysis.angle_range[0], current_fresnel_analysis.angle_range[0]],
                           y=[min(current_fresnel_analysis.measurement_data['ydata']),
                              max(current_fresnel_analysis.measurement_data['ydata'])],
                           mode='lines',
                           showlegend=False,
                           line_color='black',
                           line_dash='dash'
                           ))

            new_figure.add_trace(
                go.Scatter(x=[current_fresnel_analysis.angle_range[1], current_fresnel_analysis.angle_range[1]],
                           y=[min(current_fresnel_analysis.measurement_data['ydata']),
                              max(current_fresnel_analysis.measurement_data['ydata'])],
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

            # Check bounds structure
            if type(current_fresnel_analysis.bounds[0]) is not tuple:
                lower_bound_ = current_fresnel_analysis.bounds[0]
                upper_bound_ = current_fresnel_analysis.bounds[1]
            else:
                lower_bound_ = current_fresnel_analysis.bounds[0][0]
                upper_bound_ = current_fresnel_analysis.bounds[1][0]

            return new_figure, dash.no_update, dash.no_update, dash.no_update, analysis_options, dash.no_update, False, dash.no_update, current_fresnel_analysis.angle_range, current_fresnel_analysis.ini_guess[0], \
            lower_bound_, upper_bound_, current_fresnel_analysis.extinction_correction, 'Fit result: None', current_fresnel_analysis.sensor_object_label, exclusion_analysis_dropdown, angle_range_marks, current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 'Data path: \n' + current_fresnel_analysis.initial_data_path, dash.no_update, dash.no_update, current_fresnel_analysis.fit_offset, current_fresnel_analysis.fit_prism_k

        elif 'rename-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update

        elif 'rename-fresnel-analysis-confirm' == dash.ctx.triggered_id:

            # First remove previous fresnel analysis pickle object file
            old_path = current_session.location + r'\Analysis instances\FM{id} {name}.pickle'.format(id=current_fresnel_analysis.object_id, name=current_fresnel_analysis.name)
            os.remove(old_path)

            # Change fresnel analysis name and save new fresnel analysis pickle file and session
            current_fresnel_analysis.name = analysis_name_
            current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
            current_session.save_session()

            analysis_options = [dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                                     id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                                                     n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]

            table_header = [dash.html.Thead(
                dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td(current_session.fresnel_analysis_instances[fresnel_inst].name),
                 dash.html.Td('{layer}|{parameter}-{channel}'.format(
                     layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                     parameter=
                     current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[
                         current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                     channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                 dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for
                                           fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            x_barplot = [[current_session.fresnel_analysis_instances[
                          fresnel_inst].fitted_layer for
                      fresnel_inst in current_session.fresnel_analysis_instances],
                         [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in
                          current_session.fresnel_analysis_instances]]
            y_barplot = [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for
                         fresnel_inst in current_session.fresnel_analysis_instances]
            result_barplot_fig = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
            result_barplot_fig.update_layout(
                yaxis_title='Fitted value',
                font_family='Balto',
                font_size=19,
                margin_r=25,
                margin_l=60,
                margin_t=40,
                template='simple_white',
                uirevision=True,
                height=600)
            result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
            result_barplot_fig.update_yaxes(mirror=True, showline=True)

            return dash.no_update, fresnel_result_summary_dataframe, result_barplot_fig, dash.no_update, analysis_options, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            if len(current_session.fresnel_analysis_instances) > 1:

                # Pop out the current fresnel analysis object from the session, delete its .pickle file and make the first instance the current one
                removed = current_fresnel_analysis

                # TODO: Check that this works!
                try:
                    current_fresnel_analysis = current_session.fresnel_analysis_instances[1]
                except KeyError:  # In case the first few instances have already been removed, try the next one
                    failed = True
                    attempted = 1
                    max_attempts = len(current_session.fresnel_analysis_instances) + 1
                    while failed and not attempted == max_attempts:
                        attempted += 1
                        try:
                            current_fresnel_analysis = current_session.fresnel_analysis_instances[attempted+1]
                        except KeyError:
                            continue
                        failed = False

                current_session.remove_fresnel_analysis(removed.object_id)
                current_session.save_session()

                # Update all analysis options accordingly
                analysis_options = [
                    dbc.DropdownMenuItem(
                        'FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                        id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                        n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]

                exclusion_analysis_dropdown = [{'label': 'FM' + str(fresnel_id) + ' ' +
                                                         current_session.fresnel_analysis_instances[fresnel_id].name,
                                                'value': fresnel_id} for fresnel_id in
                                               current_session.fresnel_analysis_instances]

                if current_fresnel_analysis.fitted_result is not None:
                    result = 'Fit result: {res}'.format(res=round(current_fresnel_analysis.fitted_result[0], 4))
                else:
                    result = 'Fit result: None'

                # If the current loaded measurement data is not the same as the analysis object, use a different color
                if current_data_path != current_fresnel_analysis.initial_data_path:
                    line_color_value = '#00CC96'
                else:
                    line_color_value = '#636EFA'

                # Plot figures
                new_figure = go.Figure(go.Scatter(x=current_fresnel_analysis.measurement_data['angles'],
                                                  y=current_fresnel_analysis.measurement_data['ydata'],
                                                  mode='lines',
                                                  showlegend=False,
                                                  line_color=line_color_value
                                                  ))
                if current_fresnel_analysis.fitted_data is not None:
                    new_figure.add_trace(go.Scatter(x=current_fresnel_analysis.fitted_data['angles'],
                                                    y=current_fresnel_analysis.fitted_data['ydata'],
                                                    mode='lines',
                                                    showlegend=False,
                                                    line_color='#ef553b'
                                                    ))
                new_figure.add_trace(
                    go.Scatter(x=[current_fresnel_analysis.angle_range[0], current_fresnel_analysis.angle_range[0]],
                               y=[min(current_fresnel_analysis.measurement_data['ydata']),
                                  max(current_fresnel_analysis.measurement_data['ydata'])],
                               mode='lines',
                               showlegend=False,
                               line_color='black',
                               line_dash='dash'
                               ))
                new_figure.add_trace(
                    go.Scatter(x=[current_fresnel_analysis.angle_range[1], current_fresnel_analysis.angle_range[1]],
                               y=[min(current_fresnel_analysis.measurement_data['ydata']),
                                  max(current_fresnel_analysis.measurement_data['ydata'])],
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

                # Check bounds structure
                if type(current_fresnel_analysis.bounds[0]) is not tuple:
                    lower_bound_ = current_fresnel_analysis.bounds[0]
                    upper_bound_ = current_fresnel_analysis.bounds[1]
                else:
                    lower_bound_ = current_fresnel_analysis.bounds[0][0]
                    upper_bound_ = current_fresnel_analysis.bounds[1][0]

                table_header = [dash.html.Thead(
                    dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
                table_body = [dash.html.Tbody([dash.html.Tr(
                    [dash.html.Td(current_session.fresnel_analysis_instances[fresnel_inst].name),
                     dash.html.Td('{layer}|{parameter}-{channel}'.format(
                         layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                         parameter=current_session.fresnel_analysis_instances[
                             fresnel_inst].sensor_object.optical_parameters.columns[
                             current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                         channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                     dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))])
                                               for fresnel_inst in current_session.fresnel_analysis_instances])]
                fresnel_result_summary_dataframe = table_header + table_body

                x_barplot = [[current_session.fresnel_analysis_instances[
                          fresnel_inst].fitted_layer for
                      fresnel_inst in current_session.fresnel_analysis_instances],
                             [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in
                              current_session.fresnel_analysis_instances]]
                y_barplot = [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for
                             fresnel_inst in current_session.fresnel_analysis_instances]
                result_barplot_fig = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
                result_barplot_fig.update_layout(
                    yaxis_title='Fitted value',
                    font_family='Balto',
                    font_size=19,
                    margin_r=25,
                    margin_l=60,
                    margin_t=40,
                    template='simple_white',
                    uirevision=True,
                    height=600)
                result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
                result_barplot_fig.update_yaxes(mirror=True, showline=True)

                return new_figure, fresnel_result_summary_dataframe, result_barplot_fig, dash.no_update, analysis_options, dash.no_update, dash.no_update, False, current_fresnel_analysis.angle_range, current_fresnel_analysis.ini_guess[0], \
                    lower_bound_, upper_bound_, current_fresnel_analysis.extinction_correction, result, current_fresnel_analysis.sensor_object_label, exclusion_analysis_dropdown, dash.no_update, dash.no_update, dash.no_update, 'Data path: \n' + current_fresnel_analysis.initial_data_path, False, dash.no_update, current_fresnel_analysis.fit_offset, current_fresnel_analysis.fit_prism_k

            # If deleting the last fresnel analysis object
            else:
                try:
                    current_session.remove_fresnel_analysis(current_fresnel_analysis.object_id)
                except AttributeError:
                    pass  # There was no object at all, this will cause big problems if AttributeError can happen for other reasons though
                current_fresnel_analysis = None
                current_session.save_session()

                return figure_object, pd.DataFrame({'Analysis': [''], 'Variable': [''], 'Value': ['']}), go.Figure(go.Bar(x=[0], y=[0])), dash.no_update, [], False, False, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, [], dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-fresnel-analysis-cancel' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'batch-fresnel-analysis-start' == dash.ctx.triggered_id:
            # Load the example sensor and fresnel model into memory for pulling base settings
            example_sensor_object = current_session.sensor_instances[batch_sensor_index]
            example_analysis_object = current_session.fresnel_analysis_instances[batch_analysis_index]

            # Conditional for batch analysis radio button selection
            if batch_radio_selection == 0:  # Copy selected example layer structure

                # Use the same layer structure copied from selected example sensor object
                for file_path in batch_files:

                    # Load data from measurement file using load_csv_data
                    _, _, _, _, _, next_reflectivity_df_ = load_csv_data(path=file_path)

                    # Add copy of sensor object to session and set parameters
                    next_sensor = copy_sensor_backend(current_session, example_sensor_object)
                    try:
                        next_sensor.name = file_path.split('/')[-1][15:-10].replace('_', ' ')
                    except:
                        next_sensor.name = example_sensor_object.name

                    current_sensor = next_sensor
                    current_sensor.channel = file_path[-12:-4].replace('_', ' ')
                    TIR_angle, _, _ = TIR_determination(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'], example_analysis_object.TIR_range,
                                                        example_analysis_object.scanspeed)
                    current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(
                        np.pi / 180 * TIR_angle)
                    current_sensor.optical_parameters['n'] = current_sensor.refractive_indices
                    current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                        sensor_number=current_sensor.object_id,
                        sensor_name=current_sensor.name,
                        channel=current_sensor.channel,
                        fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                        fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])
                    current_session.save_sensor(current_sensor.object_id)

                    # Add fresnel model object to session
                    current_fresnel_analysis = add_fresnel_model_object(current_session, current_sensor,
                                                                        file_path, next_reflectivity_df_, example_analysis_object.TIR_range,
                                                                        example_analysis_object.scanspeed, current_sensor.name)
                    # Calculate angle range based on measured data
                    current_fresnel_analysis.angle_range = [
                        next_reflectivity_df_['angles'].iloc[next_reflectivity_df_['ydata'].idxmin() - auto_angle_range_points[0]],
                        next_reflectivity_df_['angles'].iloc[next_reflectivity_df_['ydata'].idxmin() + auto_angle_range_points[1]]]

                    # Set analysis options from example analysis objects
                    current_fresnel_analysis.ini_guess = example_analysis_object.ini_guess
                    current_fresnel_analysis.bounds = example_analysis_object.bounds
                    current_fresnel_analysis.polarization = example_analysis_object.polarization
                    current_fresnel_analysis.extinction_correction = example_analysis_object.extinction_correction
                    current_fresnel_analysis.y_offset = example_analysis_object.y_offset
                    current_fresnel_analysis.fit_prism_k = example_analysis_object.fit_prism_k

                    # Run calculations and modelling
                    fresnel_df = current_fresnel_analysis.model_reflectivity_trace()

                    # Update current sensor object with the fit result and prism extinction value
                    current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index] = round(
                        current_fresnel_analysis.fitted_result[0], 4)

                    if not current_fresnel_analysis.fit_prism_k:
                        current_sensor.optical_parameters.iloc[(0, 3)] = current_sensor.extinction_coefficients[0]
                    else:
                        current_sensor.optical_parameters.iloc[(0, 3)] = round(
                            current_fresnel_analysis.fitted_result[2], 5)

                    current_fresnel_analysis.sensor_object_label = 'Sensor: ' + current_sensor.sensor_table_title

                    # Save session and analysis object
                    current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
                    current_session.save_session()

            elif batch_radio_selection == 1:  # Use individual backgrounds and add new layer

                # Use the same layer structure copied from selected example sensor object
                for file_path, sensor_id in zip(batch_files, background_sensors):

                    # Load data from measurement file using load_csv_data
                    _, batch_current_data_path, _, batch_ydata_df, _, next_reflectivity_df_ = load_csv_data(path=file_path)

                    # Select background sensor
                    background_sensor_object = current_session.sensor_instances[sensor_id]

                    if batch_newlayer_radio_selection == 0:
                        # Directly modify background sensor instance
                        current_sensor = background_sensor_object

                    elif batch_newlayer_radio_selection == 1:

                        # Add copy of sensor object to session
                        current_sensor = copy_sensor_backend(current_session, background_sensor_object)

                    try:
                        current_sensor.name = file_path.split('/')[-1][15:-10].replace('_', ' ')
                    except:
                        current_sensor.name = background_sensor_object.name + ' + ' + example_sensor_object.optical_parameters.iloc[-2, 0]

                    # Add example layer row and values, also convert other parameters
                    current_sensor.optical_parameters.loc[len(current_sensor.optical_parameters)-1.5] = example_sensor_object.optical_parameters.loc[len(example_sensor_object.optical_parameters) - 2]
                    current_sensor.optical_parameters = current_sensor.optical_parameters.sort_index().reset_index(drop=True)
                    current_sensor.layer_thicknesses = current_sensor.optical_parameters['d [nm]'].to_numpy()
                    current_sensor.refractive_indices = current_sensor.optical_parameters['n'].to_numpy()
                    current_sensor.extinction_coefficients = current_sensor.optical_parameters['k'].to_numpy()

                    # Calculate TIR angle and update bulk RI
                    TIR_angle, _, _ = TIR_determination(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'],
                                                        example_analysis_object.TIR_range,
                                                        example_analysis_object.scanspeed)
                    current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(
                        np.pi / 180 * TIR_angle)
                    current_sensor.optical_parameters['n'] = current_sensor.refractive_indices

                    # Select correct variable to fit
                    current_sensor.fitted_layer_index = example_sensor_object.fitted_layer_index
                    current_sensor.fitted_var = current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index]

                    # Update sensor title
                    current_sensor.channel = file_path[-12:-4].replace('_', ' ')
                    current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                        sensor_number=current_sensor.object_id,
                        sensor_name=current_sensor.name,
                        channel=current_sensor.channel,
                        fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                        fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

                    current_session.save_sensor(current_sensor.object_id)

                    # Add fresnel model object to session
                    current_fresnel_analysis = add_fresnel_model_object(current_session, current_sensor,
                                                                        file_path, next_reflectivity_df_,
                                                                        example_analysis_object.TIR_range,
                                                                        example_analysis_object.scanspeed,
                                                                        current_sensor.name + ' (S' + str(
                                                                            current_sensor.object_id) + ') ')

                    # Calculate angle range based on measured data
                    current_fresnel_analysis.angle_range = [
                        next_reflectivity_df_['angles'].iloc[
                            next_reflectivity_df_['ydata'].idxmin() - auto_angle_range_points[0]],
                        next_reflectivity_df_['angles'].iloc[
                            next_reflectivity_df_['ydata'].idxmin() + auto_angle_range_points[1]]]

                    # Set analysis options from example analysis objects
                    current_fresnel_analysis.ini_guess = example_analysis_object.ini_guess
                    current_fresnel_analysis.bounds = example_analysis_object.bounds
                    current_fresnel_analysis.polarization = example_analysis_object.polarization
                    current_fresnel_analysis.extinction_correction = example_analysis_object.extinction_correction
                    current_fresnel_analysis.y_offset = example_analysis_object.y_offset
                    current_fresnel_analysis.fit_prism_k = example_analysis_object.fit_prism_k

                    # Run calculations and modelling
                    fresnel_df = current_fresnel_analysis.model_reflectivity_trace()

                    # Update current sensor object with the fit result and prism extinction value
                    current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index] = round(
                        current_fresnel_analysis.fitted_result[0], 4)

                    if not current_fresnel_analysis.fit_prism_k:
                        current_sensor.optical_parameters.iloc[(0, 3)] = current_sensor.extinction_coefficients[0]
                    else:
                        current_sensor.optical_parameters.iloc[(0, 3)] = round(
                            current_fresnel_analysis.fitted_result[2], 5)

                    current_fresnel_analysis.sensor_object_label = 'Sensor: ' + current_sensor.sensor_table_title

                    # Save session and analysis object
                    current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
                    current_session.save_session()

            # Fit result text
            result = 'Fit result: {res}'.format(res=round(current_fresnel_analysis.fitted_result[0], 4))

            # Plot fitted trace
            new_figure = go.Figure(go.Scatter(x=current_fresnel_analysis.measurement_data['angles'],
                                              y=current_fresnel_analysis.measurement_data['ydata'],
                                              mode='lines',
                                              showlegend=False,
                                              line_color='#636efa'
                                              ))
            new_figure.add_trace(go.Scatter(x=fresnel_df['angles'],
                                            y=fresnel_df['ydata'],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='#ef553b'
                                            ))
            new_figure.add_trace(
                go.Scatter(x=[current_fresnel_analysis.angle_range[0], current_fresnel_analysis.angle_range[0]],
                           y=[min(current_fresnel_analysis.measurement_data['ydata']),
                              max(current_fresnel_analysis.measurement_data['ydata'])],
                           mode='lines',
                           showlegend=False,
                           line_color='black',
                           line_dash='dash'
                           ))
            new_figure.add_trace(
                go.Scatter(x=[current_fresnel_analysis.angle_range[1], current_fresnel_analysis.angle_range[1]],
                           y=[min(current_fresnel_analysis.measurement_data['ydata']),
                              max(current_fresnel_analysis.measurement_data['ydata'])],
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

            # Check bounds structure
            if type(current_fresnel_analysis.bounds[0]) is not tuple:
                lower_bound_ = current_fresnel_analysis.bounds[0]
                upper_bound_ = current_fresnel_analysis.bounds[1]
            else:
                lower_bound_ = current_fresnel_analysis.bounds[0][0]
                upper_bound_ = current_fresnel_analysis.bounds[1][0]

            analysis_options = [
                dbc.DropdownMenuItem(
                    'FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                    id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                    n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]

            table_header = [dash.html.Thead(
                dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td(current_session.fresnel_analysis_instances[fresnel_inst].name),
                 dash.html.Td('{layer}|{parameter}-{channel}'.format(
                     layer=
                     current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                     parameter=
                     current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[
                         current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                     channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                 dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for
                                           fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            x_barplot = [[current_session.fresnel_analysis_instances[
                          fresnel_inst].fitted_layer for
                      fresnel_inst in current_session.fresnel_analysis_instances],
                         [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in
                          current_session.fresnel_analysis_instances]]
            y_barplot = [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for
                         fresnel_inst in current_session.fresnel_analysis_instances]
            result_barplot_fig = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
            result_barplot_fig.update_layout(
                yaxis_title='Fitted value',
                font_family='Balto',
                font_size=19,
                margin_r=25,
                margin_l=60,
                margin_t=40,
                template='simple_white',
                uirevision=True,
                height=600)
            result_barplot_fig.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
            result_barplot_fig.update_yaxes(mirror=True, showline=True)

            return new_figure, fresnel_result_summary_dataframe, result_barplot_fig, dash.no_update, analysis_options, dash.no_update, dash.no_update, dash.no_update, current_fresnel_analysis.angle_range, current_fresnel_analysis.ini_guess[0], lower_bound_, upper_bound_, dash.no_update, result, current_fresnel_analysis.sensor_object_label, dash.no_update, dash.no_update, current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 'Data path: \n' + current_fresnel_analysis.initial_data_path, dash.no_update, 'finished', dash.no_update, dash.no_update

        elif 'fresnel-reflectivity-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            fig_keys_x = ['x' + str(i) for i in range(len(figure_object.data))]
            fig_keys_y = ['y' + str(i) for i in range(len(figure_object.data))]
            fig_keys = [key for sublist in zip(fig_keys_x, fig_keys_y) for key in sublist]
            fig_values_x = []
            for i in range(len(figure_object.data)):
                fig_values_x.append(list(figure_object.data[i].x["_inputArray"].values())[:-3])
            fig_values_y = []
            for i in range(len(figure_object.data)):
                fig_values_y.append(list(figure_object.data[i].y["_inputArray"].values())[:-3])
            fig_values = [value for sublist in zip(fig_values_x, fig_values_y) for value in sublist]
            fig_df = pd.DataFrame(data={key:value for (key, value) in zip(fig_keys, fig_values)})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        # Updating the fresnel fit graph when a different model object is selected in the fresnel analysis list
        else:
            current_fresnel_analysis = current_session.fresnel_analysis_instances[
                dash.callback_context.triggered_id.index]

            if current_fresnel_analysis.fitted_result is not None:
                result = 'Fit result: {res}'.format(res=round(current_fresnel_analysis.fitted_result[0], 4))
            else:
                result = 'Fit result: None'

            # If the current loaded measurement data is not the same as the analysis object, use a different color
            if current_data_path != current_fresnel_analysis.initial_data_path:
                line_color_value = '#00CC96'
            else:
                line_color_value = '#636EFA'

            # Calculate angle range marks
            angle_range_marks = {mark_ind: str(mark_ind) for mark_ind in range(current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 1)}

            # Plot figures
            new_figure = go.Figure(go.Scatter(x=current_fresnel_analysis.measurement_data['angles'],
                                              y=current_fresnel_analysis.measurement_data['ydata'],
                                              mode='lines',
                                              showlegend=False,
                                              line_color=line_color_value
                                              ))
            if current_fresnel_analysis.fitted_data is not None:
                new_figure.add_trace(go.Scatter(x=current_fresnel_analysis.fitted_data['angles'],
                                                y=current_fresnel_analysis.fitted_data['ydata'],
                                                mode='lines',
                                                showlegend=False,
                                                line_color='#ef553b'
                                                ))
            new_figure.add_trace(go.Scatter(x=[current_fresnel_analysis.angle_range[0], current_fresnel_analysis.angle_range[0]],
                                            y=[min(current_fresnel_analysis.measurement_data['ydata']), max(current_fresnel_analysis.measurement_data['ydata'])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            new_figure.add_trace(go.Scatter(x=[current_fresnel_analysis.angle_range[1], current_fresnel_analysis.angle_range[1]],
                                            y=[min(current_fresnel_analysis.measurement_data['ydata']), max(current_fresnel_analysis.measurement_data['ydata'])],
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

            # Check bounds structure
            if type(current_fresnel_analysis.bounds[0]) is not tuple:
                lower_bound_ = current_fresnel_analysis.bounds[0]
                upper_bound_ = current_fresnel_analysis.bounds[1]
            else:
                lower_bound_ = current_fresnel_analysis.bounds[0][0]
                upper_bound_ = current_fresnel_analysis.bounds[1][0]

            return new_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, current_fresnel_analysis.angle_range, current_fresnel_analysis.ini_guess[0], \
                lower_bound_, upper_bound_, current_fresnel_analysis.extinction_correction, result, current_fresnel_analysis.sensor_object_label, dash.no_update, angle_range_marks, current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 'Data path: \n' + current_fresnel_analysis.initial_data_path, dash.no_update, dash.no_update, current_fresnel_analysis.fit_offset, current_fresnel_analysis.fit_prism_k

    @dash.callback(
        dash.Output('batch-fresnel-analysis-newlayer-radio-selection', 'style'),
        dash.Input('batch-fresnel-analysis-radio-selection', 'value'),
        prevent_initial_call=True)
    def fresnel_modelling_batch_show_newlayer_settings(radio):
        if radio == 0:
            return {'visibility': 'hidden'}
        else:
            return {'visibility': 'visible'}

    @dash.callback(
        dash.Output('batch-fresnel-analysis-example-sensor-dropdown', 'value'),
        dash.Output('batch-fresnel-analysis-example-analysis-dropdown', 'value'),
        dash.Output('batch-fresnel-analysis-background-sensors-dropdown', 'value'),
        dash.Input('batch-fresnel-analysis-button', 'n_clicks'),
        prevent_initial_call=True)
    def batch_clear_dropdowns(button):
        return None, None, None

    @dash.callback(
        dash.Output('fresnel-analysis-manual-prism-k-row', 'style'),
        dash.Output('fresnel-fit-option-extinctionslider', 'value'),
        dash.Input('fresnel-analysis-elastomer-fit', 'value'),
        dash.Input('fresnel-analysis-offset-fit', 'value'),
        dash.State('fresnel-analysis-offset-fit', 'value'),
        dash.State('fresnel-analysis-elastomer-fit', 'value'),
        prevent_initial_call=True)
    def fresnel_modelling_manual_prism_k_hide(in_elastomer_fit_flag, in_offset_fit_flag, state_offset, state_elastomer):
        if in_elastomer_fit_flag or not state_offset:
            return {'margin-bottom': '10px', 'visibility': 'hidden'}, 0
        elif not state_elastomer and in_offset_fit_flag:
            return {'margin-bottom': '10px', 'visibility': 'visible'}, dash.no_update

    @dash.callback(
        dash.Output('fresnel-analysis-elastomer-fit', 'disabled'),
        dash.Output('fresnel-analysis-elastomer-fit', 'value'),
        dash.Input('fresnel-analysis-offset-fit', 'value'),
        prevent_initial_call=True)
    def fresnel_modelling_offset_and_prism_k_boxes(offset_fit_flag):
        if offset_fit_flag:
            return False, dash.no_update
        else:
            return True, False

    @dash.callback(
        dash.Output('batch-fresnel-analysis-background-sensors-dropdown-row', 'style'),
        dash.Output('batch-fresnel-analysis-background-sensors-button-submit', 'style'),
        dash.Output('batch-fresnel-analysis-table', 'style'),
        dash.Input('batch-fresnel-analysis-radio-selection', 'value'),
        prevent_initial_call=True)
    def show_batch_fresnel_analysis_new_layer(radio_value):
        if radio_value == 0:
            return {'visibility': 'hidden'}, {'visibility': 'hidden'}, {'visibility': 'hidden', 'margin-top': '20px'}
        elif radio_value == 1:
            return {'visibility': 'visible'}, {'visibility': 'visible'}, {'visibility': 'visible', 'margin-top': '20px'}

    @dash.callback(
        dash.Output('batch-fresnel-analysis-table', 'children'),
        dash.Output('batch-fresnel-analysis-background-sensors', 'data'),
        dash.Input('batch-fresnel-analysis-background-sensors-button-submit', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-button', 'n_clicks'),
        dash.State('batch-fresnel-analysis-files', 'data'),
        dash.State('batch-fresnel-analysis-background-sensors-dropdown', 'value'),
        prevent_initial_call=True)
    def batch_analysis_submit_background_sensors(submit, batch_analysis_button, measurement_files, background_sensors):
        if 'batch-fresnel-analysis-background-sensors-button-submit' == dash.ctx.triggered_id:
            table_header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Measurement file'), dash.html.Th('Matching background sensor')]))]
            table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td(file_name.split('/')[-1]), dash.html.Td('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name)]) for file_name, sensor_id in zip(measurement_files, background_sensors)])]
            table_children = table_header + table_body
            return table_children, background_sensors

        elif 'batch-fresnel-analysis-button' == dash.ctx.triggered_id:
            table_header = [
                dash.html.Thead(dash.html.Tr([dash.html.Th('Measurement file'), dash.html.Th('Matching background sensor')]))]
            table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td(''), dash.html.Td('')])])]
            table_children = table_header + table_body
            return table_children, None

    @dash.callback(
        dash.Output('batch-fresnel-analysis-finish', 'data'),
        dash.Output('batch-fresnel-spinner', 'spinner_style', allow_duplicate=True),
        dash.Output('batch-fresnel-analysis-modal', 'is_open', allow_duplicate=True),
        dash.Input('batch-fresnel-analysis-done', 'data'),
        prevent_initial_call=True)
    def batch_pass_finish_signal(done_signal):
        if 'batch-fresnel-analysis-done' == dash.ctx.triggered_id:
            return 'finished', {'visibility': 'hidden', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}, False

    @dash.callback(
        dash.Output('batch-fresnel-analysis-modal', 'is_open'),
        dash.Output('batch-fresnel-analysis-example-sensor-dropdown', 'options'),
        dash.Output('batch-fresnel-analysis-background-sensors-dropdown', 'options'),
        dash.Output('batch-fresnel-analysis-example-analysis-dropdown', 'options'),
        dash.Output('batch-fresnel-spinner', 'spinner_style'),
        dash.Output('batch-fresnel-analysis-start', 'data'),
        dash.Input('batch-fresnel-analysis-button', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-confirm', 'n_clicks'),
        dash.Input('batch-fresnel-analysis-cancel', 'n_clicks'),
        prevent_initial_call=True)
    def batch_modal_and_start(analysis_button, confirm_button, cancel_button):

        if 'batch-fresnel-analysis-button' == dash.ctx.triggered_id:
            # Update example sensor and analysis options and open batch modal
            example_sensor_options = [
                {'label': 'S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name, 'value': sensor_id} for
                sensor_id in current_session.sensor_instances]
            example_analysis_options = [{'label': 'FM' + str(fresnel_analysis_id) + ' ' +
                                                  current_session.fresnel_analysis_instances[fresnel_analysis_id].name,
                                         'value': fresnel_analysis_id} for fresnel_analysis_id in
                                        current_session.fresnel_analysis_instances]

            return True, example_sensor_options, example_sensor_options, example_analysis_options, dash.no_update, dash.no_update

        elif 'batch-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, {'visibility': 'visible', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}, 'start'

        elif 'batch-fresnel-analysis-cancel' == dash.ctx.triggered_id:
            return False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        else:
            raise dash.exceptions.PreventUpdate


    @dash.callback(
        dash.Output('exclusion-height-sensorgram-graph', 'figure'),
        dash.Output('result-summary-exclusion-table', 'children', allow_duplicate=True),
        dash.Output('add-exclusion-height-analysis-modal', 'is_open'),
        dash.Output('exclusion-height-analysis-dropdown', 'children'),
        dash.Output('remove-exclusion-height-analysis-modal', 'is_open'),
        dash.Output('exclusion-height-sensor-label', 'children'),
        dash.Output('exclusion-height-fresnel-analysis-label', 'children'),
        dash.Output('exclusion-height-analysis-option-collapse', 'is_open'),
        dash.Output('exclusion-height-progress-collapse', 'is_open', allow_duplicate=True),
        dash.Output('exclusion-height-sensorgram-collapse', 'is_open'),
        dash.Output('exclusion-height-result-collapse', 'is_open', allow_duplicate=True),
        dash.Output('exclusion-height-result-mean-height', 'children', allow_duplicate=True),
        dash.Output('exclusion-height-result-mean-RI', 'children', allow_duplicate=True),
        dash.Output('exclusion-height-result-all-heights', 'children', allow_duplicate=True),
        dash.Output('exclusion-height-result-all-RI', 'children', allow_duplicate=True),
        dash.Output('exclusion-height-SPRvsTIR-graph', 'figure', allow_duplicate=True),
        dash.Output('exclusion-height-reflectivity-graph', 'figure', allow_duplicate=True),
        dash.Output('exclusion-height-d-n-pair-graph', 'figure', allow_duplicate=True),
        dash.Output('exclusion-height-result-pagination', 'max_value'),
        dash.Output('exclusion-height-option-lowerbound', 'value'),
        dash.Output('exclusion-height-option-upperbound', 'value'),
        dash.Output('exclusion-height-settings-injection-points', 'children'),
        dash.Output('exclusion-height-settings-buffer-points', 'children'),
        dash.Output('exclusion-height-settings-probe-points', 'children'),
        dash.Output('exclusion-height-option-resolution', 'value'),
        dash.Output('exclusion-height-analysis-label', 'children'),
        dash.Input('add-exclusion-height-analysis-button', 'n_clicks'),
        dash.Input('add-exclusion-height-analysis-confirm', 'n_clicks'),
        dash.Input({'type': 'exclusion-analysis-list', 'index': dash.ALL}, 'n_clicks'),
        dash.Input('remove-exclusion-height-analysis-button', 'n_clicks'),
        dash.Input('remove-exclusion-height-analysis-confirm', 'n_clicks'),
        dash.Input('remove-exclusion-height-analysis-cancel', 'n_clicks'),
        dash.Input('exclusion-height-sensorgram-graph', 'clickData'),
        dash.Input('exclusion-height-click-action-clear', 'n_clicks'),
        dash.Input('exclusion-height-sensorgram-save-png', 'n_clicks'),
        dash.Input('exclusion-height-sensorgram-save-svg', 'n_clicks'),
        dash.Input('exclusion-height-sensorgram-save-html', 'n_clicks'),
        dash.Input('exclusion-height-sensorgram-save-csv', 'n_clicks'),
        dash.Input('exclusion-height-SPRvsTIR-save-png', 'n_clicks'),
        dash.Input('exclusion-height-SPRvsTIR-save-svg', 'n_clicks'),
        dash.Input('exclusion-height-SPRvsTIR-save-html', 'n_clicks'),
        dash.Input('exclusion-height-SPRvsTIR-save-csv', 'n_clicks'),
        dash.Input('exclusion-height-reflectivity-save-png', 'n_clicks'),
        dash.Input('exclusion-height-reflectivity-save-svg', 'n_clicks'),
        dash.Input('exclusion-height-reflectivity-save-html', 'n_clicks'),
        dash.Input('exclusion-height-reflectivity-save-csv', 'n_clicks'),
        dash.Input('exclusion-height-d-n-pair-save-png', 'n_clicks'),
        dash.Input('exclusion-height-d-n-pair-save-svg', 'n_clicks'),
        dash.Input('exclusion-height-d-n-pair-save-html', 'n_clicks'),
        dash.Input('exclusion-height-d-n-pair-save-csv', 'n_clicks'),
        dash.Input('exclusion-height-result-pagination', 'active_page'),
        dash.Input('exclusion-height-d-n-pair-graph', 'hoverData'),
        dash.Input('exclusion-height-initialize-model', 'n_clicks'),
        dash.State('exclusion-height-analysis-name-input', 'value'),
        dash.State('exclusion-height-click-action-selector', 'value'),
        dash.State('exclusion-choose-background-dropdown', 'value'),
        dash.State('exclusion-height-sensorgram-graph', 'figure'),
        dash.State('exclusion-height-SPRvsTIR-graph', 'figure'),
        dash.State('exclusion-height-reflectivity-graph', 'figure'),
        dash.State('exclusion-height-d-n-pair-graph', 'figure'),
        dash.State('exclusion-height-result-pagination', 'active_page'),
        prevent_initial_call=True)
    def exclusion_height_analysis_control(add_exclusion, confirm_exclusion, choose_exclusion, remove_analysis_button,
                                        remove_confirm, remove_cancel, clickData, clear_points, sensorgram_png,
                                        sensorgram_svg, sensorgram_html, sensorgram_csv, SPRvsTIR_png, SPRvsTIR_svg, SPRvsTIR_html, SPRvsTIR_csv,
                                        reflectivity_save_png, reflectivity_save_svg, reflectivity_save_html,reflectivity_save_csv,
                                        dnpair_save_png, dnpair_save_svg, dnpair_save_html, dnpair_save_csv, active_page, dnpair_hoverdata, initialize_model, analysis_name,
                                        action_selected, background_selected_id, sensorgram_figure_JSON, SPRvsTIR_figure_JSON, reflectivity_figure_JSON,
                                        dnpair_figure_JSON, active_page_state):
        """
        This callback handles what happens when adding new exclusion height objects, choosing different ones, removing them and updating the sensorgram plot with selected probe points etc.
        """

        global current_session
        global current_data_path
        global current_exclusion_height_analysis
        global sensorgram_df_selection
        global ydata_df
        global time_df
        global angles_df

        if 'exclusion-height-sensorgram-graph' == dash.ctx.triggered_id:
            # Determines what happens when clicking on the sensorgram plot

            sensorgram_figure = go.Figure(sensorgram_figure_JSON)

            new_point_index = int(clickData['points'][0]['pointIndex'])
            new_point_time = float(clickData['points'][0]['x'])
            new_point_angle = float(clickData['points'][0]['y'])

            # If the user clicks on the first 20 points but the time is larger than 3 minutes, it is probably a marker click mistake
            if new_point_index < 20 and new_point_time > 3:
                raise dash.exceptions.PreventUpdate

            match action_selected:
                case 1:  # Offset data
                    current_exclusion_height_analysis.sensorgram_offset_ind = new_point_index

                    updated_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                          y=current_exclusion_height_analysis.sensorgram_data[
                                                                'SPR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'SPR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                          name='SPR angle',
                                                          line_color='#636efa'))

                    updated_figure.add_trace(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                        y=current_exclusion_height_analysis.sensorgram_data[
                                                              'TIR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'TIR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                        name='TIR angle',
                                                        line_color='#ef553b'))

                    if len(current_exclusion_height_analysis.injection_points) > 0:
                        injection_points_time = [item[1] for item in current_exclusion_height_analysis.injection_points]
                        injection_points_angle = [item[2] for item in current_exclusion_height_analysis.injection_points]

                        updated_figure.add_trace(go.Scatter(x=injection_points_time,
                                                            y=injection_points_angle,
                                                            name='Injection points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow',
                                                            marker_color='black',
                                                            marker_angle=180,
                                                            showlegend=True))

                    if len(current_exclusion_height_analysis.buffer_points) > 0:
                        buffer_points_time = [item[1] for item in current_exclusion_height_analysis.buffer_points]
                        buffer_points_angle = [item[2] for item in current_exclusion_height_analysis.buffer_points]

                        updated_figure.add_trace(go.Scatter(x=buffer_points_time,
                                                            y=buffer_points_angle,
                                                            name='Buffer points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow',
                                                            showlegend=True))

                    if len(current_exclusion_height_analysis.probe_points) > 0:
                        probe_points_time = [item[1] for item in current_exclusion_height_analysis.probe_points]
                        probe_points_angle = [item[2] for item in current_exclusion_height_analysis.probe_points]

                        updated_figure.add_trace(go.Scatter(x=probe_points_time,
                                                            y=probe_points_angle,
                                                            name='Probe points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow',
                                                            showlegend=True))

                    updated_figure.update_layout(xaxis_title=r'$\large{\text{Time [min]}}$',
                                                 yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                                                 font_family='Balto',
                                                 font_size=19,
                                                 margin_r=25,
                                                 margin_l=60,
                                                 margin_t=40,
                                                 template='simple_white',
                                                 uirevision=True)
                    updated_figure.update_xaxes(mirror=True, showline=True)
                    updated_figure.update_yaxes(mirror=True, showline=True)

                    return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

                case 2:  # Add injection points
                    current_exclusion_height_analysis.injection_points.append((new_point_index, new_point_time, new_point_angle))

                case 3:  # Add buffer points
                    current_exclusion_height_analysis.buffer_points.append((new_point_index, new_point_time, new_point_angle))

                case 4:  # Add probe points
                    current_exclusion_height_analysis.probe_points.append((new_point_index, new_point_time, new_point_angle))

            injection_points_time = [item[1] for item in current_exclusion_height_analysis.injection_points]
            injection_points_angle = [item[2] for item in current_exclusion_height_analysis.injection_points]

            buffer_points_time = [item[1] for item in current_exclusion_height_analysis.buffer_points]
            buffer_points_angle = [item[2] for item in current_exclusion_height_analysis.buffer_points]

            probe_points_time = [item[1] for item in current_exclusion_height_analysis.probe_points]
            probe_points_angle = [item[2] for item in current_exclusion_height_analysis.probe_points]

            updated_figure = go.Figure(go.Scatter(x=sensorgram_figure['data'][0]['x'],
                                                  y=sensorgram_figure['data'][0]['y'],
                                                  name='SPR angle',
                                                  line_color='#636efa'))

            updated_figure.add_trace(go.Scatter(x=sensorgram_figure['data'][1]['x'],
                                                y=sensorgram_figure['data'][1]['y'],
                                                name='TIR angle',
                                                line_color='#ef553b'))

            updated_figure.add_trace(go.Scatter(x=injection_points_time,
                                                y=injection_points_angle,
                                                name='Injection points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                marker_color='black',
                                                marker_angle=180,
                                                showlegend=True))

            updated_figure.add_trace(go.Scatter(x=buffer_points_time,
                                                y=buffer_points_angle,
                                                name='Buffer points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                showlegend=True))

            updated_figure.add_trace(go.Scatter(x=probe_points_time,
                                                y=probe_points_angle,
                                                name='Probe points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                showlegend=True))

            updated_figure.update_layout(xaxis_title=r'$\large{\text{Time [min]}}$',
                                         yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                                         font_family='Balto',
                                         font_size=19,
                                         margin_r=25,
                                         margin_l=60,
                                         margin_t=40,
                                         template='simple_white',
                                         uirevision=True)

            updated_figure.update_xaxes(mirror=True, showline=True)
            updated_figure.update_yaxes(mirror=True, showline=True)

            injection_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.injection_points]
            buffer_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.buffer_points]
            probe_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.probe_points]

            injection_time_string = '{length} injection points: {points}'.format(length=len(injection_points_time_),
                                                                                 points=injection_points_time_)

            buffer_time_string = '{length} buffer points: {points}'.format(length=len(buffer_points_time_),
                                                                           points=buffer_points_time_)

            probe_time_string = '{length} probe points: {points}'.format(length=len(probe_points_time_),
                                                                         points=probe_points_time_)

            current_session.save_session()
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)

            return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, injection_time_string, buffer_time_string, probe_time_string, dash.no_update, dash.no_update

        elif 'exclusion-height-click-action-clear' == dash.ctx.triggered_id:
            # Determines what happens when clearing the selected points (remove from graph and backend object)

            sensorgram_figure = go.Figure(sensorgram_figure_JSON)

            match action_selected:
                case 1:  # Offset data (do nothing)
                    raise dash.exceptions.PreventUpdate

                case 2:  # Clear latest injection point
                    current_exclusion_height_analysis.injection_points = []

                case 3:  # Clear latest buffer point
                    current_exclusion_height_analysis.buffer_points = []

                case 4:  # CLear latest probe point
                    current_exclusion_height_analysis.probe_points = []

            injection_points_time = [item[1] for item in current_exclusion_height_analysis.injection_points]
            injection_points_angle = [item[2] for item in current_exclusion_height_analysis.injection_points]

            buffer_points_time = [item[1] for item in current_exclusion_height_analysis.buffer_points]
            buffer_points_angle = [item[2] for item in current_exclusion_height_analysis.buffer_points]

            probe_points_time = [item[1] for item in current_exclusion_height_analysis.probe_points]
            probe_points_angle = [item[2] for item in current_exclusion_height_analysis.probe_points]

            updated_figure = go.Figure(go.Scatter(x=sensorgram_figure['data'][0]['x'],
                                                  y=sensorgram_figure['data'][0]['y'],
                                                  name='SPR angle',
                                                  line_color='#636efa'))

            updated_figure.add_trace(go.Scatter(x=sensorgram_figure['data'][1]['x'],
                                                y=sensorgram_figure['data'][1]['y'],
                                                name='TIR angle',
                                                line_color='#ef553b'))

            updated_figure.add_trace(go.Scatter(x=injection_points_time,
                                                y=injection_points_angle,
                                                name='Injection points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                marker_color='black',
                                                marker_angle=180,
                                                showlegend=True))

            updated_figure.add_trace(go.Scatter(x=buffer_points_time,
                                                y=buffer_points_angle,
                                                name='Buffer points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                showlegend=True))

            updated_figure.add_trace(go.Scatter(x=probe_points_time,
                                                y=probe_points_angle,
                                                name='Probe points',
                                                mode='markers',
                                                marker_size=14,
                                                marker_symbol='arrow',
                                                showlegend=True))

            updated_figure.update_layout(xaxis_title=r'$\large{\text{Time [min]}}$',
                                         yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                                         font_family='Balto',
                                         font_size=19,
                                         margin_r=25,
                                         margin_l=60,
                                         margin_t=40,
                                         template='simple_white',
                                         uirevision=True)

            updated_figure.update_xaxes(mirror=True, showline=True)
            updated_figure.update_yaxes(mirror=True, showline=True)

            injection_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.injection_points]
            buffer_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.buffer_points]
            probe_points_time_ = [round(item[1], 2) for item in current_exclusion_height_analysis.probe_points]

            injection_time_string = '{length} injection points: {points}'.format(length=len(injection_points_time_),
                                                                                 points=injection_points_time_)

            buffer_time_string = '{length} buffer points: {points}'.format(length=len(buffer_points_time_),
                                                                           points=buffer_points_time_)

            probe_time_string = '{length} probe points: {points}'.format(length=len(probe_points_time_),
                                                                         points=probe_points_time_)
            current_session.save_session()
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)

            return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, injection_time_string, buffer_time_string, probe_time_string, dash.no_update, dash.no_update

        elif 'exclusion-height-d-n-pair-graph' == dash.ctx.triggered_id:

            # Catch the case where the user clicks on the d_n_pair graph before active_page_state has updated
            if isinstance(active_page_state, types.NoneType):
                active_page_state = 1

            # Do not do anything if the d_n_pair_df list contains integers (i.e. is not complete)
            if any([isinstance(obj, int) for obj in current_exclusion_height_analysis.d_n_pair_dfs]):
                raise dash.exceptions.PreventUpdate

            # Calculate fresnel trace for hover data points
            buffer_angles_inj_step = current_exclusion_height_analysis.buffer_reflectivity_dfs[active_page_state-1]['angles'].to_numpy()
            buffer_reflectivity_inj_step = current_exclusion_height_analysis.buffer_reflectivity_dfs[active_page_state-1]['reflectivity'].to_numpy()
            probe_angles_inj_step = current_exclusion_height_analysis.probe_reflectivity_dfs[int((active_page_state-1)/2)]['angles'].to_numpy()
            probe_reflectivity_inj_step = current_exclusion_height_analysis.probe_reflectivity_dfs[int((active_page_state-1)/2)]['reflectivity'].to_numpy()

            buffer_RI_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state-1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'buffer RI']
            probe_RI_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state-1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'probe RI']
            height_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state-1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'height']

            if not current_exclusion_height_analysis.fit_offset:
                buffer_offset_val = current_exclusion_height_analysis.fresnel_object.y_offset
                probe_offset_val = current_exclusion_height_analysis.fresnel_object.y_offset
                buffer_prism_val = current_exclusion_height_analysis.sensor_object.extinction_coefficients[0]
                probe_prism_val = current_exclusion_height_analysis.sensor_object.extinction_coefficients[0]

            elif current_exclusion_height_analysis.fit_offset and not current_exclusion_height_analysis.fit_prism:
                buffer_offset_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state-1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'buffer offsets']
                probe_offset_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state-1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'probe offsets']
                buffer_prism_val = current_exclusion_height_analysis.sensor_object.extinction_coefficients[0]
                probe_prism_val = current_exclusion_height_analysis.sensor_object.extinction_coefficients[0]

            elif current_exclusion_height_analysis.fit_offset and current_exclusion_height_analysis.fit_prism:
                buffer_offset_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state - 1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'buffer offsets']
                probe_offset_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state - 1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'probe offsets']
                buffer_prism_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state - 1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'buffer prism k']
                probe_prism_val = current_exclusion_height_analysis.d_n_pair_dfs[active_page_state - 1].loc[dnpair_hoverdata['points'][0]['pointIndex'], 'probe prism k']

            buffer_ext_coefficients = copy.deepcopy(current_exclusion_height_analysis.sensor_object.extinction_coefficients)
            buffer_ext_coefficients[0] = buffer_prism_val

            probe_ext_coefficients = copy.deepcopy(current_exclusion_height_analysis.sensor_object.extinction_coefficients)
            probe_ext_coefficients[0] = probe_prism_val

            buffer_ref_indices = copy.deepcopy(current_exclusion_height_analysis.sensor_object.refractive_indices)
            buffer_ref_indices[-1] = current_exclusion_height_analysis.buffer_bulk_RIs[active_page_state-1]
            buffer_ref_indices[current_exclusion_height_analysis.sensor_object.fitted_layer_index[0]] = float(buffer_RI_val)

            buffer_layer_thicknesses = copy.deepcopy(current_exclusion_height_analysis.sensor_object.layer_thicknesses)
            buffer_layer_thicknesses[current_exclusion_height_analysis.sensor_object.fitted_layer_index[0]] = float(height_val)

            probe_ref_indices = copy.deepcopy(current_exclusion_height_analysis.sensor_object.refractive_indices)
            probe_ref_indices[-1] = current_exclusion_height_analysis.probe_bulk_RIs[int((active_page_state-1)/2)]
            probe_ref_indices[current_exclusion_height_analysis.sensor_object.fitted_layer_index[0]] = float(
                probe_RI_val)

            probe_layer_thicknesses = copy.deepcopy(current_exclusion_height_analysis.sensor_object.layer_thicknesses)
            probe_layer_thicknesses[current_exclusion_height_analysis.sensor_object.fitted_layer_index[0]] = float(
                height_val)

            buffer_fresnel_coefficients = fresnel_calculation(angles=buffer_angles_inj_step,
                                                              wavelength=current_exclusion_height_analysis.sensor_object.wavelength,
                                                              layer_thicknesses=buffer_layer_thicknesses,
                                                              n_re=buffer_ref_indices,
                                                              n_im=buffer_ext_coefficients,
                                                              ydata_offset=buffer_offset_val,
                                                              )
            probe_fresnel_coefficients = fresnel_calculation(angles=probe_angles_inj_step,
                                                             wavelength=current_exclusion_height_analysis.sensor_object.wavelength,
                                                             layer_thicknesses=probe_layer_thicknesses,
                                                             n_re=probe_ref_indices,
                                                             n_im=probe_ext_coefficients,
                                                             ydata_offset=probe_offset_val,
                                                             )
            
            # Plot mean reflectivity figure with fitted fresnel traces
            mean_reflectivity_figure = go.Figure(
                go.Scatter(x=buffer_angles_inj_step,
                           y=buffer_reflectivity_inj_step,
                           mode='lines',
                           name='Buffer',
                           showlegend=True,
                           line_color='#636EFA'
                           ))

            mean_reflectivity_figure.add_trace(
                go.Scatter(x=probe_angles_inj_step,
                           y=probe_reflectivity_inj_step,
                           mode='lines',
                           name='Probe',
                           showlegend=True,
                           line_color='#EF553B'
                           ))

            mean_reflectivity_figure.add_trace(
                go.Scatter(x=buffer_angles_inj_step,
                           y=buffer_fresnel_coefficients,
                           mode='lines',
                           showlegend=False,
                           line_dash='dash',
                           line_color='black'
                           ))

            mean_reflectivity_figure.add_trace(
                go.Scatter(x=probe_angles_inj_step,
                           y=probe_fresnel_coefficients,
                           mode='lines',
                           showlegend=False,
                           line_dash='dash',
                           line_color='black'
                           ))

            mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                                   yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                                   font_family='Balto',
                                                   font_size=19,
                                                   margin_r=25,
                                                   margin_l=60,
                                                   margin_t=40,
                                                   template='simple_white',
                                                   uirevision=True)
            mean_reflectivity_figure.update_xaxes(mirror=True,
                                                  showline=True)
            mean_reflectivity_figure.update_yaxes(mirror=True,
                                                  showline=True)

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, mean_reflectivity_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-exclusion-height-analysis-button' == dash.ctx.triggered_id:
            # Open add analysis name giving modal
            return dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-exclusion-height-analysis-confirm' == dash.ctx.triggered_id:

            # Add new exclusion height analysis object to session
            background_object = current_session.fresnel_analysis_instances[background_selected_id]
            current_exclusion_height_analysis = add_exclusion_height_object(current_session, background_object, sensorgram_df_selection, current_data_path, analysis_name)
            current_session.save_session()
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)

            # Calculate suggestions of lower and upper bounds for height
            lower_height_bound = float(background_object.sensor_object.layer_thicknesses[-2])
            upper_height_bound = float(background_object.sensor_object.layer_thicknesses[-2]) * 6
            current_exclusion_height_analysis.height_bounds = [lower_height_bound, upper_height_bound]

            analysis_name_ = 'EH' + str(current_exclusion_height_analysis.object_id) + ' ' + current_exclusion_height_analysis.name

            injection_time_string = '0 selected injection points '
            buffer_time_string = '0 selected buffer points '
            probe_time_string = '0 selected probe points '

            # Update choose analysis dropdown menu options
            analysis_options = [dbc.DropdownMenuItem('EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                                                                 id={'type': 'exclusion-analysis-list',
                                                                     'index': exclusion_id},
                                                                 n_clicks=0) for exclusion_id in current_session.exclusion_height_analysis_instances]

            # Update sensorgram graph
            new_sensorgram_fig = go.Figure(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                      y=current_exclusion_height_analysis.sensorgram_data['SPR angle'],
                                                      name='SPR angle',
                                                      line_color='#636efa'))

            new_sensorgram_fig.add_trace(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                    y=current_exclusion_height_analysis.sensorgram_data['TIR angle'],
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

            return new_sensorgram_fig, dash.no_update, False, analysis_options, dash.no_update, background_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, False, True, False, 'Mean exclusion height: None', 'Mean exclusion RI: None', 'All exclusion heights: None', 'All exclusion RI: None', dash.no_update, dash.no_update, dash.no_update, dash.no_update, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, 100, 'Selected analysis: ' + analysis_name_

        elif 'remove-exclusion-height-analysis-button' == dash.ctx.triggered_id:
            # Open remove analysis object confirmation modal
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-exclusion-height-analysis-confirm' == dash.ctx.triggered_id:

            if len(current_session.exclusion_height_analysis_instances) > 1:

                # Pop out the current exclusion height analysis object from the session, delete its .pickle file and make the first instance the current one
                removed = current_exclusion_height_analysis
                current_exclusion_height_analysis = current_session.exclusion_height_analysis_instances[0]
                current_session.remove_exclusion_height_analysis(removed.object_id)
                current_session.save_session()

                # Lower and upper bounds for height
                lower_height_bound = current_exclusion_height_analysis.height_bounds[0]
                upper_height_bound = current_exclusion_height_analysis.height_bounds[1]
                resolution = current_exclusion_height_analysis.d_n_pair_resolution
                analysis_name_ = 'EH' + str(
                    current_exclusion_height_analysis.object_id) + ' ' + current_exclusion_height_analysis.name

                # Update choose analysis dropdown menu options
                analysis_options = [dbc.DropdownMenuItem(
                    'EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                    id={'type': 'exclusion-analysis-list',
                        'index': exclusion_id},
                    n_clicks=0) for exclusion_id in current_session.exclusion_height_analysis_instances]

                # Update results text
                if current_exclusion_height_analysis.mean_exclusion_height_result is not None:
                    mean_result_height = 'Mean exclusion height: {res_h_mean} (std: {res_h_std})'.format(
                        res_h_mean=round(current_exclusion_height_analysis.mean_exclusion_height_result[0], 2),
                        res_h_std=round(current_exclusion_height_analysis.mean_exclusion_height_result[1], 2))
                    mean_result_RI = 'Mean exclusion RI: {res_ri_mean} (std: {res_ri_std})'.format(
                        res_ri_mean=round(current_exclusion_height_analysis.mean_exclusion_RI_result[0], 4),
                        res_ri_std=round(current_exclusion_height_analysis.mean_exclusion_RI_result[1], 4))

                else:
                    mean_result_height = 'Mean exclusion height: None'
                    mean_result_RI = 'Mean exclusion RI: None'

                if current_exclusion_height_analysis.all_exclusion_results is not None:
                    all_result_heights = 'All exclusion heights: {res_h}'.format(
                        res_h=np.round(current_exclusion_height_analysis.all_exclusion_results[0, :], decimals=2))
                    all_result_RI = 'All exclusion RI: {res_RI}'.format(
                        res_RI=np.round(current_exclusion_height_analysis.all_exclusion_results[1, :], decimals=4))
                else:
                    all_result_heights = 'All exclusion heights: None'
                    all_result_RI = 'All exclusion RI: None'

                # Update sensorgram figure to new current exclusion height object sensorgram data, also update points labels
                if current_data_path != current_exclusion_height_analysis.initial_data_path:
                    line_color_value = '#00CC96'
                else:
                    line_color_value = '#636EFA'

                new_sensorgram_fig = go.Figure(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                          y=current_exclusion_height_analysis.sensorgram_data[
                                                              'SPR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'SPR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                          name='SPR angle',
                                                          line_color=line_color_value)
                                               )

                new_sensorgram_fig.add_trace(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                        y=current_exclusion_height_analysis.sensorgram_data[
                                                            'TIR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'TIR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                        name='TIR angle',
                                                        line_color='#ef553b')
                                             )
                # Default point strings if none have been selected
                injection_time_string = '0 selected injection points '
                buffer_time_string = '0 selected buffer points '
                probe_time_string = '0 selected probe points '

                if len(current_exclusion_height_analysis.injection_points) > 0:
                    new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.injection_points],
                                                            y=[item[2] for item in current_exclusion_height_analysis.injection_points],
                                                            name='Injection points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow',
                                                            marker_color='black',
                                                            marker_angle=180)
                                                 )
                    injection_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.injection_points]
                    injection_time_string = '{length} selected injection points: {points}'.format(
                        length=len(injection_points_time),
                        points=injection_points_time)

                if len(current_exclusion_height_analysis.buffer_points) > 0:
                    new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.buffer_points],
                                                            y=[item[2] for item in current_exclusion_height_analysis.buffer_points],
                                                            name='Buffer points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow')
                                                 )
                    buffer_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.buffer_points]
                    buffer_time_string = '{length} selected buffer points: {points}'.format(length=len(buffer_points_time),
                                                                                   points=buffer_points_time)

                if len(current_exclusion_height_analysis.probe_points) > 0:
                    new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.probe_points],
                                                            y=[item[2] for item in current_exclusion_height_analysis.probe_points],
                                                            name='Probe points',
                                                            mode='markers',
                                                            marker_size=14,
                                                            marker_symbol='arrow')
                                                 )
                    probe_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.probe_points]
                    probe_time_string = '{length} selected probe points: {points}'.format(length=len(probe_points_time),
                                                                                 points=probe_points_time)

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

                # Update result figures
                if len(current_exclusion_height_analysis.SPR_vs_TIR_dfs) > 0:
                    SPRvsTIR_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['TIR angles'],
                                                           y=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['SPR angles'],
                                                           mode='lines',
                                                           showlegend=False,
                                                           line_color='#636EFA'
                                                           ))
                else:
                    SPRvsTIR_figure = go.Figure(
                        go.Scatter(x=[0],
                                   y=[0],
                                   mode='lines',
                                   showlegend=False,
                                   line_color='#636EFA'
                                   ))

                SPRvsTIR_figure.update_layout(xaxis_title=r'$\large{\text{TIR angle [ }^{\circ}\text{ ]}}$',
                                              yaxis_title=r'$\large{\text{SPR angle [ }^{\circ}\text{ ]}}$',
                                              font_family='Balto',
                                              font_size=19,
                                              margin_r=25,
                                              margin_l=60,
                                              margin_t=40,
                                              template='simple_white',
                                              uirevision=True)
                SPRvsTIR_figure.update_xaxes(mirror=True,
                                             showline=True)
                SPRvsTIR_figure.update_yaxes(mirror=True,
                                             showline=True)

                if len(current_exclusion_height_analysis.buffer_reflectivity_dfs) > 0:
                    mean_reflectivity_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['angles'],
                                                           y=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['reflectivity'],
                                                           mode='lines',
                                                           name='Buffer',
                                                           showlegend=True,
                                                           line_color='#636EFA'
                                                           ))
                    mean_reflectivity_figure.add_trace(go.Scatter(x=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['angles'],
                                                           y=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['reflectivity'],
                                                           mode='lines',
                                                           name='Probe',
                                                           showlegend=True,
                                                           line_color='#EF553B'
                                                           ))

                else:
                    mean_reflectivity_figure = go.Figure(
                        go.Scatter(x=[0],
                                   y=[0],
                                   mode='lines',
                                   showlegend=False,
                                   line_color='#636EFA'
                                   ))

                mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                              yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                              font_family='Balto',
                                              font_size=19,
                                              margin_r=25,
                                              margin_l=60,
                                              margin_t=40,
                                              template='simple_white',
                                              uirevision=True)
                mean_reflectivity_figure.update_xaxes(mirror=True,
                                             showline=True)
                mean_reflectivity_figure.update_yaxes(mirror=True,
                                             showline=True)

                if (len(current_exclusion_height_analysis.d_n_pair_dfs) > 0) and not any([isinstance(obj, int) for obj in current_exclusion_height_analysis.d_n_pair_dfs]):

                    d_n_pair_figure = go.Figure(go.Scatter(
                        x=current_exclusion_height_analysis.d_n_pair_dfs[0]['buffer RI'],
                        y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                        mode='lines',
                        name='Buffer',
                        showlegend=True,
                        line_color='#636EFA'
                    ))
                    d_n_pair_figure.add_trace(go.Scatter(
                        x=current_exclusion_height_analysis.d_n_pair_dfs[0]['probe RI'],
                        y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                        mode='lines',
                        name='Probe',
                        showlegend=True,
                        line_color='#EF553B'
                    ))

                else:
                    d_n_pair_figure = go.Figure(
                        go.Scatter(x=[0],
                                   y=[0],
                                   mode='lines',
                                   showlegend=False,
                                   line_color='#636EFA'
                                   ))

                d_n_pair_figure.update_layout(
                    xaxis_title=r'$\large{\text{Refractive index}}$',
                    yaxis_title=r'$\large{\text{Height [nm]}}$',
                    font_family='Balto',
                    font_size=19,
                    margin_r=25,
                    margin_l=60,
                    margin_t=40,
                    template='simple_white',
                    uirevision=True)
                d_n_pair_figure.update_xaxes(mirror=True,
                                                      showline=True)
                d_n_pair_figure.update_yaxes(mirror=True,
                                                      showline=True)

                # Update number of injection steps in pagination of result page
                if len(current_exclusion_height_analysis.probe_points) > 0:
                    num_injection_steps = len(current_exclusion_height_analysis.probe_points)
                else:
                    num_injection_steps = dash.no_update

                table_header = [dash.html.Thead(
                    dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean [all]'), dash.html.Th('Exclusion RI mean [all]')]))]
                table_body = [dash.html.Tbody([dash.html.Tr(
                    [dash.html.Td(current_session.exclusion_height_analysis_instances[exclusion_inst].name),
                     dash.html.Td('{mean_} {all_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[0], 2), all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2)))),
                     dash.html.Td('{mean_} {all_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0], 4), all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :], decimals=4))))])
                    for exclusion_inst in current_session.exclusion_height_analysis_instances])]
                exclusion_result_summary_dataframe = table_header + table_body

                return new_sensorgram_fig, exclusion_result_summary_dataframe, False, analysis_options, False, current_exclusion_height_analysis.fresnel_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, True, True, True, mean_result_height, mean_result_RI, all_result_heights, all_result_RI, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, num_injection_steps, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, resolution, 'Selected analysis: ' + analysis_name_

            else:
                try:
                    current_session.remove_exclusion_height_analysis(current_exclusion_height_analysis.object_id)
                except AttributeError:
                    pass  # There was no object at all

                current_exclusion_height_analysis = None
                current_session.save_session()

                table_header = [dash.html.Thead(
                    dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean [all]'),
                                  dash.html.Th('Exclusion RI mean [all]')]))]
                table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td(''), dash.html.Td(''), dash.html.Td('')])])]
                exclusion_result_summary_dataframe = table_header + table_body

                return dash.no_update, exclusion_result_summary_dataframe, dash.no_update, dash.no_update, False, 'Sensor: None', 'Fresnel background: None', False, False, False, False, 'Mean exclusion height: None', 'Mean exclusion RI: None', 'All exclusion heights: None', 'All exclusion RI: None', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, '', '', '', dash.no_update, dash.no_update

        elif 'remove-exclusion-height-analysis-cancel' == dash.ctx.triggered_id:
            # Cancel removal of exclusion height analysis object

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'exclusion-height-SPRvsTIR-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(SPRvsTIR_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(SPRvsTIR_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(SPRvsTIR_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            SPRvsTIR_fig = go.Figure(SPRvsTIR_figure_JSON)
            fig_df = pd.DataFrame(data={'TIR': list(SPRvsTIR_fig.data[0].x["_inputArray"].values())[:-3], 'SPR': list(SPRvsTIR_fig.data[0].y["_inputArray"].values())[:-3]})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(reflectivity_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(reflectivity_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(reflectivity_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            reflectivity_fig = go.Figure(reflectivity_figure_JSON)
            fig_keys_x = ['x' + str(i) for i in range(len(reflectivity_fig.data))]
            fig_keys_y = ['y' + str(i) for i in range(len(reflectivity_fig.data))]
            fig_keys = [key for sublist in zip(fig_keys_x, fig_keys_y) for key in sublist]
            fig_values_x = []
            for i in range(len(reflectivity_fig.data)):
                fig_values_x.append(list(reflectivity_fig.data[i].x["_inputArray"].values())[:-3])
            fig_values_y = []
            for i in range(len(reflectivity_fig.data)):
                fig_values_y.append(list(reflectivity_fig.data[i].y["_inputArray"].values())[:-3])
            fig_values = [value for sublist in zip(fig_values_x, fig_values_y) for value in sublist]
            fig_df = pd.DataFrame(data={key:value for (key, value) in zip(fig_keys, fig_values)})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(sensorgram_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(sensorgram_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(sensorgram_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            sensorgram_fig = go.Figure(sensorgram_figure_JSON)
            fig_keys_x = ['x' + str(i) for i in range(len(sensorgram_fig.data))]
            fig_keys_y = ['y' + str(i) for i in range(len(sensorgram_fig.data))]
            fig_keys_excl = [key for sublist in zip(fig_keys_x, fig_keys_y) for key in sublist]
            fig_values_x = []
            for i in range(len(sensorgram_fig.data)):
                if not type(sensorgram_fig.data[i].x) == tuple:
                    fig_values_x.append(list(sensorgram_fig.data[i].x["_inputArray"].values())[:-3])
                else:
                    fig_values_x.append(list(sensorgram_fig.data[i].x))
            fig_values_y = []
            for i in range(len(sensorgram_fig.data)):
                if not type(sensorgram_fig.data[i].y) == tuple:
                    fig_values_y.append(list(sensorgram_fig.data[i].y["_inputArray"].values())[:-3])
                else:
                    fig_values_y.append(list(sensorgram_fig.data[i].y))
            fig_values_excl = [value for sublist in zip(fig_values_x, fig_values_y) for value in sublist]
            for i in range(len(fig_values_excl)):
                if len(fig_values_excl[i]) < len(fig_values_excl[0]):
                    fig_values_excl[i] = fig_values_excl[i]+[0]*(len(fig_values_excl[0])-len(fig_values_excl[i]))
            fig_df = pd.DataFrame(data={key: value for (key, value) in zip(fig_keys_excl, fig_values_excl)})
            fig_df.to_csv(save_filename, sep=';')

            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(dnpair_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(dnpair_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(),
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(dnpair_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-csv' == dash.ctx.triggered_id:
            d_n_pair_fig_csv = go.Figure(dnpair_figure_JSON)
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')
            fig_df = pd.DataFrame(data={'n_buffer': list(d_n_pair_fig_csv.data[0].x["_inputArray"].values())[:-3], 'd_buffer': list(d_n_pair_fig_csv.data[0].y["_inputArray"].values())[:-3], 'n_probe': list(d_n_pair_fig_csv.data[1].x["_inputArray"].values())[:-3], 'd_probe': list(d_n_pair_fig_csv.data[1].y["_inputArray"].values())[:-3]})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-result-pagination' == dash.ctx.triggered_id:

            SPR_TIR_probe_ind = int((active_page-1)/2)

            # Update result figures
            SPRvsTIR_figure = go.Figure(
                go.Scatter(x=current_exclusion_height_analysis.SPR_vs_TIR_dfs[SPR_TIR_probe_ind]['TIR angles'],
                           y=current_exclusion_height_analysis.SPR_vs_TIR_dfs[SPR_TIR_probe_ind]['SPR angles'],
                           mode='lines',
                           showlegend=False,
                           line_color='#636EFA'
                           ))

            SPRvsTIR_figure.update_layout(xaxis_title=r'$\large{\text{TIR angle [ }^{\circ}\text{ ]}}$',
                                          yaxis_title=r'$\large{\text{SPR angle [ }^{\circ}\text{ ]}}$',
                                          font_family='Balto',
                                          font_size=19,
                                          margin_r=25,
                                          margin_l=60,
                                          margin_t=40,
                                          template='simple_white',
                                          uirevision=True)
            SPRvsTIR_figure.update_xaxes(mirror=True,
                                         showline=True)
            SPRvsTIR_figure.update_yaxes(mirror=True,
                                         showline=True)

            mean_reflectivity_figure = go.Figure(
                go.Scatter(x=current_exclusion_height_analysis.buffer_reflectivity_dfs[active_page-1]['angles'],
                           y=current_exclusion_height_analysis.buffer_reflectivity_dfs[active_page-1]['reflectivity'],
                           mode='lines',
                           name='Buffer',
                           showlegend=True,
                           line_color='#636EFA'
                           ))
            mean_reflectivity_figure.add_trace(
                go.Scatter(x=current_exclusion_height_analysis.probe_reflectivity_dfs[SPR_TIR_probe_ind]['angles'],
                           y=current_exclusion_height_analysis.probe_reflectivity_dfs[SPR_TIR_probe_ind]['reflectivity'],
                           mode='lines',
                           name='Probe',
                           showlegend=True,
                           line_color='#EF553B'
                           ))

            mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                                   yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                                   font_family='Balto',
                                                   font_size=19,
                                                   margin_r=25,
                                                   margin_l=60,
                                                   margin_t=40,
                                                   template='simple_white',
                                                   uirevision=True)
            mean_reflectivity_figure.update_xaxes(mirror=True,
                                                  showline=True)
            mean_reflectivity_figure.update_yaxes(mirror=True,
                                                  showline=True)

            if (len(current_exclusion_height_analysis.d_n_pair_dfs) > 0) and not any([isinstance(obj, int) for obj in current_exclusion_height_analysis.d_n_pair_dfs]):
                d_n_pair_figure = go.Figure(go.Scatter(
                    x=current_exclusion_height_analysis.d_n_pair_dfs[active_page-1]['buffer RI'],
                    y=current_exclusion_height_analysis.d_n_pair_dfs[active_page-1]['height'],
                    mode='lines',
                    name='Buffer',
                    showlegend=True,
                    line_color='#636EFA'
                ))
                d_n_pair_figure.add_trace(go.Scatter(
                    x=current_exclusion_height_analysis.d_n_pair_dfs[active_page-1]['probe RI'],
                    y=current_exclusion_height_analysis.d_n_pair_dfs[active_page-1]['height'],
                    mode='lines',
                    name='Probe',
                    showlegend=True,
                    line_color='#EF553B'
                ))
                d_n_pair_figure.update_layout(
                    xaxis_title=r'$\large{\text{Refractive index}}$',
                    yaxis_title=r'$\large{\text{Height [nm]}}$',
                    font_family='Balto',
                    font_size=19,
                    margin_r=25,
                    margin_l=60,
                    margin_t=40,
                    template='simple_white',
                    uirevision=True)

                d_n_pair_figure.update_xaxes(mirror=True,
                                             showline=True)
                d_n_pair_figure.update_yaxes(mirror=True,
                                             showline=True)
            else:
                d_n_pair_figure = dash.no_update

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'exclusion-height-initialize-model' == dash.ctx.triggered_id:

            # Check that appropriate points have been selected
            if len(current_exclusion_height_analysis.buffer_points) % 4 == 0 and len(current_exclusion_height_analysis.injection_points) % 2 == 0 and len(current_exclusion_height_analysis.probe_points) % 2 == 0:

                # Initializes model parameters and attributes prepping for running. Also activate run buttons and result page.
                current_exclusion_height_analysis.initialize_model(ydata_df)

                current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)
                current_session.save_session()

                SPRvsTIR_figure = go.Figure(
                    go.Scatter(x=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['TIR angles'],
                               y=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['SPR angles'],
                               mode='lines',
                               showlegend=False,
                               line_color='#636EFA'
                               ))

                SPRvsTIR_figure.update_layout(xaxis_title=r'$\large{\text{TIR angle [ }^{\circ}\text{ ]}}$',
                                              yaxis_title=r'$\large{\text{SPR angle [ }^{\circ}\text{ ]}}$',
                                              font_family='Balto',
                                              font_size=19,
                                              margin_r=25,
                                              margin_l=60,
                                              margin_t=40,
                                              template='simple_white',
                                              uirevision=True)
                SPRvsTIR_figure.update_xaxes(mirror=True,
                                             showline=True)
                SPRvsTIR_figure.update_yaxes(mirror=True,
                                             showline=True)

                mean_reflectivity_figure = go.Figure(
                    go.Scatter(x=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['angles'],
                               y=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['reflectivity'],
                               mode='lines',
                               name='Buffer',
                               showlegend=True,
                               line_color='#636EFA'
                               ))
                mean_reflectivity_figure.add_trace(
                    go.Scatter(x=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['angles'],
                               y=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['reflectivity'],
                               mode='lines',
                               name='Probe',
                               showlegend=True,
                               line_color='#EF553B'
                               ))

                mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                                       yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                                       font_family='Balto',
                                                       font_size=19,
                                                       margin_r=25,
                                                       margin_l=60,
                                                       margin_t=40,
                                                       template='simple_white',
                                                       uirevision=True)
                mean_reflectivity_figure.update_xaxes(mirror=True,
                                                      showline=True)
                mean_reflectivity_figure.update_yaxes(mirror=True,
                                                      showline=True)

                # Update number of injection steps in pagination of result page
                num_injection_steps = len(current_exclusion_height_analysis.probe_points)

                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, SPRvsTIR_figure, mean_reflectivity_figure, dash.no_update, num_injection_steps, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

            else:
                if len(current_exclusion_height_analysis.injection_points) % 2 != 0:
                    print('ERROR: Odd number of selected injection points. Need to select 2 points per injection.')
                elif len(current_exclusion_height_analysis.buffer_points) % 4 != 0:
                    print('ERROR: Wrong number of selected buffer points. Need to select 4 points per injection.')
                elif len(current_exclusion_height_analysis.probe_points) % 2 != 0:
                    print('ERROR: Odd number of selected probe points. Need to select 2 points per injection.')

                raise dash.exceptions.PreventUpdate

        else:
            # Selecting a previously existing analysis object from pattern matching callbacks

            # Select a new current exclusion height analysis object
            current_exclusion_height_analysis = current_session.exclusion_height_analysis_instances[dash.callback_context.triggered_id.index]

            # Lower and upper bounds for height
            lower_height_bound = current_exclusion_height_analysis.height_bounds[0]
            upper_height_bound = current_exclusion_height_analysis.height_bounds[1]
            resolution = current_exclusion_height_analysis.d_n_pair_resolution

            analysis_name_ = 'EH' + str(
                current_exclusion_height_analysis.object_id) + ' ' + current_exclusion_height_analysis.name

            # Update choose analysis dropdown menu options
            analysis_options = [dbc.DropdownMenuItem(
                'EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                id={'type': 'exclusion-analysis-list',
                    'index': exclusion_id},
                n_clicks=0) for exclusion_id in current_session.exclusion_height_analysis_instances]

            # Update results text
            if current_exclusion_height_analysis.mean_exclusion_height_result is not None:
                mean_result_height = 'Mean exclusion height: {res_h_mean} (std: {res_h_std})'.format(
                    res_h_mean=round(current_exclusion_height_analysis.mean_exclusion_height_result[0], 2),
                    res_h_std=round(current_exclusion_height_analysis.mean_exclusion_height_result[1], 2))
                mean_result_RI = 'Mean exclusion RI: {res_ri_mean} (std: {res_ri_std})'.format(
                    res_ri_mean=round(current_exclusion_height_analysis.mean_exclusion_RI_result[0], 4),
                    res_ri_std=round(current_exclusion_height_analysis.mean_exclusion_RI_result[1], 4))

            else:
                mean_result_height = 'Mean exclusion height: None'
                mean_result_RI = 'Mean exclusion RI: None'

            if current_exclusion_height_analysis.all_exclusion_results is not None:
                all_result_heights = 'All exclusion heights: {res_h}'.format(
                    res_h=np.round(current_exclusion_height_analysis.all_exclusion_results[0, :], decimals=2))
                all_result_RI = 'All exclusion RI: {res_RI}'.format(
                    res_RI=np.round(current_exclusion_height_analysis.all_exclusion_results[1, :], decimals=4))
            else:
                all_result_heights = 'All exclusion heights: None'
                all_result_RI = 'All exclusion RI: None'

            # Update sensorgram figure to new current exclusion height object sensorgram data
            if current_data_path != current_exclusion_height_analysis.initial_data_path:
                line_color_value = '#00CC96'
            else:
                line_color_value = '#636EFA'

            new_sensorgram_fig = go.Figure(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                      y=current_exclusion_height_analysis.sensorgram_data[
                                                          'SPR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'SPR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                      name='SPR angle',
                                                      line_color=line_color_value)
                                           )

            new_sensorgram_fig.add_trace(go.Scatter(x=current_exclusion_height_analysis.sensorgram_data['time'],
                                                    y=current_exclusion_height_analysis.sensorgram_data[
                                                        'TIR angle'] - current_exclusion_height_analysis.sensorgram_data[
                                                                'TIR angle'].iloc[current_exclusion_height_analysis.sensorgram_offset_ind],
                                                    name='TIR angle',
                                                    line_color='#ef553b')
                                         )

            # Default points string
            injection_time_string = '0 selected injection points'
            buffer_time_string = '0 selected buffer points'
            probe_time_string = '0 selected probe points'

            if len(current_exclusion_height_analysis.injection_points) > 0:
                new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.injection_points],
                                                        y=[item[2] for item in current_exclusion_height_analysis.injection_points],
                                                        name='Injection points',
                                                        mode='markers',
                                                        marker_size=14,
                                                        marker_symbol='arrow',
                                                        marker_color='black',
                                                        marker_angle=180)
                                             )
                injection_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.injection_points]
                injection_time_string = '{length} selected injection points: {points}'.format(
                    length=len(injection_points_time),
                    points=injection_points_time)

            if len(current_exclusion_height_analysis.buffer_points) > 0:
                new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.buffer_points],
                                                        y=[item[2] for item in current_exclusion_height_analysis.buffer_points],
                                                        name='Buffer points',
                                                        mode='markers',
                                                        marker_size=14,
                                                        marker_symbol='arrow')
                                             )
                buffer_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.buffer_points]
                buffer_time_string = '{length} selected buffer points: {points}'.format(length=len(buffer_points_time),
                                                                                        points=buffer_points_time)

            if len(current_exclusion_height_analysis.probe_points) > 0:
                new_sensorgram_fig.add_trace(go.Scatter(x=[item[1] for item in current_exclusion_height_analysis.probe_points],
                                                        y=[item[2] for item in current_exclusion_height_analysis.probe_points],
                                                        name='Probe points',
                                                        mode='markers',
                                                        marker_size=14,
                                                        marker_symbol='arrow')
                                             )
                probe_points_time = [round(item[1], 2) for item in current_exclusion_height_analysis.probe_points]
                probe_time_string = '{length} selected probe points: {points}'.format(length=len(probe_points_time),
                                                                                      points=probe_points_time)

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

            # Update result figures
            if len(current_exclusion_height_analysis.SPR_vs_TIR_dfs) > 0:
                SPRvsTIR_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['TIR angles'],
                                                       y=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['SPR angles'],
                                                       mode='lines',
                                                       showlegend=False,
                                                       line_color='#636EFA'
                                                       ))
            else:
                SPRvsTIR_figure = go.Figure(
                    go.Scatter(x=[0],
                               y=[0],
                               mode='lines',
                               showlegend=False,
                               line_color='#636EFA'
                               ))

            SPRvsTIR_figure.update_layout(xaxis_title=r'$\large{\text{TIR angle [ }^{\circ}\text{ ]}}$',
                                          yaxis_title=r'$\large{\text{SPR angle [ }^{\circ}\text{ ]}}$',
                                          font_family='Balto',
                                          font_size=19,
                                          margin_r=25,
                                          margin_l=60,
                                          margin_t=40,
                                          template='simple_white',
                                          uirevision=True)
            SPRvsTIR_figure.update_xaxes(mirror=True,
                                         showline=True)
            SPRvsTIR_figure.update_yaxes(mirror=True,
                                         showline=True)

            if len(current_exclusion_height_analysis.buffer_reflectivity_dfs) > 0:
                mean_reflectivity_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['angles'],
                                                       y=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['reflectivity'],
                                                       mode='lines',
                                                       name='Buffer',
                                                       showlegend=True,
                                                       line_color='#636EFA'
                                                       ))
                mean_reflectivity_figure.add_trace(go.Scatter(x=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['angles'],
                                                       y=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['reflectivity'],
                                                       mode='lines',
                                                       name='Probe',
                                                       showlegend=True,
                                                       line_color='#EF553B'
                                                       ))
            else:
                mean_reflectivity_figure = go.Figure(
                    go.Scatter(x=[0],
                               y=[0],
                               mode='lines',
                               showlegend=False,
                               line_color='#636EFA'
                               ))

            mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                          yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                          font_family='Balto',
                                          font_size=19,
                                          margin_r=25,
                                          margin_l=60,
                                          margin_t=40,
                                          template='simple_white',
                                          uirevision=True)
            mean_reflectivity_figure.update_xaxes(mirror=True,
                                         showline=True)
            mean_reflectivity_figure.update_yaxes(mirror=True,
                                         showline=True)

            if (len(current_exclusion_height_analysis.d_n_pair_dfs) > 0) and not any([isinstance(obj, int) for obj in current_exclusion_height_analysis.d_n_pair_dfs]):

                d_n_pair_figure = go.Figure(go.Scatter(
                    x=current_exclusion_height_analysis.d_n_pair_dfs[0]['buffer RI'],
                    y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                    mode='lines',
                    name='Buffer',
                    showlegend=True,
                    line_color='#636EFA'
                ))
                d_n_pair_figure.add_trace(go.Scatter(
                    x=current_exclusion_height_analysis.d_n_pair_dfs[0]['probe RI'],
                    y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                    mode='lines',
                    name='Probe',
                    showlegend=True,
                    line_color='#EF553B'
                ))

            else:
                d_n_pair_figure = go.Figure(
                    go.Scatter(x=[0],
                               y=[0],
                               mode='lines',
                               showlegend=False,
                               line_color='#636EFA'
                               ))

            d_n_pair_figure.update_layout(
                xaxis_title=r'$\large{\text{Refractive index}}$',
                yaxis_title=r'$\large{\text{Height [nm]}}$',
                font_family='Balto',
                font_size=19,
                margin_r=25,
                margin_l=60,
                margin_t=40,
                template='simple_white',
                uirevision=True)
            d_n_pair_figure.update_xaxes(mirror=True,
                                                  showline=True)
            d_n_pair_figure.update_yaxes(mirror=True,
                                                  showline=True)

            # Update number of injection steps in pagination of result page
            num_injection_steps = len(current_exclusion_height_analysis.probe_points)

            return new_sensorgram_fig, dash.no_update, False, analysis_options, False, current_exclusion_height_analysis.fresnel_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, False, True, True, mean_result_height, mean_result_RI, all_result_heights, all_result_RI, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, num_injection_steps, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, resolution, 'Selected analysis: ' + analysis_name_

    @dash.callback(
        dash.Output('exclusion-height-result-collapse', 'is_open'),
        dash.Output('result-summary-exclusion-table', 'children'),
        dash.Output('exclusion-height-result-mean-height', 'children'),
        dash.Output('exclusion-height-result-mean-RI', 'children'),
        dash.Output('exclusion-height-result-all-heights', 'children'),
        dash.Output('exclusion-height-result-all-RI', 'children'),
        dash.Output('exclusion-height-SPRvsTIR-graph', 'figure'),
        dash.Output('exclusion-height-reflectivity-graph', 'figure'),
        dash.Output('exclusion-height-d-n-pair-graph', 'figure'),
        dash.Output('exclusion-height-spinner', 'spinner_style', allow_duplicate=True),
        dash.Output('exclusion-height-abort-button', 'disabled', allow_duplicate=True),
        dash.Input('exclusion-height-run-button', 'n_clicks'),
        dash.State('exclusion-height-option-lowerbound', 'value'),
        dash.State('exclusion-height-option-upperbound', 'value'),
        dash.State('exclusion-height-option-resolution', 'value'),
        dash.State('exclusion-height-analysis-offset-refit', 'value'),
        dash.State('exclusion-height-analysis-prism-refit', 'value'),
        prevent_initial_call=True)
    def run_exclusion_height_calculations(run_button, lower_bound, upper_bound, resolution, fit_offset_flag, fit_prism_flag):
        """
        This callback runs the exclusion height calculations in the background. It is triggered by the run button, and
        updates the progress bar. It also updates the result figures when the calculations are done.

        :param set_progress: defines the progress bar value and max value
        :param run_button:
        :param check_button:
        :param abort_button:
        :param lower_bound:
        :param upper_bound:
        :param resolution:
        :return:
        """

        global current_session
        global current_exclusion_height_analysis

        if 'exclusion-height-run-button' == dash.ctx.triggered_id:

            # Set flags
            current_exclusion_height_analysis.abort_flag = False
            current_exclusion_height_analysis.fit_offset = fit_offset_flag
            current_exclusion_height_analysis.fit_prism = fit_prism_flag

            # Set resolution and height steps
            current_exclusion_height_analysis.d_n_pair_resolution = resolution
            current_exclusion_height_analysis.height_bounds[0] = lower_bound
            current_exclusion_height_analysis.height_bounds[1] = upper_bound
            current_exclusion_height_analysis.height_steps = np.linspace(lower_bound, upper_bound, resolution)

            # Overwrite previous results
            current_exclusion_height_analysis.all_exclusion_results = np.zeros((2, len(current_exclusion_height_analysis.injection_points)))
            current_exclusion_height_analysis.d_n_pair_dfs = [0] * len(current_exclusion_height_analysis.injection_points)
            current_exclusion_height_analysis.mean_exclusion_height_result = None
            current_exclusion_height_analysis.mean_exclusion_RI_result = None

            # Run exclusion height calculations
            process_all_exclusion_heights(current_exclusion_height_analysis, logical_cores)

            # Wait for all results to be in
            while 0 in current_exclusion_height_analysis.all_exclusion_results:
                time.sleep(3)

            # Calculate mean exclusion height and RI, along with standard deviation (as a tuple)
            current_exclusion_height_analysis.mean_exclusion_height_result = (np.nanmean(current_exclusion_height_analysis.all_exclusion_results[0, :]), np.nanstd(current_exclusion_height_analysis.all_exclusion_results[0, :]))
            current_exclusion_height_analysis.mean_exclusion_RI_result = (np.nanmean(current_exclusion_height_analysis.all_exclusion_results[1, :]), np.nanstd(current_exclusion_height_analysis.all_exclusion_results[1, :]))

            mean_result_height = 'Mean exclusion height: {res_h_mean} (std: {res_h_std})'.format(
                res_h_mean=round(current_exclusion_height_analysis.mean_exclusion_height_result[0], 2),
                res_h_std=round(current_exclusion_height_analysis.mean_exclusion_height_result[1], 2))
            mean_result_RI = 'Mean exclusion RI: {res_ri_mean} (std: {res_ri_std})'.format(
                res_ri_mean=round(current_exclusion_height_analysis.mean_exclusion_RI_result[0], 4),
                res_ri_std=round(current_exclusion_height_analysis.mean_exclusion_RI_result[1], 4))

            all_result_heights = 'All exclusion heights: {res_h}'.format(res_h=np.round(current_exclusion_height_analysis.all_exclusion_results[0, :], decimals=2))
            all_result_RI = 'All exclusion RI: {res_RI}'.format(res_RI=np.round(current_exclusion_height_analysis.all_exclusion_results[1, :], decimals=4))

            # Save session
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)
            current_session.save_session()

            # Update result figures
            SPRvsTIR_figure = go.Figure(go.Scatter(x=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['TIR angles'],
                                                   y=current_exclusion_height_analysis.SPR_vs_TIR_dfs[0]['SPR angles'],
                                                   mode='lines',
                                                   showlegend=False,
                                                   line_color='#636EFA'
                                                   ))
            SPRvsTIR_figure.update_layout(xaxis_title=r'$\large{\text{TIR angle [ }^{\circ}\text{ ]}}$',
                                            yaxis_title=r'$\large{\text{SPR angle [ }^{\circ}\text{ ]}}$',
                                            font_family='Balto',
                                            font_size=19,
                                            margin_r=25,
                                            margin_l=60,
                                            margin_t=40,
                                            template='simple_white',
                                            uirevision=True)
            SPRvsTIR_figure.update_xaxes(mirror=True,
                                            showline=True)
            SPRvsTIR_figure.update_yaxes(mirror=True,
                                            showline=True)

            mean_reflectivity_figure = go.Figure(
                go.Scatter(x=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['angles'],
                           y=current_exclusion_height_analysis.buffer_reflectivity_dfs[0]['reflectivity'],
                           mode='lines',
                           name='Buffer',
                           showlegend=True,
                           line_color='#636EFA'
                           ))
            mean_reflectivity_figure.add_trace(
                go.Scatter(x=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['angles'],
                           y=current_exclusion_height_analysis.probe_reflectivity_dfs[0]['reflectivity'],
                           mode='lines',
                           name='Probe',
                           showlegend=True,
                           line_color='#EF553B'
                           ))
            mean_reflectivity_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                                   yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                                   font_family='Balto',
                                                   font_size=19,
                                                   margin_r=25,
                                                   margin_l=60,
                                                   margin_t=40,
                                                   template='simple_white',
                                                   uirevision=True)
            mean_reflectivity_figure.update_xaxes(mirror=True,
                                                  showline=True)
            mean_reflectivity_figure.update_yaxes(mirror=True,
                                                  showline=True)

            d_n_pair_figure = go.Figure(go.Scatter(
                x=current_exclusion_height_analysis.d_n_pair_dfs[0]['buffer RI'],
                y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                mode='lines',
                name='Buffer',
                showlegend=True,
                line_color='#636EFA'
            ))
            d_n_pair_figure.add_trace(go.Scatter(
                x=current_exclusion_height_analysis.d_n_pair_dfs[0]['probe RI'],
                y=current_exclusion_height_analysis.d_n_pair_dfs[0]['height'],
                mode='lines',
                name='Probe',
                showlegend=True,
                line_color='#EF553B'
            ))

            d_n_pair_figure.update_layout(
                xaxis_title=r'$\large{\text{Refractive index}}$',
                yaxis_title=r'$\large{\text{Height [nm]}}$',
                font_family='Balto',
                font_size=19,
                margin_r=25,
                margin_l=60,
                margin_t=40,
                template='simple_white',
                uirevision=True)
            d_n_pair_figure.update_xaxes(mirror=True,
                                         showline=True)
            d_n_pair_figure.update_yaxes(mirror=True,
                                         showline=True)

            table_header = [dash.html.Thead(
                dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean [all]'),
                              dash.html.Th('Exclusion RI mean [all]')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td(current_session.exclusion_height_analysis_instances[exclusion_inst].name),
                 dash.html.Td('{mean_} {all_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[
                         0], 2), all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :],
                     decimals=2)))),
                 dash.html.Td('{mean_} {all_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0],
                     4), all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :],
                     decimals=4))))])
                for exclusion_inst in current_session.exclusion_height_analysis_instances])]
            exclusion_result_summary_dataframe = table_header + table_body

            return True, exclusion_result_summary_dataframe, mean_result_height, mean_result_RI, all_result_heights, all_result_RI, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, {'visibility': 'hidden', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}, True

    @dash.callback(
        dash.Output('exclusion-height-analysis-prism-refit', 'value'),
        dash.Output('exclusion-height-analysis-prism-refit', 'disabled'),
        dash.Input('exclusion-height-analysis-offset-refit', 'value'),
        prevent_initial_call=True)
    def check_exclusion_prism_fitting_checkbox(offset_value):
        if offset_value:
            return dash.no_update, False
        else:
            return False, True

    @dash.callback(
        dash.Output('exclusion-height-spinner', 'spinner_style'),
        dash.Output('exclusion-height-abort-button', 'disabled'),
        dash.Input('exclusion-height-run-button', 'n_clicks'),
        dash.Input('exclusion-height-abort-button', 'n_clicks'),
        prevent_initial_call=True)
    def exclusion_calculation_activate_spinner(run_button, abort_button):
        global current_exclusion_height_analysis

        if 'exclusion-height-run-button' == dash.ctx.triggered_id:
            return {'visibility': 'visible', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}, False

        elif 'exclusion-height-abort-button' == dash.ctx.triggered_id:
            current_exclusion_height_analysis.abort_flag = True

            return {'visibility': 'hidden', 'margin-top': '10px', 'margin-right': '10px', 'width': '2rem', 'height': '2rem'}, True

    @dash.callback(
        dash.Input('export-single-file-button', 'n_clicks'),
        prevent_initial_call=True)
    def export_result_summary(single_file):
        global current_session

        if 'export-single-file-button' == dash.ctx.triggered_id:
            fresnel_df = pd.DataFrame(
                                        {'Analysis': [current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Variable': ['{layer}|{parameter}-{channel}'.format(
                                             layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                                             parameter=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                                             channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel) for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Value': [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for fresnel_inst in current_session.fresnel_analysis_instances]})

            exclusion_df = pd.DataFrame(
                {'Analysis': [current_session.exclusion_height_analysis_instances[exclusion_inst].name
                              for exclusion_inst in
                              current_session.exclusion_height_analysis_instances],
                 'Exclusion height mean [all]': ['{mean_} {all_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[
                         0], 2), all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :],
                     decimals=2))) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                 'Exclusion RI mean [all]': ['{mean_} {all_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0],
                     4), all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :],
                     decimals=4))) for exclusion_inst in current_session.exclusion_height_analysis_instances]})

            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=os.getcwd(), file_types=[('CSV files', '*.csv')], default_extension='.csv')

            fresnel_df.to_csv(save_filename[:-4] + '_fresnel' + '.csv', sep=';')
            exclusion_df.to_csv(save_filename[:-4] + '_exclusion' + '.csv', sep=';')

    app.run(debug=True, use_reloader=False, host=session_host, port=8050)
