# This is the main file where SPRpy is initiated. It is run by executing the file in a Python interpreter.
# The webapp is built using Dash (https://dash.plotly.com/), which is a Python framework for building webapps.
import math
import sys
import time
import tomllib
import types
import dash
import dash_bootstrap_components as dbc
import pandas as pd
import copy
import plotly
import plotly.express as px
import plotly.graph_objects as go

from SPRpy_classes import *
import SPRpy_multichannel
import SPRpy_series
from __about__ import version
from splash_screen import print_banner, show_splash, open_start_menu
import threading
import webbrowser
import os

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

# Font used in all figures (CSS font stack, the first available font is used)
PLOT_FONT = "Inter, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"


def style_figure(figure, **layout_kwargs):
    """
    Applies the common SPRpy look to a plotly figure.
    """
    layout = dict(font_family=PLOT_FONT, font_size=19, margin_r=25, margin_l=60, margin_t=40, template='simple_white')
    layout.update(layout_kwargs)
    figure.update_layout(**layout)
    figure.update_xaxes(mirror=True, showline=True)
    figure.update_yaxes(mirror=True, showline=True)
    return figure


COVERAGE_EQUATION_ACTIVE = {'opacity': 1, 'border': '2px solid #198754', 'border-radius': '8px', 'padding': '4px 14px'}
COVERAGE_EQUATION_INACTIVE = {'opacity': 0.45, 'border': '2px solid transparent', 'border-radius': '8px', 'padding': '4px 14px'}


def make_sensorgram_figure(sensorgram, offset_index, S_SPR, S_TIR, layer_thickness, decay_length, drift=None, coverage_points=None):
    """
    Sensorgram figure of the response quantification tab: SPR, TIR and bulk corrected angular shifts relative to the
    point offset_index. drift: optional dict with the SPR drift 'slope' (deg/min) that is removed (relative to the first
    time point). coverage_points: optional list of selected surface coverage points ({'index', 'curve'}), shown as markers.
    """
    time_values = np.asarray(sensorgram['time'], dtype=np.float64)
    traces = sensorgram_traces(sensorgram, S_SPR, S_TIR, layer_thickness, decay_length, drift=drift)
    offset_index = offset_index if 0 <= offset_index < len(time_values) else 0
    names = {'SPR angle': 'SPR angle' + (' (drift corrected)' if drift else ''), 'TIR angle': 'TIR angle', 'Bulk corrected': 'Bulk corrected'}
    colors = {'SPR angle': '#636efa', 'TIR angle': '#ef553b', 'Bulk corrected': '#00CC96'}
    figure = go.Figure()
    for trace in SENSORGRAM_TRACES:
        figure.add_trace(go.Scatter(x=time_values, y=traces[trace] - traces[trace][offset_index], name=names[trace], line_color=colors[trace]))
    if drift and drift.get('interval'):
        figure.add_vrect(x0=drift['interval'][0], x1=drift['interval'][1], fillcolor='#636efa', opacity=0.08, line_width=0, layer='below')
    if coverage_points:
        valid = [point for point in coverage_points if 0 <= point['index'] < len(time_values) and point['curve'] < len(SENSORGRAM_TRACES)]
        figure.add_trace(go.Scatter(
            x=[time_values[point['index']] for point in valid],
            y=[traces[SENSORGRAM_TRACES[point['curve']]][point['index']] - traces[SENSORGRAM_TRACES[point['curve']]][offset_index] for point in valid],
            mode='markers+text', text=[str(number + 1) for number in range(len(valid))], textposition='top center',
            marker=dict(symbol='diamond', size=13, color='black', line=dict(color='white', width=1)),
            name='Surface coverage points', hoverinfo='skip'))
    return style_figure(figure, xaxis_title=r'$\large{\text{Time [min]}}$', yaxis_title=r'$\large{\text{Angular shift [ }^{\circ}\text{ ]}}$',
                        uirevision=True, clickmode='event+select')


def make_TIR_fit_figure(fit_details):
    """
    Figure of the TIR angle determination of one spectrum (derivative, fit and TIR angle).
    """
    figure = go.Figure(go.Scatter(x=fit_details['TIR deriv x'], y=fit_details['TIR deriv y'], mode='lines+markers',
                                  marker_size=4, name='Derivative', line_color='#636efa'))
    figure.add_trace(go.Scatter(x=fit_details['TIR deriv fit x'], y=fit_details['TIR deriv fit y'], mode='lines',
                                name='Fit', line_color='#ef553b', line_width=3))
    if np.isfinite(fit_details['TIR angle']):
        figure.add_vline(x=fit_details['TIR angle'], line_dash='dash', line_color='black',
                         annotation_text='TIR {angle:.3f}°'.format(angle=fit_details['TIR angle']), annotation_position='top left')
    return style_figure(figure, xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                        yaxis_title=r'$\large{\text{TIR angular derivative}\text{}}$')


def fresnel_summary_barplot(session):
    """
    Bar plot of the fitted values of all fresnel analyses (result summary tab).
    """
    try:
        instances = session.fresnel_analysis_instances.values()
        x_barplot = [[fm.fitted_layer for fm in instances],
                     ['S' + str(fm.sensor_object.object_id) + ' ' + fm.fitted_layer for fm in instances]]
        y_barplot = [round(fm.fitted_result[0], 3) for fm in instances]
        figure = go.Figure(go.Bar(x=x_barplot, y=y_barplot))
    except Exception:
        figure = go.Figure(go.Bar(x=[0], y=[0]))
    figure.update_layout(yaxis_title='Fitted value', font_family=PLOT_FONT, font_size=19, margin_r=25, margin_l=60, margin_t=40,
                         template='simple_white', uirevision=True, height=600, width=900)
    figure.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
    figure.update_yaxes(mirror=True, showline=True)
    return figure


def exclusion_summary_barplot(session):
    """
    Bar plot of the exclusion height results (result summary tab): the mean height of each exclusion height analysis
    with the standard deviation over its injection steps, and the height of each step as points.
    """
    names, means, stds, step_x, step_y = [], [], [], [], []
    for eh in session.exclusion_height_analysis_instances.values():
        heights = np.asarray(eh.all_exclusion_results[0, :], dtype=np.float64) if np.ndim(eh.all_exclusion_results) == 2 else np.array([np.nan])
        name = 'EH{id} {name}'.format(id=eh.object_id, name=eh.name)
        names.append(name)
        finite = heights[np.isfinite(heights)]
        means.append(float(np.mean(finite)) if len(finite) else np.nan)
        stds.append(float(np.std(finite)) if len(finite) else np.nan)
        step_x += [name] * len(finite)
        step_y += [float(h) for h in finite]
    figure = go.Figure(go.Bar(x=names, y=means, error_y=dict(type='data', array=stds, visible=True), name='Mean ± std',
                              marker_color='#636efa'))
    figure.add_trace(go.Scatter(x=step_x, y=step_y, mode='markers', name='Injection steps', marker=dict(color='#ef553b', size=7)))
    figure.update_layout(yaxis_title='Exclusion height [nm]', font_family=PLOT_FONT, font_size=19, margin_r=25, margin_l=60, margin_t=40,
                         template='simple_white', uirevision=True, height=600, width=900, legend=dict(orientation='h', y=1.08))
    figure.update_xaxes(mirror=True, showline=True, autotickangles=[0, -90])
    figure.update_yaxes(mirror=True, showline=True)
    return figure


def fresnel_summary_table(session):
    """
    Children of the fresnel modelling table in the result summary tab.
    """
    instances = session.fresnel_analysis_instances
    header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
    body = [dash.html.Tbody([dash.html.Tr([
        dash.html.Td('FM' + str(fm.object_id) + ' ' + fm.name),
        dash.html.Td('S' + str(fm.sensor_object.object_id) + ' ' + fm.sensor_object.name),
        dash.html.Td('{layer}|{parameter}-{channel}'.format(layer=fm.fitted_layer, parameter=fm.sensor_object.optical_parameters.columns[fm.fitted_layer_index[1]], channel=fm.sensor_object.channel)),
        dash.html.Td(round(fm.fitted_result[0], 3))]) for fm in instances.values()])]
    return header + body


def exclusion_summary_table(session):
    """
    Children of the exclusion height table in the result summary tab.
    """
    header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean'), dash.html.Th('Exclusion height all'), dash.html.Th('Exclusion RI mean'), dash.html.Th('Exclusion RI all')]))]
    body = [dash.html.Tbody([dash.html.Tr([
        dash.html.Td(eh.name),
        dash.html.Td('{mean_} {all_}'.format(mean_=round(eh.mean_exclusion_height_result[0], 2), all_=str(np.round(eh.all_exclusion_results[0, :], decimals=2)))),
        dash.html.Td(str(np.round(eh.all_exclusion_results[0, :], decimals=2))),
        dash.html.Td(str(round(eh.mean_exclusion_RI_result[0], 4))),
        dash.html.Td(str(np.round(eh.all_exclusion_results[1, :], decimals=4)))]) for eh in session.exclusion_height_analysis_instances.values()])]
    return header + body


def channel_options(data_path, current_path, channel_spots):
    """
    Checklist options with the channels of a measurement, with their measurement spot (the channel of the current
    analysis is always included).
    """
    options, values = [], []
    for channel in SPRpy_multichannel.find_channel_files(data_path, channel_spots):
        is_current = channel['path'] == current_path
        label = channel['label'] + (' (spot {s})'.format(s=channel['spot']) if channel['spot'] is not None else '')
        options.append({'label': label + (', this analysis' if is_current else ''), 'value': channel['path'], 'disabled': is_current})
        values.append(channel['path'])
    return options, values


def by_spot(analyses, channel_spots):
    """
    Channel analyses sorted by measurement spot and channel number.
    """
    def key(analysis):
        spot = SPRpy_multichannel.channel_spot(analysis.initial_data_path, channel_spots)
        return (spot is None, spot or 0, SPRpy_multichannel.channel_number(analysis.initial_data_path) or 0)
    return sorted(analyses, key=key)


def spot_label(analysis, channel_spots):
    spot = SPRpy_multichannel.channel_spot(analysis.initial_data_path, channel_spots)
    return str(spot) if spot is not None else '-'


WAVELENGTH_NOTE = ('Channels of the same spot measure the same part of the sample and should agree, apart from the different '
                   'sensitivity of each wavelength (670 nm is generally most accurate for thin layers, longer wavelengths for thick '
                   'layers). Differences between spots can indicate an inhomogeneous sample. The spots are set with channel_spots in config.toml.')


def small_table(columns, rows):
    return dbc.Table([dash.html.Thead(dash.html.Tr([dash.html.Th(column) for column in columns])),
                      dash.html.Tbody([dash.html.Tr([dash.html.Td(value) for value in row]) for row in rows])],
                     bordered=True, size='sm', striped=True, style={'margin-top': '10px', 'font-size': '0.9rem'})


def fit_result_text(fresnel_analysis):
    """
    Fit result of a fresnel analysis, with the angle of the fitted model dip relative to the measured SPR angle.
    """
    if fresnel_analysis.fitted_result is None:
        return 'Fit result: None'
    text = 'Fit result: {res}'.format(res=round(fresnel_analysis.fitted_result[0], 4))
    mismatch = getattr(fresnel_analysis, 'dip_angle_mismatch', np.nan)
    if mismatch is not None and np.isfinite(mismatch):
        text += ' (model dip {m:+.1f} mdeg from the measured SPR angle{matched})'.format(
            m=mismatch * 1000, matched=', dip angle matched' if getattr(fresnel_analysis, 'match_dip_angle', False) else '')
    return text


def dip_log_text(fresnel_analysis):
    mismatch = getattr(fresnel_analysis, 'dip_angle_mismatch', np.nan)
    if mismatch is None or not np.isfinite(mismatch):
        return ''
    return ', model dip {m:+.2f} mdeg from the measured SPR angle{matched}'.format(
        m=mismatch * 1000, matched=' (matched)' if getattr(fresnel_analysis, 'match_dip_angle', False) else '')


def exclusion_height_result_texts(exclusion_height_analysis):
    """
    Result texts of an exclusion height analysis: mean height, mean RI, all heights and all RIs.
    """
    eh = exclusion_height_analysis
    if np.isnan(eh.mean_exclusion_height_result[0]):
        mean_height = 'Mean exclusion height: nan'
        mean_RI = 'Mean exclusion RI: nan'
    else:
        mean_height = 'Mean exclusion height: {mean} nm (std: {std})'.format(mean=round(eh.mean_exclusion_height_result[0], 2),
                                                                          std=round(eh.mean_exclusion_height_result[1], 2))
        sensitivity = np.asarray(getattr(eh, 'height_sensitivity', [np.nan]), dtype=np.float64)
        if np.any(np.isfinite(sensitivity)):
            mean_height += ', changes by {s} nm per 1 mdeg error in the TIR shift between buffer and probe'.format(s=round(float(np.nanmean(sensitivity)), 2))
        mean_RI = 'Mean exclusion RI: {mean} (std: {std})'.format(mean=round(eh.mean_exclusion_RI_result[0], 4),
                                                                 std=round(eh.mean_exclusion_RI_result[1], 4))

    if np.array_equal(eh.all_exclusion_results, np.array([[np.nan, np.nan], [np.nan, np.nan]]), equal_nan=True):
        all_heights = 'All exclusion heights: nan'
        all_RI = 'All exclusion RI: nan'
    else:
        all_heights = 'All exclusion heights: {h}'.format(h=np.round(eh.all_exclusion_results[0, :], decimals=2))
        TIR_shifts = getattr(eh, 'TIR_shifts', [])
        if len(TIR_shifts) > 0:
            all_heights += ' (TIR shifts buffer to probe [mdeg]: {shifts}, {method})'.format(
                shifts=np.round(np.asarray(TIR_shifts, dtype=np.float64) * 1000, 1),
                method='from TIR edge overlay' if getattr(eh, 'bulk_RI_method', 'TIR_angle') == 'edge_overlay' else 'from separate TIR angles')
        all_RI = 'All exclusion RI: {RI}'.format(RI=np.round(eh.all_exclusion_results[1, :], decimals=4))
    return mean_height, mean_RI, all_heights, all_RI


def make_SPR_fit_figure(fit_details):
    """
    Figure of the SPR angle determination of one spectrum (data around the dip, fit and SPR angle).
    """
    figure = go.Figure(go.Scatter(x=fit_details['SPR data x'], y=fit_details['SPR data y'], mode='lines+markers',
                                  marker_size=4, name='Reflectivity', line_color='#636efa'))
    figure.add_trace(go.Scatter(x=fit_details['SPR fit x'], y=fit_details['SPR fit y'], mode='lines',
                                name='Fit', line_color='#ef553b', line_width=3))
    if np.isfinite(fit_details['SPR angle']):
        figure.add_vline(x=fit_details['SPR angle'], line_dash='dash', line_color='black',
                         annotation_text='SPR {angle:.3f}°'.format(angle=fit_details['SPR angle']), annotation_position='top right')
    return style_figure(figure, xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                        yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$')


if __name__ == '__main__':

    # Read configuration parameters (missing entries are filled in from config.toml.defaults)
    config = load_config('config.toml', 'config.toml.defaults')
    show_start_menu = config.get('show_start_menu', True)
    show_splash_screen = config.get('show_splash_screen', True)
    open_browser = config.get('open_browser', True)
    print_banner(version)

    # Access individual parameters as variables
    ask_for_previous_session = config["ask_for_previous_session"]
    if config["default_data_folder"] == '':
        default_data_folder = os.path.expanduser('~')
    else:
        if os.path.exists(config["default_data_folder"]):
            default_data_folder = config["default_data_folder"]
        else:
            print('Warning! Custom data folder in config.toml does not exist.')
            default_data_folder = os.path.expanduser('~')
    if config['default_session_folder'] != '':
        default_session_folder = config['default_session_folder']
    else:
        default_session_folder = None
    session_host = config["session_host"]
    default_sensor_values = config["default_sensor_values"]
    channel_spots = config.get('channel_spots', [[1, 2], [3, 4]])
    default_dn_dC = config.get('dn_dC', 0.18)
    evanescent_decay_length = config["evanescent_decay_length"]
    instrument_SPR_sensitivity = config["instrument_SPR_sensitivity"]
    instrument_TIR_sensitivity = config["instrument_TIR_sensitivity"]
    TIR_default_parameters = config['TIR_fitting_parameters']
    SPR_default_parameters = config['SPR_fitting_parameters']

    # Exclusion height calculations no longer use parallel processes (max_logical_cores in config.toml is not used)
    logical_cores = None

    # Start menu: open a measurement, reopen a session or only convert .spr2 files (closing it quits SPRpy).
    # Without the start menu (show_start_menu = false in config.toml, or no display), the console prompts are used.
    start_choice = None
    loading_splash = None
    if show_start_menu:
        start_menu = open_start_menu(version, default_sensor_values, data_folder=default_data_folder,
                                     session_folder=default_session_folder or os.path.join(os.getcwd(), 'SPRpy sessions'))
        if start_menu is not None:
            start_choice = start_menu.run()
            if start_choice is None:
                print('SPRpy was closed from the start menu.')
                sys.exit(0)
            loading_splash = start_menu  # The menu window shows the loading status from here on

    load_session_flag = False
    session_file = None
    if start_choice is not None:
        if start_choice['action'] == 'session':
            load_session_flag = True
            session_file = start_choice['path']

    elif ask_for_previous_session is True:

        session_prompt = str(input(
            r'Would you like to load a previous session? Type "y" for yes, or simply skip by pressing enter.'))

        if session_prompt == 'y' or session_prompt == 'Y' or session_prompt == '"y"' or session_prompt == '\'y\'':
            session_file = select_file(r'Choose a previous session file', prompt_folder=default_session_folder)
            if not session_file:
                print('No session file selected. SPRpy was closed.')
                sys.exit(0)
            load_session_flag = True

    if load_session_flag:
            print('Loading previous session, please wait...')
            if loading_splash is not None:
                loading_splash.status('Loading session {name}...'.format(name=os.path.basename(os.path.dirname(session_file))))

            with open(session_file, 'rb') as file:
                current_session = pickle.load(file)

            # Assert that the session file is of the same version as the current SPRpy version
            if current_session.version != version and current_session.version.split('.')[0] == version.split('.')[0]:
                print('Note: the session file was created with SPRpy ' + current_session.version + ' (now ' + version + '). Its analyses keep '
                      'the settings they were made with (new options, e.g. matching the dip angle, can be switched on for them).')
            elif current_session.version != version and current_session.version.split('.')[0] == '1' and version.split('.')[0] == '2':
                print('Note: the session file was created with SPRpy ' + current_session.version + '. Sessions from SPRpy 1.x can be '
                      'opened with SPRpy 2: their analyses keep the methods and settings they were made with (e.g. the legacy TIR '
                      'method, or separate TIR angles for the exclusion height) until they are run again with new settings.')
            elif current_session.version != version:
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

            # Compatibility fixes for sessions created with older SPRpy versions
            for analysis_instance in current_session.fresnel_analysis_instances.values():
                if not hasattr(analysis_instance, 'fitted_layer_index'):
                    analysis_instance.fitted_layer_index = copy.deepcopy(analysis_instance.sensor_object.fitted_layer_index)
                if not hasattr(analysis_instance, 'fitted_layer'):
                    analysis_instance.fitted_layer = copy.deepcopy(analysis_instance.sensor_object.optical_parameters.iloc[analysis_instance.fitted_layer_index[0], 0])
                if getattr(analysis_instance, 'fitted_result', None) is None:
                    analysis_instance.fitted_result = np.array([np.nan, np.nan, np.nan])
                if not hasattr(analysis_instance, 'update_bulk_RI_from_TIR'):
                    analysis_instance.update_bulk_RI_from_TIR = True
            ensure_SPR_TIR_fitting_parameters(current_session.SPR_TIR_fitting_parameters, TIR_default_parameters)

            # Make sure the location and name of the session file is updated
            current_session.location = os.path.dirname(session_file)
            current_session.name = current_session.location.split('/')[-1]

            # Load measurement data
            try:
                current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(
                    path=current_session.current_data_path, default_data_folder=default_data_folder)
            except FileNotFoundError:
                print('The measurement file of the session was not found: ' + str(current_session.current_data_path))
                try:
                    current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(prompt='Select the original data file matching '+current_session.current_data_path)
                except FileSelectionCancelled:
                    print('No measurement file selected. SPRpy was closed.')
                    sys.exit(0)

            if loading_splash is None and show_splash_screen:
                loading_splash = show_splash(version, default_sensor_values)
            if loading_splash is not None:
                loading_splash.status('Calculating sensorgram of {file}...'.format(file=current_data_path.split('/')[-1]))
            sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, current_session.SPR_TIR_fitting_parameters)

            # Offset to start at 0 degrees at 0 minutes
            sensorgram_df_selection = copy.deepcopy(sensorgram_df)
            sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                                   sensorgram_df_selection['SPR angle'][0]
            sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                                   sensorgram_df_selection['TIR angle'][0]

            # Calculate bulk correction
            corrected_sensorgram_df_selection = sensorgram_df_selection['SPR angle'] - sensorgram_df_selection['TIR angle']*instrument_SPR_sensitivity[current_data_path[-9:-6]]/instrument_TIR_sensitivity*math.exp(-2*0/evanescent_decay_length[current_data_path[-9:-6]])

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
            current_session.add_log('Reopened session with SPRpy v{v}. Measurement: {file} ({n} scans)'.format(
                v=version, file=current_data_path.split('/')[-1], n=ydata_df.shape[0]))

    # If no previous session data was loaded
    if not load_session_flag:

        # Measurement data chosen in the start menu, or prompt the user for it
        print('Please wait...')
        try:
            current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(
                path=start_choice['path'] if start_choice is not None else False, default_data_folder=default_data_folder)
        except FileSelectionCancelled:
            print('No measurement file selected. SPRpy was closed.')
            sys.exit(0)

        # SPR and TIR fitting parameters from config.toml, and the TIR range (liquid or air) and scan speed dependent
        # settings from the measurement
        SPR_TIR_fitting_parameters = default_SPR_TIR_fitting_parameters(TIR_default_parameters, SPR_default_parameters)
        update_TIR_settings_for_measurement(SPR_TIR_fitting_parameters, scanspeed, angles_df, ydata_df, TIR_default_parameters)

        # Create initial session
        current_session = Session(version, SPR_TIR_fitting_parameters, directory=default_session_folder, current_data_path=current_data_path)

        if loading_splash is None and show_splash_screen:
            loading_splash = show_splash(version, default_sensor_values)
        if loading_splash is not None:
            loading_splash.status('Calculating sensorgram of {file}...'.format(file=current_data_path.split('/')[-1]))

        # Calculate sensorgram
        sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, current_session.SPR_TIR_fitting_parameters)

        # Offset to start at 0 degrees at 0 minutes
        sensorgram_df_selection = copy.deepcopy(sensorgram_df)
        sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                               sensorgram_df_selection['SPR angle'][0]
        sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                               sensorgram_df_selection['TIR angle'][0]

        # Calculate bulk correction
        corrected_sensorgram_df_selection = sensorgram_df_selection['SPR angle'] - sensorgram_df_selection[
            'TIR angle'] * instrument_SPR_sensitivity[current_data_path[-9:-6]] / instrument_TIR_sensitivity * math.exp(
            -2 * 0 / evanescent_decay_length[current_data_path[-9:-6]])

        # Add sensor object based on chosen measurement data
        current_sensor = add_sensor_backend(current_session, current_data_path, default_sensor_values)

        # Calculate TIR angle and update current_sensor.refractive_indices accordingly
        TIR_angle, _, _, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], current_session.SPR_TIR_fitting_parameters)
        current_sensor.refractive_indices[-1] = current_sensor.refractive_indices[0] * np.sin(np.pi / 180 * TIR_angle)
        current_sensor.optical_parameters.replace(current_sensor.optical_parameters['n'].iloc[-1], current_sensor.refractive_indices[-1], inplace=True)
        current_sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
            sensor_number=current_sensor.object_id,
            sensor_name=current_sensor.name,
            channel=current_sensor.channel,
            fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
            fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])

        current_session.add_log('New session with SPRpy v{v}. Measurement: {file} ({n} scans, scan speed {speed}, TIR range {TIR_range} deg). '
                                'TIR method: {TIR_method}, SPR method: {SPR_method}. Added sensor S{id} {name}'.format(
                                    v=version, file=current_data_path.split('/')[-1], n=ydata_df.shape[0], speed=scanspeed,
                                    TIR_range=current_session.SPR_TIR_fitting_parameters['TIR range'],
                                    TIR_method=current_session.SPR_TIR_fitting_parameters['TIR method'],
                                    SPR_method=current_session.SPR_TIR_fitting_parameters['SPR method'],
                                    id=current_sensor.object_id, name=current_sensor.name))
        current_session.save_session()
        current_session.save_sensor(current_sensor.object_id)

        # Add empty analysis objects
        current_fresnel_analysis = None
        current_exclusion_height_analysis = None

    # Log every exported figure or table in the session log
    _save_file_dialog = save_file

    def save_file(*args, **kwargs):
        saved_file = _save_file_dialog(*args, **kwargs)
        if not saved_file:
            raise dash.exceptions.PreventUpdate  # Save dialog cancelled
        current_session.add_log('Exported ' + saved_file, save=True)
        return saved_file

    # A cancelled file dialog when loading data in the GUI leaves everything as it was
    _load_csv_data = load_csv_data

    def load_csv_data(*args, **kwargs):
        try:
            return _load_csv_data(*args, **kwargs)
        except FileSelectionCancelled:
            raise dash.exceptions.PreventUpdate

    if loading_splash is not None:
        loading_splash.status('Building the user interface...')

    # Dash app
    app = dash.Dash(name='SPRpy', title='SPRpy', external_stylesheets=[dash_app_theme, 'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap'], compress=True)
    app._favicon = 'icon.ico'
    # app.config.suppress_callback_exceptions = True  # NOTE: Comment out this line for debugging purposes if callbacks do not fire when supposed to

    # Dash figures
    reflectivity_fig = px.line(reflectivity_df, x='angles', y='ydata')
    reflectivity_fig.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                   yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                   font_family=PLOT_FONT,
                                   font_size=19,
                                   margin_r=25,
                                   margin_l=60,
                                   margin_t=40,
                                   template='simple_white')
    reflectivity_fig.update_xaxes(mirror=True, showline=True)
    reflectivity_fig.update_yaxes(mirror=True, showline=True)

    sensorgram_fig = make_sensorgram_figure(sensorgram_df_selection, 0, instrument_SPR_sensitivity[current_data_path[-9:-6]],
                                            instrument_TIR_sensitivity, 0, evanescent_decay_length[current_data_path[-9:-6]])

    initial_fit_details = sensorgram_fit_details(angles_df, ydata_df.iloc[-1, :], current_session.SPR_TIR_fitting_parameters)
    TIR_fitting_fig = make_TIR_fit_figure(initial_fit_details)
    SPR_fitting_fig = make_SPR_fit_figure(initial_fit_details)

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
        font_family=PLOT_FONT,
        font_size=19,
        margin_r=25,
        margin_l=60,
        margin_t=40,
        template='simple_white',
        uirevision=True)
    d_n_pair_fig.update_xaxes(mirror=True, showline=True)
    d_n_pair_fig.update_yaxes(mirror=True, showline=True)

    result_barplot_fig = fresnel_summary_barplot(current_session)

    # Dash webapp layout
    app.layout = dash.html.Div([

        # Top bar
        dbc.Navbar(
            dbc.Container([
                dash.html.Div([
                    dash.html.Img(src='static/images/icon.svg', height='40px', className='me-2'),
                    dbc.NavbarBrand('SPRpy', className='sprpy-brand'),
                    dbc.Badge('v' + version, color='light', text_color='primary', pill=True, className='ms-2'),
                ], className='d-flex align-items-center'),
                dash.html.Div('Multi-parameter surface plasmon resonance analysis', className='sprpy-tagline d-none d-lg-block'),
                dbc.Button('Illustrations', id='illustrations-toggle', color='light', outline=True, size='sm', n_clicks=0,
                           title='Show or hide the illustrations of the methods'),
            ], fluid=True, className='d-flex justify-content-between align-items-center'),
            color='primary', dark=True, className='sprpy-navbar'),

        # Illustrations of the methods
        dbc.Collapse(
            dash.html.Div([
                dash.html.Figure([dash.html.Img(src='static/images/SPR_principle.svg'),
                                  dash.html.Figcaption('Surface plasmon resonance')]),
                dash.html.Figure([dash.html.Img(src='static/images/fresnel_material.svg'),
                                  dash.html.Figcaption('Fresnel modelling')]),
                dash.html.Figure([dash.html.Img(src='static/images/SPR_angular_spectrum.svg'),
                                  dash.html.Figcaption('Angular spectrum')]),
                dash.html.Figure([dash.html.Img(src='static/images/non-interacting_height_probe.PNG'),
                                  dash.html.Figcaption('Non-interacting height probing')]),
            ], className='sprpy-hero'),
            id='illustrations-collapse', is_open=True),

        dbc.Container([
        dbc.Row([
            # Session log
            dbc.Col(dbc.Card([
                dbc.CardHeader(dash.html.H5('{name_} - Session log'.format(name_=current_session.name),
                                            id='session-title', className='sprpy-card-title')),
                dbc.CardBody([
                    dash.dcc.Textarea(
                        id='console',
                        value=current_session.log,
                        readOnly=True,
                        className='sprpy-console'),
                    dash.dcc.Interval(id='session-log-refresh', interval=2000),
                    dash.dcc.Store(id='console-scroll-dummy'),
                    dbc.InputGroup(
                        [
                            dbc.Input(id='test-input', value='', type='text', placeholder='Write a note or a new session name...'),
                            dbc.Button('Add note to log', id='submit-button', n_clicks=0, color='info'),
                            dbc.Button('Rename session', id='rename-session-button', n_clicks=0, color='warning'),
                        ], size='sm', className='mt-2'),
                ]),
            ], className='sprpy-card h-100'), xl=5, className='mb-3'),

            # File and sensor controls
            dbc.Col(dbc.Card([
                dbc.CardHeader(dash.html.H5('Measurement and sensor', className='sprpy-card-title')),
                dbc.CardBody([
        dash.html.Div(['Current measurement file:    ', current_data_path.split('/')[-1]],
                      id='datapath-textfield',
                      className='sprpy-datapath'),
        dash.html.Div([
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
        ], className='sprpy-toolbar'),

        # Sensor datatable
            dash.html.Div([
                dash.html.H4(['S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
                    sensor_number=current_sensor.object_id,
                    sensor_name=current_sensor.name,
                    channel=current_sensor.channel,
                    fitted_layer=current_sensor.optical_parameters.iloc[current_sensor.fitted_layer_index[0], 0],
                    fitted_param=current_sensor.optical_parameters.columns[current_sensor.fitted_layer_index[1]])
                ], id='sensor-table-title', className='sprpy-sensor-title'),
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
                                              style_cell={'textAlign': 'center', 'fontFamily': 'inherit', 'padding': '6px'},
                                              style_table={'overflowX': 'auto'}),
                ]),
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
                ], size='sm', className='sprpy-table-buttons'),
            ], className='sprpy-sensor-table'),
                ]),
            ], className='sprpy-card h-100'), xl=7, className='mb-3'),
        ], className='g-3'),

            # Default sensor values
            dash.html.Div([
                dbc.Collapse(
                    dbc.Card(
                        dbc.CardBody(
                            dbc.Table.from_dataframe(pd.DataFrame(
                                default_sensor_values,
                            ), size='sm', striped=True, bordered=True, hover=True, responsive=True)
                        ), className='sprpy-card mb-3'),
                    id='default-values-collapse',
                    is_open=False)
            ], className='sprpy-default-values'),

        # Analysis tabs
        dbc.Card([
            dbc.CardBody([
            dbc.Tabs([

                # Response quantification tab
                dbc.Tab([
                    dash.html.Div([
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
                            dash.html.Div([
                                dash.dcc.Graph(id='quantification-sensorgram-graph',
                                               figure=sensorgram_fig,
                                               mathjax=True),
                                dbc.ButtonGroup([
                                    dbc.Switch(
                                        id='hover-selection-switch',
                                        label='Stop mouse hover updates',
                                        value=False),
                                    dbc.Switch(label='Show TIR/SPR fitting parameters',
                                               id='quantification-show-SPR-TIR-fit-options-switch',
                                               value=False),
                                    dbc.Switch(label=dash.html.Span('Select surface coverage points',
                                                                    title='When on, clicking the SPR or bulk corrected trace of the sensorgram selects the two '
                                                                          'points between which the surface coverage is calculated (see Surface coverage below), '
                                                                          'instead of setting the zero point of the sensorgram.',
                                                                    className='sprpy-help'),
                                               id='quantification-coverage-switch',
                                               value=False),
                                    dbc.DropdownMenu(
                                        id='sensorgram-save-dropdown',
                                        label='Save as...',
                                        color='info',
                                        children=[dbc.DropdownMenuItem('.PNG', id='quantification-sensorgram-save-png', n_clicks=0),
                                                  dbc.DropdownMenuItem('.SVG', id='quantification-sensorgram-save-svg', n_clicks=0),
                                                  dbc.DropdownMenuItem('.HTML', id='quantification-sensorgram-save-html', n_clicks=0),
                                                  dbc.DropdownMenuItem('.csv', id='quantification-sensorgram-save-csv', n_clicks=0)],
                                        style={'margin-left': '10%'}),
                                ], style={'margin-left': '10%'}),
                                dash.dcc.Store(id='sensorgram-drift-store', storage_type='memory', data=None),
                                dash.dcc.Store(id='quantification-offset-store', storage_type='memory', data=0),
                                dash.dcc.Store(id='quantification-coverage-points', storage_type='memory', data=[]),
                                dbc.InputGroup([
                                    dbc.InputGroupText(dash.html.Span('SPR drift correction, baseline from', title='Removes a linear drift of the SPR angle, fitted to a stable baseline (e.g. buffer before the injections). Drag a box over the baseline in the sensorgram (box select tool) or type the interval, then click "Correct drift". The TIR angle is not corrected.', className='sprpy-help')),
                                    dbc.Input(id='sensorgram-drift-start', type='number', placeholder='start', debounce=True, style={'max-width': '110px'}),
                                    dbc.InputGroupText('to'),
                                    dbc.Input(id='sensorgram-drift-end', type='number', placeholder='end', debounce=True, style={'max-width': '110px'}),
                                    dbc.InputGroupText('min'),
                                    dbc.Button('Correct drift', id='sensorgram-drift-apply', n_clicks=0, color='primary',
                                               title='Fit a straight line to the SPR angle within the baseline interval and remove its slope from the whole sensorgram.'),
                                    dbc.Button('Reset', id='sensorgram-drift-reset', n_clicks=0, color='secondary', outline=True,
                                               title='Remove the drift correction.'),
                                ], size='sm', className='sprpy-drift-group'),
                                dash.html.Div(id='sensorgram-drift-status', className='sprpy-hint'),
                            ], style={'width': '60%'}),
                        ], style={'display': 'flex', 'justify-content': 'center', 'margin-bottom': '20px'}),
                        dash.html.Div([  # TODO: Add layout options for changing TIR and SPR fit settings. Also add toggle for changing TIR angle algorithm to be more similar to Bionavis approach?
                            dbc.Collapse(
                                dbc.Card([
                                    dash.html.Div([
                                        dash.html.Div([
                                            dash.dcc.Graph(id='quantification-TIR-fit-graph',
                                                           figure=TIR_fitting_fig,
                                                           mathjax=True),
                                        ], style={'width': '45%'}),
                                        dash.html.Div([
                                            dash.dcc.Graph(id='quantification-SPR-fit-graph',
                                                           figure=SPR_fitting_fig,
                                                           mathjax=True),
                                        ], style={'width': '45%'}),
                                    ], style={'display': 'flex', 'justify-content': 'center'}),
                                    dbc.CardBody(
                                        dbc.Form([
                                            dash.html.H6('TIR angle', className='sprpy-form-heading'),
                                            dbc.Row([
                                                dbc.Col([
                                                    dbc.InputGroup([
                                                        dbc.InputGroupText('Method'),
                                                        dbc.Select(id='TIR-fit-option-method',
                                                                   options=[{'label': label, 'value': key} for key, label in TIR_METHODS.items()],
                                                                   value=current_session.SPR_TIR_fitting_parameters.get('TIR method', 'legacy')),
                                                    ], size='sm'),
                                                ], md=5),
                                                dbc.Col([
                                                    dbc.InputGroup([
                                                        dbc.InputGroupText('TIR angle range [°]'),
                                                        dbc.Input(id='TIR-fit-option-range-low',
                                                                  value=current_session.SPR_TIR_fitting_parameters['TIR range'][0],
                                                                  type='number', step=0.1),
                                                        dbc.Input(id='TIR-fit-option-range-high',
                                                                  value=current_session.SPR_TIR_fitting_parameters['TIR range'][1],
                                                                  type='number', step=0.1),
                                                    ], size='sm'),
                                                ], md=5),
                                            ], className='g-2 mb-2'),
                                            dash.html.Div([
                                                dbc.InputGroup([
                                                    dbc.InputGroupText(dash.html.Span('Derivative filter width [°]', title='Width of the Savitzky-Golay derivative filter. Wider = less noise, but distorts the TIR edge. Default 0.15', className='sprpy-help')),
                                                    dbc.Input(id='TIR-fit-option-smooth-width',
                                                              value=current_session.SPR_TIR_fitting_parameters.get('TIR smooth width deg', 0.15),
                                                              type='number', step=0.005, min=0.01),
                                                    dbc.InputGroupText(dash.html.Span('Peak fit half-width [°]', title='Half-width of the parabola fitted to the top of the derivative peak. Default 0.045', className='sprpy-help')),
                                                    dbc.Input(id='TIR-fit-option-fit-halfwidth',
                                                              value=current_session.SPR_TIR_fitting_parameters.get('TIR fit halfwidth deg', 0.045),
                                                              type='number', step=0.005, min=0.005),
                                                ], size='sm'),
                                            ], id='TIR-fit-sg-options', className='mb-2'),
                                            dash.html.Div([
                                                dbc.InputGroup([
                                                    dbc.InputGroupText('Smoothening window size'),
                                                    dbc.Input(id='TIR-fit-option-window',
                                                              value=current_session.SPR_TIR_fitting_parameters['TIR window count'],
                                                              type='number'),
                                                    dbc.InputGroupText('Nr of fit points'),
                                                    dbc.Input(id='TIR-fit-option-points',
                                                              value=current_session.SPR_TIR_fitting_parameters['TIR fit points'],
                                                              type='number'),
                                                    dbc.InputGroupText('Points below and above peak'),
                                                    dbc.Input(id='TIR-fit-option-below-peak',
                                                              value=current_session.SPR_TIR_fitting_parameters['points_below_TIR_peak'],
                                                              type='number'),
                                                    dbc.Input(id='TIR-fit-option-above-peak',
                                                              value=current_session.SPR_TIR_fitting_parameters['points_above_TIR_peak'],
                                                              type='number')
                                                ], size='sm'),
                                            ], id='TIR-fit-legacy-options', className='mb-2'),
                                            dash.html.H6('SPR angle', className='sprpy-form-heading'),
                                            dbc.Row([
                                                dbc.Col([
                                                    dbc.InputGroup([
                                                        dbc.InputGroupText('Method'),
                                                        dbc.Select(id='SPR-fit-option-method',
                                                                   options=[{'label': label, 'value': key} for key, label in SPR_METHODS.items()],
                                                                   value=current_session.SPR_TIR_fitting_parameters.get('SPR method', 'legacy')),
                                                    ], size='sm'),
                                                ], md=5),
                                                dbc.Col([
                                                    dash.html.Div([
                                                        dbc.InputGroup([
                                                            dbc.InputGroupText(dash.html.Span('Dip level', title='Fraction of the dip depth included in the fit window (0.1-0.5 all work well). The window adapts to the width of the dip. Default 0.2', className='sprpy-help')),
                                                            dbc.Input(id='SPR-fit-option-dip-level',
                                                                      value=current_session.SPR_TIR_fitting_parameters.get('SPR dip level', 0.2),
                                                                      type='number', step=0.05, min=0.02, max=0.9),
                                                        ], size='sm'),
                                                    ], id='SPR-fit-dip-options'),
                                                ], md=5),
                                            ], className='g-2 mb-2'),
                                            dash.html.Div([
                                                dbc.InputGroup([
                                                    dbc.InputGroupText('Nr of fit points'),
                                                    dbc.Input(id='SPR-fit-option-points',
                                                              value=current_session.SPR_TIR_fitting_parameters['SPR fit points'],
                                                              type='number'),
                                                    dbc.InputGroupText('Points below and above dip'),
                                                    dbc.Input(id='SPR-fit-option-below-peak',
                                                              value=current_session.SPR_TIR_fitting_parameters['sensorgram_angle_range_points'][0],
                                                              type='number'),
                                                    dbc.Input(id='SPR-fit-option-above-peak',
                                                              value=current_session.SPR_TIR_fitting_parameters['sensorgram_angle_range_points'][1],
                                                              type='number')
                                                ], size='sm'),
                                            ], id='SPR-fit-legacy-options', className='mb-2'),
                                            dbc.Button('Apply TIR/SPR fit options',
                                                       id='quantification-apply-fitting-SPR-TIR-button',
                                                       n_clicks=0,
                                                       color='success',
                                                       className='mt-2',
                                                       title='Recalculate the TIR and SPR sensorgrams with these settings')
                                        ], style={'margin-bottom': '10px'}),
                                    )
                                ]),
                                id='quantification-TIR-SPR-fit-collapse', is_open=False)
                        ], style={'margin': 'auto', 'width': '100%'}),
                        dash.html.Div([
                            dash.html.H4('Bulk correction parameters', style={'text-align': 'center'}),
                            dbc.Form([
                                dash.dcc.Markdown(
                                    '''
                                    $$\\Delta\\theta_{SPR}^{*}=\\Delta\\theta_{SPR}-\\Delta\\theta_{TIR}\\frac{S_{SPR}}{S_{TIR}}\\Large{e^{\\frac{-2d}{\\delta}}}$$
                                    ''',
                                    mathjax=True,
                                    style={'text-align': 'center'}),
                                dbc.Row([
                                    dbc.Col([
                                        dash.dcc.Markdown(
                                            '$d[nm]=$',
                                            mathjax=True)
                                    ], style={'padding-top': '7px', 'padding-left': '10px'}, width='auto'),
                                    dbc.Col([
                                        dbc.Input(id='sensorgram-correction-layer-thickness', value=0, type='number', style={'width': '110px'})
                                    ], width='auto'),
                                    dbc.Col([
                                        dash.dcc.Markdown(
                                            '$S_{SPR}[deg/RIU]=$',
                                            mathjax=True)
                                    ], style={'padding-top': '7px', 'padding-left': '10px'}, width='auto'),
                                    dbc.Col([
                                        dbc.Input(id='sensorgram-correction-layer-S_SPR',
                                                  value=instrument_SPR_sensitivity[current_data_path[-9:-6]],
                                                  type='number', style={'width': '110px'})
                                    ], width='auto'),
                                    dbc.Col([
                                        dash.dcc.Markdown(
                                            '$S_{TIR}[deg/RIU]=$',
                                            mathjax=True)
                                    ], style={'padding-top': '7px', 'padding-left': '10px'}, width='auto'),
                                    dbc.Col([
                                        dbc.Input(id='sensorgram-correction-layer-S_TIR',
                                                  value=instrument_TIR_sensitivity,
                                                  type='number', style={'width': '110px'})
                                    ], width='auto'),
                                    dbc.Col([
                                        dash.dcc.Markdown(
                                            r'$\delta[nm]=$',
                                            mathjax=True)
                                    ], style={'padding-top': '7px', 'padding-left': '10px'}, width='auto'),
                                    dbc.Col([
                                        dbc.Input(id='sensorgram-correction-layer-decay-length',
                                                  value=evanescent_decay_length[current_data_path[-9:-6]],
                                                  type='number', style={'width': '110px'})
                                    ], width='auto')
                                ], style={'textAlign': 'center'})
                            ])
                        ], style={'padding-top':'50px', 'margin': 'auto', 'width': '70%'}),
                        dash.html.Div([
                            dash.html.H4('Surface coverage', style={'text-align': 'center'}),
                            dash.html.Div([
                                dash.html.Div([
                                    dash.html.Div('Thin film approximation (d = 0)', id='coverage-equation-thin-label', className='sprpy-hint'),
                                    dash.dcc.Markdown(
                                        '''
                                        $$m_{S}=\\frac{\\delta\\Delta\\theta_{SPR}}{2\\frac{\\partial n}{\\partial C}S_{SPR}}$$
                                        ''',
                                        mathjax=True),
                                ], id='coverage-equation-thin', style=COVERAGE_EQUATION_ACTIVE),
                                dash.html.Div([
                                    dash.html.Div('Layer thickness d from the bulk correction parameters (d > 0)', id='coverage-equation-thick-label', className='sprpy-hint'),
                                    dash.dcc.Markdown(
                                        '''
                                        $$m_{S}=\\frac{d\\Delta\\theta_{SPR}}{\\frac{\\partial n}{\\partial C}S_{SPR}\\left(1-\\Large{e^{\\frac{-2d}{\\delta(\\lambda)}}}\\right)}$$
                                        ''',
                                        mathjax=True),
                                ], id='coverage-equation-thick', style=COVERAGE_EQUATION_INACTIVE),
                            ], style={'display': 'flex', 'justify-content': 'center', 'gap': '40px', 'text-align': 'center'}),
                            dbc.Row([
                                dbc.Col([
                                    dash.dcc.Markdown(r'$\frac{\partial n}{\partial C}[cm^3/g]=$', mathjax=True)
                                ], style={'padding-top': '2px', 'padding-left': '10px'}, width='auto'),
                                dbc.Col([
                                    dbc.Input(id='coverage-dn-dC', value=default_dn_dC, type='number', step=0.001, min=0, style={'width': '150px'})
                                ], width='auto'),
                                dbc.Col([
                                    dbc.Button('Clear points', id='coverage-clear-points', n_clicks=0, color='warning', outline=True,
                                               title='Clear the selected surface coverage points.')
                                ], width='auto'),
                            ], style={'textAlign': 'center', 'justify-content': 'center'}),
                            dash.html.Div('Turn on "Select surface coverage points" below the sensorgram and click two points on the SPR or bulk corrected trace '
                                          '(e.g. before an injection and after rinsing). The shift between them is used for the surface coverage, with '
                                          'δ and S_SPR from the bulk correction parameters above.',
                                          id='coverage-result', className='sprpy-hint', style={'text-align': 'center', 'margin-top': '10px'}),
                        ], style={'padding-top': '30px', 'padding-bottom': '30px', 'margin': 'auto', 'width': '70%'}),
                    ], id='quantification-tab-content')
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
                                        dbc.Button('Layer series',
                                                   id='series-button',
                                                   n_clicks=0,
                                                   color='primary',
                                                   title='Fit repeated samples measured after each added layer (e.g. clean SiO2, APS, SMCC, PEG), for all channels at once.'),
                                        dbc.Modal([
                                            dbc.ModalHeader(dbc.ModalTitle('Layer series')),
                                            dbc.ModalBody([
                                                dash.html.P(['Fits repeated samples that were measured after each added layer, for all channels (wavelengths) at once. '
                                                             'The sample and step of each measurement are read from the file names: ',
                                                             dash.html.Code('<date>-<time>_<sample>_<step>-L<channel>_<wavelength>nm.csv'),
                                                             ' (e.g. 20260116-135311_SiO2_1_APS-L3_670nm.csv: sample SiO2_1, step APS), and the steps are ordered by '
                                                             'measurement time. The first step is fitted with the base sensor below (its fitted variable, e.g. the SiO2 '
                                                             'thickness), and each following step adds one layer whose thickness is fitted with the earlier layers fixed at '
                                                             'the results of the same sample and channel.'], style={'font-size': '0.9rem'}),
                                                dbc.Button('Choose measurement files', id='series-choose-files', n_clicks=0, color='secondary'),
                                                dash.dcc.Store(id='series-files', storage_type='memory', data=None),
                                                dash.html.Div(id='series-detected', className='sprpy-hint', style={'margin-top': '10px'}),
                                                dbc.Row([
                                                    dbc.Col(dbc.Label('Base sensor (first step)'), width='auto'),
                                                    dbc.Col(dash.dcc.Dropdown(id='series-base-sensor', options=[], clearable=False), width=5),
                                                    dbc.Col(dash.html.Div(id='series-base-sensor-info', className='sprpy-hint'), width=True),
                                                ], className='g-2', style={'margin-top': '15px', 'align-items': 'center'}),
                                                dbc.Label('Layer added in each following step (layers listed in config.toml get n and k of each wavelength from there):',
                                                          style={'margin-top': '15px'}),
                                                dash.dash_table.DataTable(data=[],
                                                                          columns=[{'name': 'Step', 'id': 'Step', 'type': 'text', 'editable': False},
                                                                                   {'name': 'Layer', 'id': 'Layer', 'type': 'text'},
                                                                                   {'name': 'n', 'id': 'n', 'type': 'numeric'},
                                                                                   {'name': 'k', 'id': 'k', 'type': 'numeric'},
                                                                                   {'name': 'd guess [nm]', 'id': 'd guess [nm]', 'type': 'numeric'}],
                                                                          editable=True,
                                                                          id='series-recipe',
                                                                          style_header={'backgroundColor': '#446e9b', 'color': 'white', 'fontWeight': 'bold'},
                                                                          style_cell={'textAlign': 'center'}),
                                                dbc.Checklist(id='series-options',
                                                              options=[{'label': 'Fit reflectivity offset', 'value': 'offset'},
                                                                       {'label': 'Fit prism k', 'value': 'prism'},
                                                                       {'label': 'Bulk RI from TIR angle', 'value': 'bulk'},
                                                                       {'label': 'Match dip angle', 'value': 'dip'}],
                                                              value=['offset', 'prism', 'bulk', 'dip'], inline=True, style={'margin-top': '15px'}),
                                                dbc.ButtonGroup([
                                                    dbc.Button('Run layer series', id='series-run', n_clicks=0, color='success'),
                                                    dbc.Button('Export results (.csv)', id='series-export', n_clicks=0, color='info', disabled=True),
                                                ], style={'margin-top': '15px'}),
                                                dash.dcc.Store(id='series-results', storage_type='memory', data=None),
                                                dash.dcc.Store(id='series-export-signal', storage_type='memory', data=None),
                                                dbc.Spinner(dash.html.Div(id='series-result'), color='success', spinner_style={'margin-top': '10px'}),
                                            ]),
                                            dbc.ModalFooter(dbc.Button('Close', id='series-close', n_clicks=0, color='secondary')),
                                        ],
                                            id='series-modal',
                                            size='xl',
                                            is_open=False,
                                            backdrop='static',
                                            keyboard=False,
                                            scrollable=True),
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
                                                            dash.dcc.RangeSlider(value=[reflectivity_df['angles'].iloc[0], reflectivity_df['angles'].iloc[-1]],
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
                                                        dbc.Col([
                                                            dbc.Checkbox(
                                                                id='fresnel-analysis-bulk-RI-from-TIR',
                                                                label="Bulk RI from TIR angle?",
                                                                value=True),
                                                            dbc.Tooltip('When checked (default), the refractive index of the bulk layer is calculated from the TIR angle of the measurement '
                                                                        'every time the model is run. Uncheck it to use the bulk refractive index entered in the sensor table instead.',
                                                                        target='fresnel-analysis-bulk-RI-from-TIR', placement='top'),
                                                        ], width='auto', style={'padding-top': '20px'}),
                                                        dbc.Col([
                                                            dbc.Checkbox(
                                                                id='fresnel-analysis-match-dip',
                                                                label="Match dip angle?",
                                                                value=True),
                                                            dbc.Tooltip('When checked (default), the model dip is placed at the measured SPR angle (which carries the information '
                                                                        'about the layer), and the reflectivity offset and prism k describe the depth and width of the dip. '
                                                                        'Otherwise the reflectivity is fitted over the whole angle range, which can leave the model dip some '
                                                                        'millidegrees from the measured one, depending on the angle range (as in SPRpy <= 2.0.0).',
                                                                        target='fresnel-analysis-match-dip', placement='top'),
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
                    ], id='fresnel-tab-content', style={'display': 'flex', 'justify-content': 'center'}),
                            dash.html.Div([
                                dbc.Card(dbc.CardBody([
                                    dash.html.H5('All wavelengths'),
                                    dash.html.P('Repeats the current (fitted) fresnel model for the other channels of this measurement: same scan, layers and fit '
                                                'settings, with n and k of each wavelength from config.toml (layers not listed there keep their values). '
                                                'A sensor and a fresnel analysis are added for each channel, and reused when fitting again.',
                                                style={'font-size': '0.9rem'}),
                                    dbc.Checklist(id='fresnel-all-wavelengths-channels', options=[], value=[], inline=True),
                                    dash.html.Div('Select or add a fresnel analysis to list the channels of its measurement.', id='fresnel-all-wavelengths-channels-hint', className='sprpy-hint'),
                                    dbc.Button('Fit all wavelengths', id='fresnel-all-wavelengths-run', color='success', n_clicks=0, style={'margin-top': '10px'}),
                                    dbc.Spinner(dash.html.Div(id='fresnel-all-wavelengths-result'), color='success', spinner_style={'margin-top': '10px'})
                                ]))
                            ], style={'margin': '25px 5% 10px 5%'})
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
                                                            dbc.Label('Refractive index initial guess and fitting bounds (x0, min, max)', width='auto'),
                                                            dbc.Col([
                                                                dbc.InputGroup([
                                                                    dbc.Input(id='exclusion-height-option-RI-initialguess',
                                                                              value=float(1.35),
                                                                              type='number'),
                                                                    dbc.Input(id='exclusion-height-option-RI-lowerbound',
                                                                              value=float(1.2),
                                                                              type='number'),
                                                                    dbc.Input(id='exclusion-height-option-RI-upperbound',
                                                                              value=float(1.6),
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
                                                            dbc.Col([
                                                                dbc.Checkbox(
                                                                    id='exclusion-height-match-dip',
                                                                    label=dash.html.Span('Match dip angle?', className='sprpy-help',
                                                                                         title='When checked (default), the model dip of each buffer and probe spectrum is placed at its '
                                                                                               'measured SPR angle when fitting the RI at each height (as in SPR_WetHeight.m), instead of '
                                                                                               'fitting the reflectivity over the whole angle range.'),
                                                                    value=True)
                                                            ], width='auto'),
                                                        ]),
                                                        dbc.Row([
                                                            dbc.Label(dash.html.Span('Bulk RI change from',
                                                                                     title='The exclusion height depends strongly on the bulk RI change between buffer and probe '
                                                                                           '(about 1 nm per mdeg of TIR shift). Overlaying the TIR edges of buffer and probe gives a '
                                                                                           'shift that does not depend on the TIR angle settings. Applied when initializing the model.',
                                                                                     className='sprpy-help'), width='auto'),
                                                            dbc.Col([
                                                                dbc.Select(id='exclusion-height-bulk-RI-method',
                                                                           options=[{'label': label, 'value': key} for key, label in EXCLUSION_BULK_RI_METHODS.items()],
                                                                           value='edge_overlay')
                                                            ], width=7)
                                                        ], style={'margin-top': '10px', 'margin-bottom': '10px'}),
                                                        dbc.Row([
                                                            dbc.Label(dash.html.Span('Height search',
                                                                                     title='Adaptive: a coarse sweep (31 heights between the height bounds) finds where the buffer and '
                                                                                           'probe curves cross, then the curves are calculated with "Resolution" points in a narrow '
                                                                                           'window around the crossing, which is located within 0.001 nm with a root finder. A wide '
                                                                                           'height range (e.g. up to 300 nm for thick brushes) therefore costs little time, and the '
                                                                                           'curves are well resolved around the crossing (30-100 points recommended). Fixed grid: all '
                                                                                           '"Resolution" height steps between the bounds are calculated and the crossing is '
                                                                                           'interpolated between them (as in SPRpy <= 1.2).',
                                                                                     className='sprpy-help'), width='auto'),
                                                            dbc.Col([
                                                                dbc.Select(id='exclusion-height-search-method',
                                                                           options=[{'label': label, 'value': key} for key, label in EXCLUSION_HEIGHT_SEARCH_METHODS.items()],
                                                                           value='adaptive')
                                                            ], width=7)
                                                        ], style={'margin-bottom': '10px'}),
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

                        dash.html.Div([
                            dbc.Card(dbc.CardBody([
                                dash.html.H5('All wavelengths'),
                                dash.html.P('Runs this analysis (same points and settings as above) for each selected channel, with the fresnel '
                                            'background of that channel. First fit the fresnel background for all wavelengths in the fresnel '
                                            'modelling tab.',
                                            style={'font-size': '0.9rem'}),
                                dbc.Checklist(id='exclusion-all-wavelengths-channels', options=[], value=[], inline=True),
                                dash.html.Div('Select or add an exclusion height analysis to list the channels of its measurement.', id='exclusion-all-wavelengths-channels-hint', className='sprpy-hint'),
                                dbc.Button('Run all wavelengths', id='exclusion-all-wavelengths-run', color='success', n_clicks=0, style={'margin-top': '10px'}),
                                dbc.Spinner(dash.html.Div(id='exclusion-all-wavelengths-result'), color='success', spinner_style={'margin-top': '10px'})
                            ]))
                        ], style={'margin': '25px 5% 10px 5%'}),
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
                                    # table_header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
                                    #             table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td('FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td('S'+ str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name), dash.html.Td('{layer}|{parameter}-{channel}'.format(
                                    #
                                    dbc.Table.from_dataframe(pd.DataFrame(
                                        {'Analysis': ['FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Sensor': ['S'+ str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name for fresnel_inst in current_session.fresnel_analysis_instances],
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
                                         'Exclusion height mean': ['{mean_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[0], 2)) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion height all': ['{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2))) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion RI mean': ['{mean_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0], 4)) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion RI all': ['{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :], decimals=4))) for exclusion_inst in current_session.exclusion_height_analysis_instances]
                                         }),
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
                                        ),
                                        dash.dcc.Store(id='summary-plots-save-signal', storage_type='memory', data=None),
                                        dash.dcc.Graph(id='summary-exclusion-barplot-graph',
                                                       figure=exclusion_summary_barplot(current_session),
                                                       mathjax=True, style={'margin-top': '30px'}),
                                        dbc.DropdownMenu(
                                            id='exclusion-barplot-save-dropdown',
                                            label='Save as...',
                                            color='info',
                                            children=[
                                                dbc.DropdownMenuItem('.PNG', id='exclusion-barplot-save-png', n_clicks=0),
                                                dbc.DropdownMenuItem('.SVG', id='exclusion-barplot-save-svg', n_clicks=0),
                                                dbc.DropdownMenuItem('.HTML', id='exclusion-barplot-save-html', n_clicks=0),
                                                dbc.DropdownMenuItem('.csv', id='exclusion-barplot-save-csv', n_clicks=0)],
                                            style={'margin-left': '299px'}
                                        )
                                    ]),
                                ]),
                            ], style={'margin-top': '20px', 'display': 'flex', 'justify-content': 'center'}),
                        ], style={'margin-top': '1.9rem', 'margin-left': '5%'}),
                    ], id='summary-tab-content')
                ], label='Result summary and export', tab_id='summary-tab', style={'margin-top': '10px'}),
            ], id='analysis-tabs', active_tab='quantification-tab', className='sprpy-tabs'),
            ]),
        ], className='sprpy-card mb-4'),
        ], fluid=True, className='sprpy-main'),
    ], className='sprpy-app')

    # Show or hide the illustrations at the top of the page
    @dash.callback(
        dash.Output('illustrations-collapse', 'is_open'),
        dash.Input('illustrations-toggle', 'n_clicks'),
        dash.State('illustrations-collapse', 'is_open'),
        prevent_initial_call=True)
    def toggle_illustrations(n_clicks, is_open):
        return not is_open

    # Show new automatic session log messages
    @dash.callback(
        dash.Output('console', 'value', allow_duplicate=True),
        dash.Input('session-log-refresh', 'n_intervals'),
        dash.State('console', 'value'),
        prevent_initial_call=True)
    def refresh_session_log(n_intervals, shown_log):
        if shown_log == current_session.log:
            raise dash.exceptions.PreventUpdate
        return current_session.log

    # Keep the session log scrolled to the latest message
    dash.clientside_callback(
        """
        function(value) {
            const console_element = document.getElementById('console');
            if (console_element) {
                setTimeout(function() { console_element.scrollTop = console_element.scrollHeight; }, 50);
            }
            return window.dash_clientside.no_update;
        }
        """,
        dash.Output('console-scroll-dummy', 'data'),
        dash.Input('console', 'value'))

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

            old_name = current_session.name
            current_session.update_name_and_location(state1)
            current_session.add_log('Renamed session "{old}" to "{new}"'.format(old=old_name, new=state1))
            current_session.save_session()
            new_name = '{name_} - Session log'.format(name_=current_session.name)

            return current_session.log, new_name, ''


    # TODO: Include logic for updating fitting parameters for TIR and SPR angle when calculating sensorgram. Also to select TIR fitting algorithm (implement something similar to Bionavis)
    # Load in new measurement data and send a Store signal to other callbacks to update appropriately
    @dash.callback(
        dash.Output('loaded-new-measurement', 'data', allow_duplicate=True),
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
        global corrected_sensorgram_df_selection
        global TIR_default_parameters
        global default_data_folder

        if 'load-data' == dash.ctx.triggered_id:
            # Get the folder location of the last loaded datafile
            previous_path = os.path.dirname(current_data_path)

            # Load measurement data and update session current data path
            current_data_path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(default_data_folder=previous_path)
            current_session.current_data_path = current_data_path

            # Choose the TIR range (liquid or air) from the TIR edge in the data, and the scan speed dependent TIR settings
            update_TIR_settings_for_measurement(current_session.SPR_TIR_fitting_parameters, scanspeed, angles_df, ydata_df, TIR_default_parameters)
            current_session.add_log('Loaded measurement {file} ({n} scans, scan speed {speed}, TIR range {TIR_range} deg)'.format(
                file=current_data_path, n=ydata_df.shape[0], speed=scanspeed, TIR_range=current_session.SPR_TIR_fitting_parameters['TIR range']))

            current_session.save_session()

            sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, current_session.SPR_TIR_fitting_parameters)

            # Offset to start at 0 degrees at 0 minutes
            sensorgram_df_selection = copy.deepcopy(sensorgram_df)
            sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                                   sensorgram_df_selection['SPR angle'][0]
            sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                                   sensorgram_df_selection['TIR angle'][0]

            # Calculate bulk correction
            corrected_sensorgram_df_selection = sensorgram_df_selection['SPR angle'] - sensorgram_df_selection[
                'TIR angle'] * instrument_SPR_sensitivity[current_data_path[-9:-6]] / instrument_TIR_sensitivity * math.exp(-2 * 0 / evanescent_decay_length[current_data_path[-9:-6]])

            return 'signal', ['Current measurement file:    ', current_data_path.split('/')[-1]], dash.no_update

        elif 'batch-fresnel-analysis-choose-files' == dash.ctx.triggered_id:
            print('Select the measurement data files (.csv)')
            batch_data_paths_ = select_files('Select the measurement data files', prompt_folder=os.path.dirname(current_data_path), file_types=[('CSV files', '*.csv')])
            if not batch_data_paths_:
                raise dash.exceptions.PreventUpdate  # File dialog cancelled
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

        if 'new-sensor-gold' == dash.ctx.triggered_id:
            current_sensor = add_sensor_backend(current_session, current_data_path, default_sensor_values, sensor_metal='Au')
            current_sensor.name = 'Gold sensor'
            current_session.add_log('Added new gold sensor S{id} ({channel})'.format(id=current_sensor.object_id, channel=current_sensor.channel))

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], current_session.SPR_TIR_fitting_parameters)
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
            current_sensor = add_sensor_backend(current_session, current_data_path, default_sensor_values, sensor_metal='SiO2')
            current_sensor.name = 'Glass sensor'
            current_session.add_log('Added new glass sensor S{id} ({channel})'.format(id=current_sensor.object_id, channel=current_sensor.channel))

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], current_session.SPR_TIR_fitting_parameters)
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
            current_sensor = add_sensor_backend(current_session, current_data_path, default_sensor_values, sensor_metal='Pd')
            current_sensor.name = 'Palladium sensor'
            current_session.add_log('Added new palladium sensor S{id} ({channel})'.format(id=current_sensor.object_id, channel=current_sensor.channel))

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], current_session.SPR_TIR_fitting_parameters)
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
            current_sensor = add_sensor_backend(current_session, current_data_path, default_sensor_values, sensor_metal='Pt')
            current_sensor.name = 'Platinum sensor'
            current_session.add_log('Added new platinum sensor S{id} ({channel})'.format(id=current_sensor.object_id, channel=current_sensor.channel))

            # Calculate TIR angle and bulk refractive index from measured data
            TIR_angle, _, _, _, _ = TIR_determination(reflectivity_df['angles'], reflectivity_df['ydata'], current_session.SPR_TIR_fitting_parameters)
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
            current_session.add_log('Copied sensor S{old} {name} to S{new}'.format(old=current_sensor.object_id, name=current_sensor.name, new=new_sensor.object_id))
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
            old_path = current_session.location + '/Sensors/S{id} {name}.pickle'.format(id=current_sensor.object_id,
                                                                                         name=current_sensor.name)
            os.remove(old_path)

            # Change sensor name and save new sensor pickle file and session
            current_session.add_log('Renamed sensor S{id} "{old}" to "{new}"'.format(id=current_sensor.object_id, old=current_sensor.name, new=sensor_name_))
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
                current_session.add_log('Removed sensor S{id} {name}'.format(id=removed.object_id, name=removed.name))
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
            current_session.add_log('Saved edited values of sensor S{id} {name}: {layers}'.format(
                id=current_sensor.object_id, name=current_sensor.name,
                layers='; '.join('{l} d={d} n={n} k={k}'.format(l=row['Layers'], d=row['d [nm]'], n=row['n'], k=row['k']) for row in table_rows)))

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
            current_session.add_log('Fitted variable of sensor S{id} set to {layer}|{param}'.format(
                id=current_sensor.object_id, layer=current_sensor.optical_parameters.iloc[active_cell['row'], 0],
                param=current_sensor.optical_parameters.columns[active_cell['column']]))

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
        global sensorgram_df
        global current_sensor
        global reflectivity_df
        global current_data_path

        figure_object = go.Figure(figure_JSON)

        # Update based on hover over sensorgram figure
        if 'quantification-sensorgram-graph' == dash.ctx.triggered_id:

            # # First make sure no other traces has been added and the very first value is ignored
            # if figure_object.data.__len__() == 1:

            # Then also make sure lock hover switch is set to inactive
            if lock_hover is False:
                time_index = hoverData['points'][0]['pointIndex']
                SPR_angle = float(sensorgram_df['SPR angle'][time_index])
                TIR_angle = float(sensorgram_df['TIR angle'][time_index])
                reflectivity_df['ydata'] = ydata_df.loc[time_index+1]

                new_figure = go.Figure(go.Scatter(x=reflectivity_df['angles'],
                                                  y=reflectivity_df['ydata'],
                                                  mode='lines',
                                                  showlegend=False,
                                                  line_color='#636efa'))
                new_figure.add_trace(go.Scatter(x=[SPR_angle, SPR_angle],
                                                y=[reflectivity_df['ydata'].min(axis=0)-0.05, reflectivity_df['ydata'].min(axis=0)+0.05],
                                                mode='lines',
                                                showlegend=False,
                                                line_color='#636efa'
                                                ))
                new_figure.add_trace(go.Scatter(x=[TIR_angle, TIR_angle],
                                                y=[reflectivity_df['ydata'].max(axis=0)-0.07, reflectivity_df['ydata'].max(axis=0)],
                                                mode='lines',
                                                showlegend=False,
                                                line_color='#ef553b',
                                                ))
                new_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                         yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                         font_family=PLOT_FONT,
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
            # else:
            #     return dash.no_update

        # This adds a trace to the reflectivity plot from a separate measurement file. The trace data is not stored.
        elif 'quantification-reflectivity-add-data-trace' == dash.ctx.triggered_id:
            _, _, _, _, _, trace_reflectivity_df = load_csv_data(default_data_folder=os.path.dirname(current_data_path))
            figure_object.add_trace(go.Scatter(x=trace_reflectivity_df['angles'],
                                               y=trace_reflectivity_df['ydata'],
                                               mode='lines',
                                               showlegend=True))
            return figure_object

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
                                               showlegend=True))
            return figure_object

        # Clear added traces from the reflectivity plot
        elif 'quantification-reflectivity-clear-traces' == dash.ctx.triggered_id:
            new_figure = go.Figure(go.Scatter(x=reflectivity_df['angles'],
                                              y=reflectivity_df['ydata'],
                                              mode='lines',
                                              showlegend=False,
                                              line_color='#636efa'))
            new_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                     yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                     font_family=PLOT_FONT,
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
                                     font_family=PLOT_FONT,
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
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
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
    @dash.callback(
        dash.Output('quantification-sensorgram-graph', 'figure'),
        dash.Output('quantification-offset-store', 'data'),
        dash.Input('quantification-sensorgram-save-png', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-svg', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-html', 'n_clicks'),
        dash.Input('quantification-sensorgram-save-csv', 'n_clicks'),
        dash.Input('quantification-sensorgram-graph', 'clickData'),
        dash.Input('loaded-new-measurement', 'data'),
        dash.Input('sensorgram-correction-layer-thickness', 'value'),
        dash.Input('sensorgram-correction-layer-S_SPR', 'value'),
        dash.Input('sensorgram-correction-layer-S_TIR', 'value'),
        dash.Input('sensorgram-correction-layer-decay-length', 'value'),
        dash.Input('sensorgram-drift-store', 'data'),
        dash.Input('quantification-coverage-points', 'data'),
        dash.State('sensorgram-correction-layer-thickness', 'value'),
        dash.State('sensorgram-correction-layer-S_SPR', 'value'),
        dash.State('sensorgram-correction-layer-S_TIR', 'value'),
        dash.State('sensorgram-correction-layer-decay-length', 'value'),
        dash.State('quantification-sensorgram-graph', 'figure'),
        dash.State('quantification-coverage-switch', 'value'),
        dash.State('quantification-offset-store', 'data'),
        prevent_initial_call=True)  # Adding this fixed a weird bug with graph not updating after firing clickData callbacks
    def update_sensorgram_quantification_tab(save_png, save_svg, save_html, save_csv, clickData, data_update, layer_thickness, S_SPR, S_TIR, decay_length, drift, coverage_points, layer_thickness_state, S_SPR_state, S_TIR_state, decay_length_state, figure_JSON, coverage_mode, offset_state):

        figure_object = go.Figure(figure_JSON)
        global sensorgram_df_selection
        global current_data_path

        if 'loaded-new-measurement' == dash.ctx.triggered_id:
            return make_sensorgram_figure(sensorgram_df_selection, 0, S_SPR_state, S_TIR_state, layer_thickness_state, decay_length_state, drift=None), 0

        elif 'quantification-sensorgram-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'quantification-sensorgram-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
            fig_df = pd.DataFrame(data={'Time': list(figure_object.data[0].x["_inputArray"].values())[:-3], 'SPR': list(figure_object.data[0].y["_inputArray"].values())[:-3], 'TIR': list(figure_object.data[1].y["_inputArray"].values())[:-3], 'Bulk corrected': list(figure_object.data[2].y["_inputArray"].values())[:-3]})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        else:
            offset_index = offset_state or 0
            # A click sets the zero point, unless surface coverage points are being selected (the figure is then redrawn
            # with the selected points, which arrive in the same update) or a coverage point marker was clicked
            if 'quantification-sensorgram-graph.clickData' in dash.ctx.triggered_prop_ids and not coverage_mode \
                    and clickData['points'][0].get('curveNumber', 0) < len(SENSORGRAM_TRACES):
                offset_index = clickData['points'][0]['pointIndex']
            if offset_index >= len(sensorgram_df_selection):
                offset_index = 0

            return make_sensorgram_figure(sensorgram_df_selection, offset_index, S_SPR_state, S_TIR_state, layer_thickness_state, decay_length_state,
                                          drift=drift, coverage_points=coverage_points), offset_index

    # SPR drift correction in the response quantification tab
    @dash.callback(
        dash.Output('sensorgram-drift-store', 'data'),
        dash.Output('sensorgram-drift-start', 'value'),
        dash.Output('sensorgram-drift-end', 'value'),
        dash.Output('sensorgram-drift-status', 'children'),
        dash.Input('sensorgram-drift-apply', 'n_clicks'),
        dash.Input('sensorgram-drift-reset', 'n_clicks'),
        dash.Input('quantification-sensorgram-graph', 'selectedData'),
        dash.Input('loaded-new-measurement', 'data'),
        dash.State('sensorgram-drift-start', 'value'),
        dash.State('sensorgram-drift-end', 'value'),
        prevent_initial_call=True)
    def update_sensorgram_drift_correction(apply_button, reset_button, selectedData, new_measurement, drift_start, drift_end):

        if 'quantification-sensorgram-graph' == dash.ctx.triggered_id:
            # Fill in the baseline interval from a box selection in the sensorgram
            if not selectedData:
                raise dash.exceptions.PreventUpdate
            if selectedData.get('range') and selectedData['range'].get('x'):
                selected_times = selectedData['range']['x']
            elif selectedData.get('points'):
                selected_times = [point['x'] for point in selectedData['points']]
            else:
                raise dash.exceptions.PreventUpdate
            return dash.no_update, round(min(selected_times), 2), round(max(selected_times), 2), 'Baseline interval selected, click "Correct drift" to apply.'

        elif 'sensorgram-drift-apply' == dash.ctx.triggered_id:
            if drift_start is None or drift_end is None:
                return dash.no_update, dash.no_update, dash.no_update, 'Select a baseline interval first (box select in the sensorgram, or type the start and end time).'
            SPR_angles = sensorgram_df['SPR angle'].replace(0, np.nan)  # Spectra where no SPR angle was found are 0
            slope = linear_drift(sensorgram_df['time'], SPR_angles, [drift_start, drift_end])
            if slope is None:
                return dash.no_update, dash.no_update, dash.no_update, 'Fewer than two sensorgram points in the baseline interval, drift correction not applied.'
            current_session.add_log(
                'SPR drift correction of {slope:.3g} mdeg/min (baseline {start}-{end} min) applied to {file}'.format(
                    slope=slope * 1000, start=min(drift_start, drift_end), end=max(drift_start, drift_end), file=current_data_path.split('/')[-1]))
            return ({'slope': slope, 'interval': [min(drift_start, drift_end), max(drift_start, drift_end)]}, dash.no_update, dash.no_update,
                    'SPR drift of {slope:.3g} mdeg/min removed (fitted between {start} and {end} min).'.format(
                        slope=slope * 1000, start=min(drift_start, drift_end), end=max(drift_start, drift_end)))

        # Reset, or a new measurement was loaded
        return None, dash.no_update, dash.no_update, ''

    # Surface coverage between two selected sensorgram points
    @dash.callback(
        dash.Output('quantification-coverage-points', 'data'),
        dash.Output('coverage-result', 'children'),
        dash.Output('coverage-equation-thin', 'style'),
        dash.Output('coverage-equation-thick', 'style'),
        dash.Output('coverage-equation-thin-label', 'children'),
        dash.Output('coverage-equation-thick-label', 'children'),
        dash.Input('quantification-sensorgram-graph', 'clickData'),
        dash.Input('coverage-clear-points', 'n_clicks'),
        dash.Input('loaded-new-measurement', 'data'),
        dash.Input('coverage-dn-dC', 'value'),
        dash.Input('sensorgram-correction-layer-thickness', 'value'),
        dash.Input('sensorgram-correction-layer-S_SPR', 'value'),
        dash.Input('sensorgram-correction-layer-S_TIR', 'value'),
        dash.Input('sensorgram-correction-layer-decay-length', 'value'),
        dash.Input('sensorgram-drift-store', 'data'),
        dash.State('quantification-coverage-switch', 'value'),
        dash.State('quantification-coverage-points', 'data'),
        prevent_initial_call=True)
    def update_surface_coverage(clickData, clear_clicks, new_measurement, dn_dC, layer_thickness, S_SPR, S_TIR, decay_length, drift, coverage_mode, points):
        """
        Selects two points on the SPR or bulk corrected trace of the sensorgram and calculates the surface coverage
        m_S = delta * delta_theta_SPR / (2 * dn/dC * S_SPR) from the shift between them.
        """
        points = list(points or [])
        trigger = dash.ctx.triggered_id

        # Thin film approximation for d = 0, otherwise the layer thickness d of the bulk correction parameters
        try:
            thickness = float(layer_thickness or 0)
        except (TypeError, ValueError):
            thickness = 0.0
        thin_film = thickness <= 0
        equation_outputs = (COVERAGE_EQUATION_ACTIVE if thin_film else COVERAGE_EQUATION_INACTIVE,
                            COVERAGE_EQUATION_INACTIVE if thin_film else COVERAGE_EQUATION_ACTIVE,
                            'Thin film approximation (d = 0)' + (' - used' if thin_film else ''),
                            'Layer thickness d = {d:g} nm (from the bulk correction parameters)'.format(d=thickness) + (' - used' if not thin_film else '')
                            if not thin_film else 'Layer thickness d from the bulk correction parameters (d > 0)')

        if trigger in ('coverage-clear-points', 'loaded-new-measurement'):
            return ([], 'Points cleared. Turn on "Select surface coverage points" and click two points on the SPR or bulk corrected trace.') + equation_outputs

        if trigger == 'quantification-sensorgram-graph':
            if not coverage_mode or not clickData:
                raise dash.exceptions.PreventUpdate
            clicked = clickData['points'][0]
            curve = clicked.get('curveNumber', 0)
            if curve >= len(SENSORGRAM_TRACES):
                raise dash.exceptions.PreventUpdate  # Click on a selected point marker
            if SENSORGRAM_TRACES[curve] == 'TIR angle':
                return (dash.no_update, 'Select the points on the SPR or bulk corrected trace (not the TIR trace).') + equation_outputs
            point = {'index': int(clicked['pointIndex']), 'curve': int(curve)}
            points = [point] if len(points) >= 2 else points + [point]
            if len(points) == 1:
                return (points, 'Point 1 selected at {t:.2f} min, now select point 2.'.format(t=float(sensorgram_df_selection['time'].iloc[point['index']]))) + equation_outputs

        if len(points) < 2:
            return (dash.no_update, dash.no_update) + equation_outputs

        # Shift between the two points on the trace of the first point (recalculated with the current settings)
        traces = sensorgram_traces(sensorgram_df_selection, S_SPR, S_TIR, layer_thickness, decay_length, drift=drift)
        trace_name = SENSORGRAM_TRACES[points[0]['curve']]
        trace = traces[trace_name]
        if max(points[0]['index'], points[1]['index']) >= len(trace):
            return ([], 'The selected points do not belong to the current measurement, select them again.') + equation_outputs
        times = [float(sensorgram_df_selection['time'].iloc[point['index']]) for point in points]
        delta_theta = float(trace[points[1]['index']] - trace[points[0]['index']])
        note = '' if points[0]['curve'] == points[1]['curve'] else ' (point 2 was on another trace, the {name} trace of point 1 is used)'.format(name=trace_name)
        try:
            coverage = surface_coverage(delta_theta, float(decay_length), float(dn_dC), float(S_SPR), thickness=thickness)
        except (TypeError, ValueError, ZeroDivisionError):
            return (points, 'Δθ = {d:.2f} mdeg ({name}, {t1:.2f} to {t2:.2f} min). Enter dn/dC, δ and S_SPR to calculate the surface coverage.'.format(
                d=delta_theta * 1000, name=trace_name, t1=times[0], t2=times[1])) + equation_outputs

        result = 'Δθ = {d:.2f} mdeg ({name}, {t1:.2f} to {t2:.2f} min) → m_S = {m:.1f} ng/cm² ({model}, δ = {delta} nm, dn/dC = {dndc} cm³/g, S_SPR = {S} deg/RIU){note}'.format(
            d=delta_theta * 1000, name={'SPR angle': 'SPR', 'Bulk corrected': 'bulk corrected'}.get(trace_name, trace_name), t1=times[0], t2=times[1], m=coverage, delta=decay_length, dndc=dn_dC, S=S_SPR, note=note,
            model='thin film approximation' if thin_film else 'layer thickness d = {d:g} nm'.format(d=thickness))
        if trigger == 'quantification-sensorgram-graph':
            current_session.add_log('Surface coverage {file}: {result}'.format(file=current_data_path.split('/')[-1], result=result))
        return (points, result) + equation_outputs

    # Show the settings of the chosen TIR and SPR angle determination methods
    @dash.callback(
        dash.Output('TIR-fit-sg-options', 'style'),
        dash.Output('TIR-fit-legacy-options', 'style'),
        dash.Output('SPR-fit-dip-options', 'style'),
        dash.Output('SPR-fit-legacy-options', 'style'),
        dash.Input('TIR-fit-option-method', 'value'),
        dash.Input('SPR-fit-option-method', 'value'))
    def show_fit_method_options(TIR_method, SPR_method):
        shown, hidden = {'display': 'block'}, {'display': 'none'}
        return (shown if TIR_method == 'savitzky_golay' else hidden, hidden if TIR_method == 'savitzky_golay' else shown,
                shown if SPR_method == 'dip_level' else hidden, hidden if SPR_method == 'dip_level' else shown)

    # Show the TIR range chosen for a newly loaded measurement in the fit options
    @dash.callback(
        dash.Output('TIR-fit-option-range-low', 'value'),
        dash.Output('TIR-fit-option-range-high', 'value'),
        dash.Input('loaded-new-measurement', 'data'),
        prevent_initial_call=True)
    def update_TIR_range_inputs(signal):
        return current_session.SPR_TIR_fitting_parameters['TIR range'][0], current_session.SPR_TIR_fitting_parameters['TIR range'][1]

    @dash.callback(
        dash.Output('fresnel-analysis-bulk-RI-from-TIR', 'value'),
        dash.Output('fresnel-analysis-match-dip', 'value'),
        dash.Input({'type': 'fresnel-analysis-list', 'index': dash.ALL}, 'n_clicks'),
        dash.Input('add-fresnel-analysis-confirm', 'n_clicks'),
        prevent_initial_call=True)
    def update_bulk_RI_flag_checkbox(selected_fresnel_object, add_confirm):
        if 'add-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            return True, True
        if not dash.ctx.triggered or not dash.ctx.triggered[0]['value']:
            raise dash.exceptions.PreventUpdate
        analysis_object = current_session.fresnel_analysis_instances.get(dash.ctx.triggered_id['index'])
        if analysis_object is None:
            raise dash.exceptions.PreventUpdate
        return getattr(analysis_object, 'update_bulk_RI_from_TIR', True), getattr(analysis_object, 'match_dip_angle', False)

    @dash.callback(
        dash.Output('sensorgram-correction-layer-S_SPR', 'value'),
        dash.Output('sensorgram-correction-layer-decay-length', 'value'),
        dash.Input('loaded-new-measurement', 'data'),
        prevent_initial_call=True)
    def update_bulk_correction_parameters(signal):
        return instrument_SPR_sensitivity[current_data_path[-9:-6]], evanescent_decay_length[current_data_path[-9:-6]]


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
        dash.State('fresnel-analysis-bulk-RI-from-TIR', 'value'),
        dash.State('fresnel-analysis-match-dip', 'value'),
        prevent_initial_call=True)
    def update_reflectivity_fresnel_graph(run_model, add_button, add_confirm_button, remove_button, remove_confirm, remove_cancel, rangeslider_inp,
                                          selected_fresnel_object, save_png, save_svg, save_html, save_csv, rename_button, rename_confirm, batch_start_signal, analysis_name, figure_JSON, rangeslider_state, ini_guess,
                                          lower_bound, upper_bound,
                                          extinction_correction, analysis_name_, batch_files, background_sensors, batch_radio_selection, batch_newlayer_radio_selection, batch_sensor_index, batch_analysis_index, offset_fit_flag, elastomer_fit_flag, polarization_factor, bulk_RI_from_TIR_flag, match_dip_flag):

        global current_fresnel_analysis
        global current_data_path
        global current_session
        global reflectivity_df
        global current_sensor

        figure_object = go.Figure(figure_JSON)

        if 'fresnel-fit-option-rangeslider' == dash.ctx.triggered_id:

            # # First check if model has been run previously, then include model data before adding angle range lines
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
            new_figure.add_trace(go.Scatter(x=[rangeslider_inp[0], rangeslider_inp[0]],  # Adding angle range lines
                                            y=[min(list(figure_object.data[0]['y']['_inputArray'].values())[:-4]), max(list(figure_object.data[0]['y']['_inputArray'].values())[:-4])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            new_figure.add_trace(go.Scatter(x=[rangeslider_inp[1], rangeslider_inp[1]],  # Adding angle range lines
                                            y=[min(list(figure_object.data[0]['y']['_inputArray'].values())[:-4]), max(list(figure_object.data[0]['y']['_inputArray'].values())[:-4])],
                                            mode='lines',
                                            showlegend=False,
                                            line_color='black',
                                            line_dash='dash'
                                            ))
            # Updating layout
            new_figure.update_layout(xaxis_title=r'$\large{\text{Incident angle [ }^{\circ}\text{ ]}}$',
                                     yaxis_title=r'$\large{\text{Reflectivity [a.u.]}}$',
                                     font_family=PLOT_FONT,
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
            current_fresnel_analysis.update_bulk_RI_from_TIR = bool(bulk_RI_from_TIR_flag)
            current_fresnel_analysis.match_dip_angle = bool(match_dip_flag)

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

            # Remember which scan was fitted (used when repeating the model for all wavelengths)
            if current_fresnel_analysis.initial_data_path == current_data_path:
                current_fresnel_analysis.scan_index = SPRpy_multichannel.find_scan_index(ydata_df, current_fresnel_analysis.measurement_data['ydata'])

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
            current_session.add_log('FM{id} {name}: fitted {layer}|{param} = {value:.4f}{extra} (sensor S{sensor}, {file}, angle range {low:.3f}-{high:.3f} deg)'.format(
                id=current_fresnel_analysis.object_id, name=current_fresnel_analysis.name, layer=current_fresnel_analysis.fitted_layer,
                param=current_fresnel_analysis.sensor_object.optical_parameters.columns[current_fresnel_analysis.fitted_layer_index[1]],
                value=current_fresnel_analysis.fitted_result[0],
                extra=(', offset {o:.4f}'.format(o=current_fresnel_analysis.fitted_result[1]) if current_fresnel_analysis.fit_offset and len(current_fresnel_analysis.fitted_result) > 1 else '')
                      + (', prism k {k:.5f}'.format(k=current_fresnel_analysis.fitted_result[2]) if current_fresnel_analysis.fit_prism_k and len(current_fresnel_analysis.fitted_result) > 2 else '')
                      + ('' if getattr(current_fresnel_analysis, 'update_bulk_RI_from_TIR', True) else ', bulk RI from sensor table') + dip_log_text(current_fresnel_analysis),
                sensor=current_fresnel_analysis.sensor_object.object_id, file=current_fresnel_analysis.initial_data_path.split('/')[-1],
                low=current_fresnel_analysis.angle_range[0], high=current_fresnel_analysis.angle_range[1]), save=True)

            # Fit result text
            result = fit_result_text(current_fresnel_analysis)

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
                                     font_family=PLOT_FONT,
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

            table_header = [dash.html.Thead(dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td('FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td('S'+ str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name), dash.html.Td('{layer}|{parameter}-{channel}'.format(
                                             layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                                             parameter=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                                             channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)), dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            result_barplot_fig = fresnel_summary_barplot(current_session)

            return new_figure, fresnel_result_summary_dataframe, result_barplot_fig, 'finished', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, result, current_fresnel_analysis.sensor_object_label, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            current_fresnel_analysis = add_fresnel_model_object(current_session, current_sensor, current_data_path, reflectivity_df, analysis_name)
            current_session.add_log('Added fresnel analysis FM{id} {name} (sensor S{sensor}, {file})'.format(
                id=current_fresnel_analysis.object_id, name=analysis_name, sensor=current_sensor.object_id, file=current_data_path.split('/')[-1]))

            # Calculate initial intensity offset from data
            FR_y = fresnel_calculation(
                angles=reflectivity_df['angles'].iloc[[reflectivity_df['ydata'].to_numpy().argmin()]],
                wavelength=current_fresnel_analysis.sensor_object.wavelength,
                layer_thicknesses=current_fresnel_analysis.sensor_object.layer_thicknesses,
                n_re=current_fresnel_analysis.sensor_object.refractive_indices,
                n_im=current_fresnel_analysis.sensor_object.extinction_coefficients,
                ydata_type='R',
                polarization=current_fresnel_analysis.polarization)
            current_fresnel_analysis.y_offset = reflectivity_df['ydata'].min() - FR_y[0]

            current_fresnel_analysis.ini_guess = np.array([float(current_sensor.fitted_var), current_fresnel_analysis.y_offset])
            current_fresnel_analysis.bounds = [(current_fresnel_analysis.ini_guess[0] / 4, -np.inf), (current_fresnel_analysis.ini_guess[0] * 2, np.inf)]
            current_fresnel_analysis.angle_range = angle_range_around_minimum(reflectivity_df['angles'], reflectivity_df['ydata'],
                                                                             *current_session.SPR_TIR_fitting_parameters['Fresnel_angle_range_points'],
                                                                             min_angle=current_session.SPR_TIR_fitting_parameters['TIR range'][1])

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
                                     font_family=PLOT_FONT,
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
            old_path = current_session.location + '/Analysis instances/FM{id} {name}.pickle'.format(id=current_fresnel_analysis.object_id, name=current_fresnel_analysis.name)
            os.remove(old_path)

            # Change fresnel analysis name and save new fresnel analysis pickle file and session
            current_session.add_log('Renamed fresnel analysis FM{id} "{old}" to "{new}"'.format(id=current_fresnel_analysis.object_id, old=current_fresnel_analysis.name, new=analysis_name_))
            current_fresnel_analysis.name = analysis_name_
            current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
            current_session.save_session()

            analysis_options = [dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                                     id={'type': 'fresnel-analysis-list', 'index': fresnel_id},
                                                     n_clicks=0) for fresnel_id in current_session.fresnel_analysis_instances]

            table_header = [dash.html.Thead(dash.html.Tr(
                [dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td('FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td(
                    'S' + str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' +
                    current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name),
                 dash.html.Td('{layer}|{parameter}-{channel}'.format(
                     layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                     parameter=
                     current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[
                         current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                     channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                 dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for
                                           fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            result_barplot_fig = fresnel_summary_barplot(current_session)

            return dash.no_update, fresnel_result_summary_dataframe, result_barplot_fig, dash.no_update, analysis_options, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-fresnel-analysis-button' == dash.ctx.triggered_id:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-fresnel-analysis-confirm' == dash.ctx.triggered_id:
            if len(current_session.fresnel_analysis_instances) > 1:

                # Pop out the current fresnel analysis object from the session, delete its .pickle file and make the first instance the current one
                removed = current_fresnel_analysis
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
                current_session.add_log('Removed fresnel analysis FM{id} {name}'.format(id=removed.object_id, name=removed.name))
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
                    result = fit_result_text(current_fresnel_analysis)
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
                                         font_family=PLOT_FONT,
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

                table_header = [dash.html.Thead(dash.html.Tr(
                    [dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'),
                     dash.html.Th('Value')]))]
                table_body = [dash.html.Tbody([dash.html.Tr(
                    [dash.html.Td('FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td(
                        'S' + str(
                            current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' +
                        current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name),
                     dash.html.Td('{layer}|{parameter}-{channel}'.format(
                         layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                         parameter=current_session.fresnel_analysis_instances[
                             fresnel_inst].sensor_object.optical_parameters.columns[
                             current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                         channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                     dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))])
                                               for fresnel_inst in current_session.fresnel_analysis_instances])]
                fresnel_result_summary_dataframe = table_header + table_body

                result_barplot_fig = fresnel_summary_barplot(current_session)

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
                    TIR_angle, _, _, _, _ = TIR_determination(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'], current_session.SPR_TIR_fitting_parameters)
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
                                                                        file_path, next_reflectivity_df_, current_sensor.name)
                    # Calculate angle range based on measured data
                    current_fresnel_analysis.angle_range = angle_range_around_minimum(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'],
                                                                                     *current_session.SPR_TIR_fitting_parameters['Fresnel_angle_range_points'],
                                                                                     min_angle=current_session.SPR_TIR_fitting_parameters['TIR range'][1])

                    # Set analysis options from example analysis objects
                    current_fresnel_analysis.ini_guess = example_analysis_object.ini_guess
                    current_fresnel_analysis.bounds = example_analysis_object.bounds
                    current_fresnel_analysis.polarization = example_analysis_object.polarization
                    current_fresnel_analysis.extinction_correction = example_analysis_object.extinction_correction
                    current_fresnel_analysis.y_offset = example_analysis_object.y_offset
                    current_fresnel_analysis.fit_prism_k = example_analysis_object.fit_prism_k
                    current_fresnel_analysis.update_bulk_RI_from_TIR = getattr(example_analysis_object, 'update_bulk_RI_from_TIR', True)
                    current_fresnel_analysis.match_dip_angle = getattr(example_analysis_object, 'match_dip_angle', False)

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

                    current_session.add_log('Batch analysis: ' + 'FM{id} {name}: fitted {layer}|{param} = {value:.4f}{extra} (sensor S{sensor}, {file}, angle range {low:.3f}-{high:.3f} deg)'.format(
                            id=current_fresnel_analysis.object_id, name=current_fresnel_analysis.name, layer=current_fresnel_analysis.fitted_layer,
                            param=current_fresnel_analysis.sensor_object.optical_parameters.columns[current_fresnel_analysis.fitted_layer_index[1]],
                            value=current_fresnel_analysis.fitted_result[0],
                            extra=(', offset {o:.4f}'.format(o=current_fresnel_analysis.fitted_result[1]) if current_fresnel_analysis.fit_offset and len(current_fresnel_analysis.fitted_result) > 1 else '')
                                  + (', prism k {k:.5f}'.format(k=current_fresnel_analysis.fitted_result[2]) if current_fresnel_analysis.fit_prism_k and len(current_fresnel_analysis.fitted_result) > 2 else '')
                                  + ('' if getattr(current_fresnel_analysis, 'update_bulk_RI_from_TIR', True) else ', bulk RI from sensor table') + dip_log_text(current_fresnel_analysis),
                            sensor=current_fresnel_analysis.sensor_object.object_id, file=current_fresnel_analysis.initial_data_path.split('/')[-1],
                            low=current_fresnel_analysis.angle_range[0], high=current_fresnel_analysis.angle_range[1]))

                    # Save session and analysis object
                    current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
                    current_session.save_session()

            elif batch_radio_selection == 1:  # Use individual backgrounds and add new layer

                # Use the same layer structure copied from selected example sensor object
                for file_path, sensor_id in zip(batch_files, background_sensors):

                    # Load data from measurement file using load_csv_data
                    *_, next_reflectivity_df_ = load_csv_data(path=file_path)

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
                    TIR_angle, _, _, _, _ = TIR_determination(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'], current_session.SPR_TIR_fitting_parameters)
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
                                                                        current_sensor.name + ' (S' + str(
                                                                            current_sensor.object_id) + ') ')

                    # Calculate angle range based on measured data
                    current_fresnel_analysis.angle_range = angle_range_around_minimum(next_reflectivity_df_['angles'], next_reflectivity_df_['ydata'],
                                                                                     *current_session.SPR_TIR_fitting_parameters['Fresnel_angle_range_points'],
                                                                                     min_angle=current_session.SPR_TIR_fitting_parameters['TIR range'][1])

                    # Set analysis options from example analysis objects
                    current_fresnel_analysis.ini_guess = example_analysis_object.ini_guess
                    current_fresnel_analysis.bounds = example_analysis_object.bounds
                    current_fresnel_analysis.polarization = example_analysis_object.polarization
                    current_fresnel_analysis.extinction_correction = example_analysis_object.extinction_correction
                    current_fresnel_analysis.y_offset = example_analysis_object.y_offset
                    current_fresnel_analysis.fit_prism_k = example_analysis_object.fit_prism_k
                    current_fresnel_analysis.update_bulk_RI_from_TIR = getattr(example_analysis_object, 'update_bulk_RI_from_TIR', True)
                    current_fresnel_analysis.match_dip_angle = getattr(example_analysis_object, 'match_dip_angle', False)

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

                    current_session.add_log('Batch analysis: ' + 'FM{id} {name}: fitted {layer}|{param} = {value:.4f}{extra} (sensor S{sensor}, {file}, angle range {low:.3f}-{high:.3f} deg)'.format(
                            id=current_fresnel_analysis.object_id, name=current_fresnel_analysis.name, layer=current_fresnel_analysis.fitted_layer,
                            param=current_fresnel_analysis.sensor_object.optical_parameters.columns[current_fresnel_analysis.fitted_layer_index[1]],
                            value=current_fresnel_analysis.fitted_result[0],
                            extra=(', offset {o:.4f}'.format(o=current_fresnel_analysis.fitted_result[1]) if current_fresnel_analysis.fit_offset and len(current_fresnel_analysis.fitted_result) > 1 else '')
                                  + (', prism k {k:.5f}'.format(k=current_fresnel_analysis.fitted_result[2]) if current_fresnel_analysis.fit_prism_k and len(current_fresnel_analysis.fitted_result) > 2 else '')
                                  + ('' if getattr(current_fresnel_analysis, 'update_bulk_RI_from_TIR', True) else ', bulk RI from sensor table') + dip_log_text(current_fresnel_analysis),
                            sensor=current_fresnel_analysis.sensor_object.object_id, file=current_fresnel_analysis.initial_data_path.split('/')[-1],
                            low=current_fresnel_analysis.angle_range[0], high=current_fresnel_analysis.angle_range[1]))

                    # Save session and analysis object
                    current_session.save_fresnel_analysis(current_fresnel_analysis.object_id)
                    current_session.save_session()

            # Fit result text
            result = fit_result_text(current_fresnel_analysis)

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
                                     font_family=PLOT_FONT,
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

            table_header = [dash.html.Thead(dash.html.Tr(
                [dash.html.Th('Analysis'), dash.html.Th('Sensor'), dash.html.Th('Variable'), dash.html.Th('Value')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td('FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name), dash.html.Td(
                    'S' + str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' +
                    current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name),
                 dash.html.Td('{layer}|{parameter}-{channel}'.format(
                     layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                     parameter=
                     current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[
                         current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                     channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel)),
                 dash.html.Td(round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3))]) for
                                           fresnel_inst in current_session.fresnel_analysis_instances])]
            fresnel_result_summary_dataframe = table_header + table_body

            result_barplot_fig = fresnel_summary_barplot(current_session)

            return new_figure, fresnel_result_summary_dataframe, result_barplot_fig, dash.no_update, analysis_options, dash.no_update, dash.no_update, dash.no_update, current_fresnel_analysis.angle_range, current_fresnel_analysis.ini_guess[0], lower_bound_, upper_bound_, dash.no_update, result, current_fresnel_analysis.sensor_object_label, dash.no_update, dash.no_update, current_fresnel_analysis.measurement_data['angles'].iloc[0].astype('int'), current_fresnel_analysis.measurement_data['angles'].iloc[-1].astype('int')+1, 'Data path: \n' + current_fresnel_analysis.initial_data_path, dash.no_update, 'finished', dash.no_update, dash.no_update

        elif 'fresnel-reflectivity-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure_object, save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(figure_object, save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(figure_object, save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'fresnel-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
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
                result = fit_result_text(current_fresnel_analysis)
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
                                     font_family=PLOT_FONT,
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
        dash.Output('quantification-TIR-SPR-fit-collapse', 'is_open'),
        dash.Output('quantification-TIR-fit-graph', 'figure'),
        dash.Output('quantification-SPR-fit-graph', 'figure'),
        dash.Output('loaded-new-measurement', 'data'),
        dash.Input('quantification-show-SPR-TIR-fit-options-switch', 'value'),
        dash.Input('quantification-sensorgram-graph', 'hoverData'),
        dash.Input('quantification-apply-fitting-SPR-TIR-button', 'n_clicks'),
        dash.State('hover-selection-switch', 'value'),
        dash.State('TIR-fit-option-range-low', 'value'),
        dash.State('TIR-fit-option-range-high', 'value'),
        dash.State('TIR-fit-option-window', 'value'),
        dash.State('TIR-fit-option-points', 'value'),
        dash.State('TIR-fit-option-below-peak', 'value'),
        dash.State('TIR-fit-option-above-peak', 'value'),
        dash.State('SPR-fit-option-points', 'value'),
        dash.State('SPR-fit-option-below-peak', 'value'),
        dash.State('SPR-fit-option-above-peak', 'value'),
        dash.State('quantification-show-SPR-TIR-fit-options-switch', 'value'),
        dash.State('TIR-fit-option-method', 'value'),
        dash.State('TIR-fit-option-smooth-width', 'value'),
        dash.State('TIR-fit-option-fit-halfwidth', 'value'),
        dash.State('SPR-fit-option-method', 'value'),
        dash.State('SPR-fit-option-dip-level', 'value'),
        prevent_initial_call=True)
    def SPR_TIR_fitting_parameters_update(fit_show_switch, hoverData, run_button, hover_selection_switch, TIR_range_low, TIR_range_high, TIR_window, TIR_fit_points, TIR_below_peak, TIR_above_peak, SPR_fit_points, SPR_below_peak, SPR_above_peak, fit_show_switch_state,
                                          TIR_method, TIR_smooth_width, TIR_fit_halfwidth, SPR_method, SPR_dip_level):

        global current_session
        global sensorgram_df
        global sensorgram_df_selection
        global corrected_sensorgram_df_selection
        global ydata_df
        global reflectivity_df
        global scanspeed

        if 'quantification-show-SPR-TIR-fit-options-switch' == dash.ctx.triggered_id:
            return fit_show_switch, dash.no_update, dash.no_update, dash.no_update

        # Applying the fit settings and updating  the session object
        elif 'quantification-apply-fitting-SPR-TIR-button' == dash.ctx.triggered_id:
            current_session.SPR_TIR_fitting_parameters['TIR range'] = [float(TIR_range_low), float(TIR_range_high)]
            current_session.SPR_TIR_fitting_parameters['TIR window count'] = int(TIR_window)
            current_session.SPR_TIR_fitting_parameters['TIR fit points'] = int(TIR_fit_points)
            current_session.SPR_TIR_fitting_parameters['points_below_TIR_peak'] = int(TIR_below_peak)
            current_session.SPR_TIR_fitting_parameters['points_above_TIR_peak'] = int(TIR_above_peak)
            current_session.SPR_TIR_fitting_parameters['SPR fit points'] = int(SPR_fit_points)
            current_session.SPR_TIR_fitting_parameters['sensorgram_angle_range_points'] = [int(SPR_below_peak), int(SPR_above_peak)]
            current_session.SPR_TIR_fitting_parameters['TIR method'] = TIR_method
            current_session.SPR_TIR_fitting_parameters['TIR smooth width deg'] = float(TIR_smooth_width)
            current_session.SPR_TIR_fitting_parameters['TIR fit halfwidth deg'] = float(TIR_fit_halfwidth)
            current_session.SPR_TIR_fitting_parameters['SPR method'] = SPR_method
            current_session.SPR_TIR_fitting_parameters['SPR dip level'] = float(SPR_dip_level)
            fit_parameters = current_session.SPR_TIR_fitting_parameters
            if TIR_method == 'savitzky_golay':
                TIR_settings = 'filter width {w} deg, peak fit half-width {h} deg'.format(w=fit_parameters['TIR smooth width deg'], h=fit_parameters['TIR fit halfwidth deg'])
            else:
                TIR_settings = 'window {w}, points below/above peak {b}/{a}'.format(w=fit_parameters['TIR window count'], b=fit_parameters['points_below_TIR_peak'], a=fit_parameters['points_above_TIR_peak'])
            if SPR_method == 'dip_level':
                SPR_settings = 'dip level {d}'.format(d=fit_parameters['SPR dip level'])
            else:
                SPR_settings = 'points below/above dip {b}/{a}'.format(b=fit_parameters['sensorgram_angle_range_points'][0], a=fit_parameters['sensorgram_angle_range_points'][1])
            current_session.add_log('Recalculated sensorgram with TIR method {TIR_method} (range {TIR_range} deg, {TIR_settings}) and SPR method {SPR_method} ({SPR_settings})'.format(
                TIR_method=TIR_method, TIR_range=fit_parameters['TIR range'], TIR_settings=TIR_settings, SPR_method=SPR_method, SPR_settings=SPR_settings))

            # Select active TIR fitting parameters based on scanspeed
            if scanspeed <= 5:
                current_session.SPR_TIR_fitting_parameters['window_count_scanspeeds_1_5'] = int(TIR_window)
                current_session.SPR_TIR_fitting_parameters['points_above_TIR_peak_scanspeed_1_5'] = int(TIR_above_peak)
                current_session.SPR_TIR_fitting_parameters['points_below_TIR_peak_scanspeed_1_5'] = int(TIR_below_peak)
            else:
                current_session.SPR_TIR_fitting_parameters['window_count_scanspeeds_10'] = int(TIR_window)
                current_session.SPR_TIR_fitting_parameters['points_above_TIR_peak_scanspeed_10'] = int(TIR_above_peak)
                current_session.SPR_TIR_fitting_parameters['points_below_TIR_peak_scanspeed_10'] = int(TIR_below_peak)

            current_session.save_session()

            sensorgram_df = calculate_sensorgram(time_df, angles_df, ydata_df, current_session.SPR_TIR_fitting_parameters)

            # Offset to start at 0 degrees at 0 minutes
            sensorgram_df_selection = copy.deepcopy(sensorgram_df)
            sensorgram_df_selection['SPR angle'] = sensorgram_df_selection['SPR angle'] - \
                                                   sensorgram_df_selection['SPR angle'][0]
            sensorgram_df_selection['TIR angle'] = sensorgram_df_selection['TIR angle'] - \
                                                   sensorgram_df_selection['TIR angle'][0]

            # Calculate bulk correction
            corrected_sensorgram_df_selection = sensorgram_df_selection['SPR angle'] - sensorgram_df_selection[
                'TIR angle'] * instrument_SPR_sensitivity[current_data_path[-9:-6]] / instrument_TIR_sensitivity * math.exp(-2 * 0 / evanescent_decay_length[current_data_path[-9:-6]])

            return dash.no_update, dash.no_update, dash.no_update, 'signal'

        # When hovering over data in the sensorgram plot, update the TIR and SPR fitting graphs accordingly
        elif 'quantification-sensorgram-graph' == dash.ctx.triggered_id:

            if not hover_selection_switch and fit_show_switch_state:
                time_index = hoverData['points'][0]['pointIndex']
                fit_details = sensorgram_fit_details(angles_df, ydata_df.iloc[time_index, :], current_session.SPR_TIR_fitting_parameters)
                TIR_fitting_figure = make_TIR_fit_figure(fit_details)
                SPR_fitting_figure = make_SPR_fit_figure(fit_details)

                return dash.no_update, TIR_fitting_figure, SPR_fitting_figure, dash.no_update
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
        dash.Output('exclusion-height-option-RI-initialguess', 'value'),
        dash.Output('exclusion-height-option-RI-lowerbound', 'value'),
        dash.Output('exclusion-height-option-RI-upperbound', 'value'),
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
        dash.State('exclusion-height-option-RI-initialguess', 'value'),
        dash.State('exclusion-height-option-RI-lowerbound', 'value'),
        dash.State('exclusion-height-option-RI-upperbound', 'value'),
        dash.State('exclusion-height-bulk-RI-method', 'value'),
        prevent_initial_call=True)
    def exclusion_height_analysis_control(add_exclusion, confirm_exclusion, choose_exclusion, remove_analysis_button,
                                        remove_confirm, remove_cancel, clickData, clear_points, sensorgram_png,
                                        sensorgram_svg, sensorgram_html, sensorgram_csv, SPRvsTIR_png, SPRvsTIR_svg, SPRvsTIR_html, SPRvsTIR_csv,
                                        reflectivity_save_png, reflectivity_save_svg, reflectivity_save_html,reflectivity_save_csv,
                                        dnpair_save_png, dnpair_save_svg, dnpair_save_html, dnpair_save_csv, active_page, dnpair_hoverdata, initialize_model, analysis_name,
                                        action_selected, background_selected_id, sensorgram_figure_JSON, SPRvsTIR_figure_JSON, reflectivity_figure_JSON,
                                        dnpair_figure_JSON, active_page_state, RI_initial_guess_state, lower_RI_bound_state, upper_RI_bound_state, bulk_RI_method_state):
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
                                                 font_family=PLOT_FONT,
                                                 font_size=19,
                                                 margin_r=25,
                                                 margin_l=60,
                                                 margin_t=40,
                                                 template='simple_white',
                                                 uirevision=True)
                    updated_figure.update_xaxes(mirror=True, showline=True)
                    updated_figure.update_yaxes(mirror=True, showline=True)

                    return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

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
                                         font_family=PLOT_FONT,
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

            injection_time_string = '{length} injection points: {points}'.format(length=len(injection_points_time), points=injection_points_time_)

            buffer_time_string = '{length} buffer points: {points}'.format(length=len(buffer_points_time), points=buffer_points_time_)

            probe_time_string = '{length} probe points: {points}'.format(length=len(probe_points_time), points=probe_points_time_)

            current_session.save_session()
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)

            return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, injection_time_string, buffer_time_string, probe_time_string, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

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
                                         font_family=PLOT_FONT,
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

            return updated_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, injection_time_string, buffer_time_string, probe_time_string, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

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
                                                   font_family=PLOT_FONT,
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

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, mean_reflectivity_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-exclusion-height-analysis-button' == dash.ctx.triggered_id:
            # Open add analysis name giving modal
            return dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'add-exclusion-height-analysis-confirm' == dash.ctx.triggered_id:

            # Add new exclusion height analysis object to session
            background_object = current_session.fresnel_analysis_instances[background_selected_id]
            current_exclusion_height_analysis = add_exclusion_height_object(current_session, background_object, sensorgram_df_selection, current_data_path, analysis_name)
            current_exclusion_height_analysis.RI_initial_guess = RI_initial_guess_state
            current_exclusion_height_analysis.RI_bounds = [lower_RI_bound_state, upper_RI_bound_state]
            current_session.add_log('Added exclusion height analysis EH{id} {name} (fresnel background FM{fm}, {file})'.format(
                id=current_exclusion_height_analysis.object_id, name=analysis_name, fm=background_object.object_id, file=current_data_path.split('/')[-1]))
            current_session.save_session()
            current_session.save_exclusion_height_analysis(current_exclusion_height_analysis.object_id)

            # Calculate suggestions of lower and upper bounds for height
            lower_height_bound = float(background_object.sensor_object.layer_thicknesses[-2])
            upper_height_bound = max(300.0, float(background_object.sensor_object.layer_thicknesses[-2]) * 12)  # Wide range, cheap with the adaptive search
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
                                             font_family=PLOT_FONT,
                                             font_size=19,
                                             margin_r=25,
                                             margin_l=60,
                                             margin_t=40,
                                             template='simple_white',
                                             uirevision=True)
            new_sensorgram_fig.update_xaxes(mirror=True, showline=True)
            new_sensorgram_fig.update_yaxes(mirror=True, showline=True)

            return new_sensorgram_fig, dash.no_update, False, analysis_options, dash.no_update, background_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, False, True, False, 'Mean exclusion height: None', 'Mean exclusion RI: None', 'All exclusion heights: None', 'All exclusion RI: None', dash.no_update, dash.no_update, dash.no_update, dash.no_update, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, 100, 'Selected analysis: ' + analysis_name_, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-exclusion-height-analysis-button' == dash.ctx.triggered_id:
            # Open remove analysis object confirmation modal
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-exclusion-height-analysis-confirm' == dash.ctx.triggered_id:

            if len(current_session.exclusion_height_analysis_instances) > 1:

                # Pop out the current exclusion height analysis object from the session, delete its .pickle file and make the first instance the current one
                removed = current_exclusion_height_analysis
                try:
                    current_exclusion_height_analysis = current_session.exclusion_height_analysis_instances[1]
                except KeyError:  # In case the first few instances have already been removed, try the next one
                    failed = True
                    attempted = 1
                    max_attempts = len(current_session.exclusion_height_analysis_instances) + 1
                    while failed and not attempted == max_attempts:
                        attempted += 1
                        try:
                            current_exclusion_height_analysis = current_session.exclusion_height_analysis_instances[attempted+1]
                        except KeyError:
                            continue
                        failed = False
                current_session.remove_exclusion_height_analysis(removed.object_id)
                current_session.add_log('Removed exclusion height analysis EH{id} {name}'.format(id=removed.object_id, name=removed.name))
                current_session.save_session()

                # Lower and upper bounds for height and RI
                lower_height_bound = current_exclusion_height_analysis.height_bounds[0]
                upper_height_bound = current_exclusion_height_analysis.height_bounds[1]
                resolution = current_exclusion_height_analysis.d_n_pair_resolution
                RI_initial_guess = current_exclusion_height_analysis.RI_initial_guess
                lower_RI_bound = current_exclusion_height_analysis.RI_bounds[0]
                upper_RI_bound = current_exclusion_height_analysis.RI_bounds[1]

                analysis_name_ = 'EH' + str(
                    current_exclusion_height_analysis.object_id) + ' ' + current_exclusion_height_analysis.name

                # Update choose analysis dropdown menu options
                analysis_options = [dbc.DropdownMenuItem(
                    'EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                    id={'type': 'exclusion-analysis-list',
                        'index': exclusion_id},
                    n_clicks=0) for exclusion_id in current_session.exclusion_height_analysis_instances]

                # Update results text
                mean_result_height, mean_result_RI, all_result_heights, all_result_RI = exclusion_height_result_texts(current_exclusion_height_analysis)

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
                                                 font_family=PLOT_FONT,
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
                                              font_family=PLOT_FONT,
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
                                              font_family=PLOT_FONT,
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
                    font_family=PLOT_FONT,
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
                    dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean'), dash.html.Th('Exclusion height all'), dash.html.Th('Exclusion RI mean'), dash.html.Th('Exclusion RI all')]))]
                table_body = [dash.html.Tbody([dash.html.Tr(
                    [dash.html.Td(current_session.exclusion_height_analysis_instances[exclusion_inst].name),
                     dash.html.Td('{mean_} {all_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[0], 2), all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2)))),
                     dash.html.Td('{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2)))),
                     dash.html.Td('{mean_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0], 4))),
                     dash.html.Td('{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :], decimals=4))))])
                    for exclusion_inst in current_session.exclusion_height_analysis_instances])]
                exclusion_result_summary_dataframe = table_header + table_body

                return new_sensorgram_fig, exclusion_result_summary_dataframe, False, analysis_options, False, current_exclusion_height_analysis.fresnel_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, True, True, True, mean_result_height, mean_result_RI, all_result_heights, all_result_RI, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, num_injection_steps, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, resolution, 'Selected analysis: ' + analysis_name_, RI_initial_guess, lower_RI_bound, upper_RI_bound

            else:
                try:
                    current_session.remove_exclusion_height_analysis(current_exclusion_height_analysis.object_id)
                except AttributeError:
                    pass  # There was no object at all

                current_exclusion_height_analysis = None
                current_session.save_session()

                table_header = [dash.html.Thead(
                    dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean'), dash.html.Th('Exclusion height all'),
                                  dash.html.Th('Exclusion RI mean'), dash.html.Th('Exclusion RI all')]))]
                table_body = [dash.html.Tbody([dash.html.Tr([dash.html.Td(''), dash.html.Td(''), dash.html.Td(''), dash.html.Td(''), dash.html.Td('')])])]
                exclusion_result_summary_dataframe = table_header + table_body

                return dash.no_update, exclusion_result_summary_dataframe, dash.no_update, '', False, 'Sensor: None', 'Fresnel background: None', False, False, False, False, 'Mean exclusion height: None', 'Mean exclusion RI: None', 'All exclusion heights: None', 'All exclusion RI: None', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, '', '', '', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'remove-exclusion-height-analysis-cancel' == dash.ctx.triggered_id:
            # Cancel removal of exclusion height analysis object

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'exclusion-height-SPRvsTIR-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(SPRvsTIR_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(SPRvsTIR_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(SPRvsTIR_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-SPRvsTIR-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
            SPRvsTIR_fig = go.Figure(SPRvsTIR_figure_JSON)
            fig_df = pd.DataFrame(data={'TIR': list(SPRvsTIR_fig.data[0].x["_inputArray"].values())[:-3], 'SPR': list(SPRvsTIR_fig.data[0].y["_inputArray"].values())[:-3]})
            fig_df.to_csv(save_filename, sep=';')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-png' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(reflectivity_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(reflectivity_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(reflectivity_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-reflectivity-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
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
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(sensorgram_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(sensorgram_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(sensorgram_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-sensorgram-save-csv' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
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
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('PNG files', '*.png')], default_extension='.png')
            plotly.io.write_image(go.Figure(dnpair_figure_JSON), save_filename, format='png')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-svg' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('SVG files', '*.svg')], default_extension='.svg')
            plotly.io.write_image(go.Figure(dnpair_figure_JSON), save_filename, format='svg')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-html' == dash.ctx.triggered_id:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(go.Figure(dnpair_figure_JSON), save_filename, include_mathjax='cdn')
            raise dash.exceptions.PreventUpdate

        elif 'exclusion-height-d-n-pair-save-csv' == dash.ctx.triggered_id:
            d_n_pair_fig_csv = go.Figure(dnpair_figure_JSON)
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')
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
                                          font_family=PLOT_FONT,
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
                                                   font_family=PLOT_FONT,
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
                    font_family=PLOT_FONT,
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

            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        elif 'exclusion-height-initialize-model' == dash.ctx.triggered_id:

            # Check that appropriate points have been selected
            if len(current_exclusion_height_analysis.buffer_points) % 4 == 0 and len(current_exclusion_height_analysis.injection_points) % 2 == 0 and len(current_exclusion_height_analysis.probe_points) % 2 == 0:

                # Initializes model parameters and attributes prepping for running. Also activate run buttons and result page.
                current_exclusion_height_analysis.initialize_model(ydata_df, bulk_RI_method=bulk_RI_method_state or 'edge_overlay')
                current_session.add_log('EH{id} {name}: initialized model with {n} injection(s), TIR shifts buffer to probe {shifts} mdeg ({method}), buffer bulk RI {buffer}, probe bulk RI {probe}'.format(
                    id=current_exclusion_height_analysis.object_id, name=current_exclusion_height_analysis.name,
                    n=len(current_exclusion_height_analysis.probe_points) // 2,
                    shifts=[round(float(shift) * 1000, 1) for shift in current_exclusion_height_analysis.TIR_shifts],
                    method=EXCLUSION_BULK_RI_METHODS[current_exclusion_height_analysis.bulk_RI_method],
                    buffer=[round(float(ri), 5) for ri in current_exclusion_height_analysis.buffer_bulk_RIs],
                    probe=[round(float(ri), 5) for ri in current_exclusion_height_analysis.probe_bulk_RIs]))

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
                                              font_family=PLOT_FONT,
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
                                                       font_family=PLOT_FONT,
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

                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True, dash.no_update, True, dash.no_update, dash.no_update, dash.no_update, dash.no_update, SPRvsTIR_figure, mean_reflectivity_figure, dash.no_update, num_injection_steps, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

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

            # Lower and upper bounds for height and RI fitting
            lower_height_bound = current_exclusion_height_analysis.height_bounds[0]
            upper_height_bound = current_exclusion_height_analysis.height_bounds[1]
            resolution = current_exclusion_height_analysis.d_n_pair_resolution
            RI_initial_guess = current_exclusion_height_analysis.RI_initial_guess
            lower_RI_bound = current_exclusion_height_analysis.RI_bounds[0]
            upper_RI_bound = current_exclusion_height_analysis.RI_bounds[1]

            analysis_name_ = 'EH' + str(
                current_exclusion_height_analysis.object_id) + ' ' + current_exclusion_height_analysis.name

            # Update choose analysis dropdown menu options
            analysis_options = [dbc.DropdownMenuItem(
                'EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                id={'type': 'exclusion-analysis-list',
                    'index': exclusion_id},
                n_clicks=0) for exclusion_id in current_session.exclusion_height_analysis_instances]

            # Update results text
            mean_result_height, mean_result_RI, all_result_heights, all_result_RI = exclusion_height_result_texts(current_exclusion_height_analysis)

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
                                             font_family=PLOT_FONT,
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
                                          font_family=PLOT_FONT,
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
                                          font_family=PLOT_FONT,
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
                font_family=PLOT_FONT,
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

            return new_sensorgram_fig, dash.no_update, False, analysis_options, False, current_exclusion_height_analysis.fresnel_object.sensor_object_label, current_exclusion_height_analysis.fresnel_object_label, True, False, True, True, mean_result_height, mean_result_RI, all_result_heights, all_result_RI, SPRvsTIR_figure, mean_reflectivity_figure, d_n_pair_figure, num_injection_steps, lower_height_bound, upper_height_bound, injection_time_string, buffer_time_string, probe_time_string, resolution, 'Selected analysis: ' + analysis_name_, RI_initial_guess, lower_RI_bound, upper_RI_bound

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
        dash.State('exclusion-height-option-RI-lowerbound', 'value'),
        dash.State('exclusion-height-option-RI-upperbound', 'value'),
        dash.State('exclusion-height-option-RI-initialguess', 'value'),
        dash.State('exclusion-height-option-resolution', 'value'),
        dash.State('exclusion-height-analysis-offset-refit', 'value'),
        dash.State('exclusion-height-analysis-prism-refit', 'value'),
        dash.State('exclusion-height-search-method', 'value'),
        dash.State('exclusion-height-match-dip', 'value'),
        prevent_initial_call=True)
    def run_exclusion_height_calculations(run_button, lower_height_bound, upper_height_bound, lower_RI_bound, upper_RI_bound, RI_initial_guess, resolution, fit_offset_flag, fit_prism_flag, height_search, match_dip_flag):
        """
        This callback runs the exclusion height calculations in the background. It is triggered by the run button, and
        updates the progress bar. It also updates the result figures when the calculations are done.

        :param set_progress: defines the progress bar value and max value
        :param run_button:
        :param check_button:
        :param abort_button:
        :param lower_height_bound:
        :param upper_height_bound:
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
            current_exclusion_height_analysis.height_search = height_search or 'adaptive'
            current_exclusion_height_analysis.match_dip_angle = bool(match_dip_flag)

            # Set resolution, height steps and RI bounds and guess
            current_exclusion_height_analysis.d_n_pair_resolution = resolution
            current_exclusion_height_analysis.height_bounds[0] = lower_height_bound
            current_exclusion_height_analysis.height_bounds[1] = upper_height_bound
            current_exclusion_height_analysis.height_steps = np.linspace(lower_height_bound, upper_height_bound, resolution)
            current_exclusion_height_analysis.RI_bounds = [lower_RI_bound, upper_RI_bound]
            current_exclusion_height_analysis.RI_initial_guess = RI_initial_guess

            # Run exclusion height calculations (overwrites previous results) and calculate the mean results
            run_exclusion_height_calculation(current_exclusion_height_analysis, logical_cores)

            mean_result_height, mean_result_RI, all_result_heights, all_result_RI = exclusion_height_result_texts(current_exclusion_height_analysis)

            current_session.add_log('EH{id} {name}: {height}, {RI}, {all_heights} (heights {heights_from}-{heights_to} nm, {resolution} steps, {search} search, RI bounds {RI_bounds}, fit offset {offset}, fit prism k {prism}, dip angle matched {dip})'.format(
                id=current_exclusion_height_analysis.object_id, name=current_exclusion_height_analysis.name,
                height=mean_result_height, RI=mean_result_RI, all_heights=all_result_heights, heights_from=lower_height_bound, heights_to=upper_height_bound, resolution=resolution,
                RI_bounds=current_exclusion_height_analysis.RI_bounds, offset=fit_offset_flag, prism=fit_prism_flag,
                search=current_exclusion_height_analysis.height_search, dip=current_exclusion_height_analysis.match_dip_angle))

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
                                            font_family=PLOT_FONT,
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
                                                   font_family=PLOT_FONT,
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
                font_family=PLOT_FONT,
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
                dash.html.Tr([dash.html.Th('Analysis'), dash.html.Th('Exclusion height mean'),
                              dash.html.Th('Exclusion height all'), dash.html.Th('Exclusion RI mean'),
                              dash.html.Th('Exclusion RI all')]))]
            table_body = [dash.html.Tbody([dash.html.Tr(
                [dash.html.Td(current_session.exclusion_height_analysis_instances[exclusion_inst].name),
                 dash.html.Td('{mean_} {all_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[
                         0], 2), all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :],
                     decimals=2)))),
                 dash.html.Td('{all_}'.format(all_=str(np.round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :],
                     decimals=2)))),
                 dash.html.Td('{mean_}'.format(mean_=round(
                     current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0],
                     4))),
                 dash.html.Td('{all_}'.format(all_=str(np.round(
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

    channel_data_cache = {}  # Loaded data of the channels of the measurement {path: MeasurementChannel}

    def cached_channels():
        # Reuse the data of the current measurement file
        if current_data_path not in channel_data_cache:
            channel = SPRpy_multichannel.MeasurementChannel.__new__(SPRpy_multichannel.MeasurementChannel)
            channel.path, channel.scanspeed, channel.time_df, channel.angles_df, channel.ydata_df = current_data_path, scanspeed, time_df, angles_df, ydata_df
            channel.label = SPRpy_multichannel.channel_label(current_data_path)
            channel_data_cache[current_data_path] = channel
        return channel_data_cache

    @dash.callback(
        dash.Output('fresnel-all-wavelengths-channels-hint', 'children'),
        dash.Output('fresnel-all-wavelengths-channels', 'options'),
        dash.Output('fresnel-all-wavelengths-channels', 'value'),
        dash.Output('fresnel-all-wavelengths-result', 'children'),
        dash.Output('fresnel-analysis-dropdown', 'children', allow_duplicate=True),
        dash.Output('exclusion-choose-background-dropdown', 'options', allow_duplicate=True),
        dash.Output('chosen-sensor-dropdown', 'children', allow_duplicate=True),
        dash.Output('result-summary-fresnel-table', 'children', allow_duplicate=True),
        dash.Output('summary-fresnel-barplot-graph', 'figure', allow_duplicate=True),
        dash.Input('fresnel-fit-result', 'children'),
        dash.Input('fresnel-all-wavelengths-run', 'n_clicks'),
        dash.State('fresnel-all-wavelengths-channels', 'value'),
        prevent_initial_call=True)
    def fresnel_all_wavelengths(fit_result, run_clicks, selected_paths):
        """
        Repeats the current fresnel model for the other channels (wavelengths) of the measurement.
        """
        template = globals().get('current_fresnel_analysis')
        no_change = [dash.no_update] * 5

        if 'fresnel-fit-result' == dash.ctx.triggered_id:
            # Channels of the measurement of the selected analysis
            if template is None:
                raise dash.exceptions.PreventUpdate
            options, values = channel_options(template.initial_data_path, template.initial_data_path, channel_spots)
            hint = '' if len(options) > 1 else 'No other channels found for {file} (expected file names like "name-L2_980nm.csv" in the same folder).'.format(
                file=os.path.basename(template.initial_data_path))
            return [hint, options, values, dash.no_update] + no_change

        if template is None or template.fitted_data is None:
            return [dash.no_update, dash.no_update, dash.no_update, dbc.Alert('Run the fresnel model of this analysis first.', color='warning')] + no_change

        paths = [template.initial_data_path] + [path for path in (selected_paths or []) if path != template.initial_data_path]
        try:
            fresnel_objects, notes = SPRpy_multichannel.fit_fresnel_all_channels(current_session, template, default_sensor_values,
                                                                                 channel_paths=paths, channels=cached_channels())
        except (FileNotFoundError, ValueError, KeyError) as error:
            return [dash.no_update, dash.no_update, dash.no_update, dbc.Alert('Fitting all wavelengths failed: {e}'.format(e=error), color='danger')] + no_change
        current_session.save_session()

        parameter = template.sensor_object.optical_parameters.columns[template.fitted_layer_index[1]]
        rows = []
        for fm in by_spot(fresnel_objects, channel_spots):
            rows.append([spot_label(fm, channel_spots), fm.sensor_object.channel, 'FM{id} {name}'.format(id=fm.object_id, name=fm.name), 'S{id}'.format(id=fm.sensor_object.object_id),
                         round(float(fm.fitted_result[0]), 4),
                         round(float(fm.fitted_result[1]), 4) if fm.fit_offset and len(fm.fitted_result) > 1 else '-',
                         round(float(fm.fitted_result[2]), 5) if fm.fit_prism_k and len(fm.fitted_result) > 2 else '-',
                         round(float(fm.sensor_object.refractive_indices[-1]), 5),
                         '{a:.2f}-{b:.2f}'.format(a=fm.angle_range[0], b=fm.angle_range[1]),
                         ', '.join(SPRpy_multichannel.channel_fit_warnings(fm)) or '-'])
        spot_lines = SPRpy_multichannel.spot_summary_lines(fresnel_objects, channel_spots, lambda fm: float(fm.fitted_result[0]), decimals=4)
        current_session.add_log('Fresnel modelling for all wavelengths from FM{id} ({layer}|{parameter}, scan {scan}): {spots}'.format(
            id=template.object_id, layer=template.fitted_layer, parameter=parameter, scan=getattr(template, 'scan_index', '?'), spots='; '.join(spot_lines)), save=True)
        children = [small_table(['Spot', 'Channel', 'Analysis', 'Sensor', '{layer}|{parameter}'.format(layer=template.fitted_layer, parameter=parameter),
                                 'Offset', 'Prism k', 'Bulk RI', 'Angle range', 'Check'], rows),
                    dash.html.Ul([dash.html.Li(line) for line in spot_lines], style={'font-size': '0.9rem', 'margin-bottom': '4px'}),
                    dash.html.P('Scan {scan}. {note}'.format(scan=getattr(template, 'scan_index', '?'), note=WAVELENGTH_NOTE), style={'font-size': '0.8rem', 'color': '#6c757d'})]
        flagged = [fm.sensor_object.channel for fm in fresnel_objects if SPRpy_multichannel.channel_fit_warnings(fm)]
        if flagged:
            children.append(dbc.Alert('{channels}: the fit is at a bound or needs a large offset, so the model does not describe the measured dip '
                                      'well. Check the optical constants for this wavelength in config.toml (gold and prism in particular): the '
                                      'fit then compensates with a wrong layer thickness, which also affects exclusion heights calculated with this '
                                      'background.'.format(channels=', '.join(flagged)), color='warning', style={'font-size': '0.85rem'}))
        if notes:
            children.append(dash.html.Ul([dash.html.Li(note) for note in notes], style={'font-size': '0.8rem', 'color': '#6c757d'}))

        fresnel_dropdown = [dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                                 id={'type': 'fresnel-analysis-list', 'index': fresnel_id}, n_clicks=0)
                            for fresnel_id in current_session.fresnel_analysis_instances]
        background_options = [{'label': 'FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name, 'value': fresnel_id}
                              for fresnel_id in current_session.fresnel_analysis_instances]
        sensor_options = [dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                               id={'type': 'sensor-list', 'index': sensor_id}, n_clicks=0)
                          for sensor_id in current_session.sensor_instances]
        return [dash.no_update, dash.no_update, dash.no_update, children, fresnel_dropdown, background_options, sensor_options, fresnel_summary_table(current_session),
                fresnel_summary_barplot(current_session)]

    @dash.callback(
        dash.Output('exclusion-all-wavelengths-channels-hint', 'children'),
        dash.Output('exclusion-all-wavelengths-channels', 'options'),
        dash.Output('exclusion-all-wavelengths-channels', 'value'),
        dash.Output('exclusion-all-wavelengths-result', 'children'),
        dash.Output('exclusion-height-analysis-dropdown', 'children', allow_duplicate=True),
        dash.Output('result-summary-exclusion-table', 'children', allow_duplicate=True),
        dash.Input('exclusion-height-analysis-label', 'children'),
        dash.Input('exclusion-all-wavelengths-run', 'n_clicks'),
        dash.State('exclusion-all-wavelengths-channels', 'value'),
        dash.State('exclusion-height-option-lowerbound', 'value'),
        dash.State('exclusion-height-option-upperbound', 'value'),
        dash.State('exclusion-height-option-RI-lowerbound', 'value'),
        dash.State('exclusion-height-option-RI-upperbound', 'value'),
        dash.State('exclusion-height-option-RI-initialguess', 'value'),
        dash.State('exclusion-height-option-resolution', 'value'),
        dash.State('exclusion-height-analysis-offset-refit', 'value'),
        dash.State('exclusion-height-analysis-prism-refit', 'value'),
        dash.State('exclusion-height-bulk-RI-method', 'value'),
        dash.State('exclusion-height-search-method', 'value'),
        dash.State('exclusion-height-match-dip', 'value'),
        prevent_initial_call=True)
    def exclusion_height_all_wavelengths(analysis_label, run_clicks, selected_paths, lower_height_bound, upper_height_bound,
                                         lower_RI_bound, upper_RI_bound, RI_initial_guess, resolution, fit_offset_flag, fit_prism_flag, bulk_RI_method, height_search, match_dip_flag):
        """
        Runs the current exclusion height analysis for the other channels (wavelengths) of the measurement.
        """
        template = globals().get('current_exclusion_height_analysis')
        no_change = [dash.no_update] * 2

        if 'exclusion-height-analysis-label' == dash.ctx.triggered_id:
            # Channels of the measurement of the selected analysis
            if template is None:
                raise dash.exceptions.PreventUpdate
            options, values = channel_options(template.initial_data_path, template.initial_data_path, channel_spots)
            hint = '' if len(options) > 1 else 'No other channels found for {file} (expected file names like "name-L2_980nm.csv" in the same folder).'.format(
                file=os.path.basename(template.initial_data_path))
            return [hint, options, values, ''] + no_change

        if template is None:
            return [dash.no_update, dash.no_update, dash.no_update, dbc.Alert('Add an exclusion height analysis first.', color='warning')] + no_change
        if len(template.probe_points) == 0 or len(template.buffer_points) != 2 * len(template.probe_points) or len(template.injection_points) % 2 != 0:
            return [dash.no_update, dash.no_update, dash.no_update, dbc.Alert('Select injection, buffer and probe points first (two buffer intervals per probe interval).', color='warning')] + no_change

        # Same settings as the run button
        template.fit_offset = fit_offset_flag
        template.fit_prism = fit_prism_flag
        template.height_search = height_search or 'adaptive'
        template.match_dip_angle = bool(match_dip_flag)
        template.d_n_pair_resolution = resolution
        template.height_bounds = [lower_height_bound, upper_height_bound]
        template.height_steps = np.linspace(lower_height_bound, upper_height_bound, resolution)
        template.RI_bounds = [lower_RI_bound, upper_RI_bound]
        template.RI_initial_guess = RI_initial_guess

        paths = [template.initial_data_path] + [path for path in (selected_paths or []) if path != template.initial_data_path]
        try:
            exclusion_objects = SPRpy_multichannel.exclusion_height_all_channels(current_session, template, logical_cores, channel_paths=paths,
                                                                                 channels=cached_channels(), bulk_RI_method=bulk_RI_method or 'edge_overlay')
        except (FileNotFoundError, ValueError, KeyError, TypeError) as error:
            return [dash.no_update, dash.no_update, dash.no_update, dbc.Alert('Running all wavelengths failed: {e}'.format(e=error), color='danger')] + no_change
        current_session.save_session()

        rows, warnings = [], []
        for eh in by_spot(exclusion_objects, channel_spots):
            heights = np.asarray(eh.all_exclusion_results[0, :], dtype=np.float64)
            RIs = np.asarray(eh.all_exclusion_results[1, :], dtype=np.float64)
            sensitivity = np.asarray(getattr(eh, 'height_sensitivity', [np.nan]), dtype=np.float64)
            rows.append([spot_label(eh, channel_spots), eh.sensor_object.channel, 'EH{id} {name}'.format(id=eh.object_id, name=eh.name),
                         '{h:.2f} ± {s:.2f}'.format(h=eh.mean_exclusion_height_result[0], s=eh.mean_exclusion_height_result[1]),
                         str(np.round(heights, 2)), '{n:.4f}'.format(n=eh.mean_exclusion_RI_result[0]),
                         '{s:.2f}'.format(s=np.nanmean(sensitivity)) if np.any(np.isfinite(sensitivity)) else '-',
                         str(np.round(np.asarray(eh.TIR_shifts, dtype=np.float64) * 1000, 1))])
            at_bound = np.isclose(heights, eh.height_steps[0], atol=1e-6) | np.isclose(RIs, eh.RI_bounds[0], atol=1e-4) | np.isclose(RIs, eh.RI_bounds[1], atol=1e-4)
            if np.any(np.isnan(heights)) or np.any(at_bound):
                warnings.append('{channel}: no reliable crossing for {n} of {total} steps (no crossing within the height bounds, e.g. raise the upper '
                                'height bound, or height or RI at the bounds).'.format(channel=eh.sensor_object.channel, total=len(heights),
                                                                                       n=int(np.count_nonzero(np.isnan(heights) | at_bound))))
        spot_lines = SPRpy_multichannel.spot_summary_lines(exclusion_objects, channel_spots, lambda eh: float(eh.mean_exclusion_height_result[0]), unit='nm')
        current_session.add_log('Exclusion height for all wavelengths from EH{id}: {spots}'.format(id=template.object_id, spots='; '.join(spot_lines)), save=True)

        children = [small_table(['Spot', 'Channel', 'Analysis', 'Height [nm]', 'All heights [nm]', 'RI', 'nm per mdeg', 'TIR shifts [mdeg]'], rows),
                    dash.html.Ul([dash.html.Li(line) for line in spot_lines], style={'font-size': '0.9rem', 'margin-bottom': '4px'}),
                    dash.html.P('{note} "nm per mdeg" is how much the height changes per millidegree error in the TIR shift (lower is more robust). '
                                'The TIR shift only depends on the bulk solution, so it should be nearly the same in all channels. '
                                'Select an analysis above to see its curves.'.format(note=WAVELENGTH_NOTE), style={'font-size': '0.8rem', 'color': '#6c757d'})]
        if warnings:
            children.append(dbc.Alert(dash.html.Ul([dash.html.Li(warning) for warning in warnings]), color='warning', style={'font-size': '0.85rem'}))

        exclusion_dropdown = [dbc.DropdownMenuItem('EH' + str(exclusion_id) + ' ' + current_session.exclusion_height_analysis_instances[exclusion_id].name,
                                                   id={'type': 'exclusion-analysis-list', 'index': exclusion_id}, n_clicks=0)
                              for exclusion_id in current_session.exclusion_height_analysis_instances]
        return [dash.no_update, dash.no_update, dash.no_update, children, exclusion_dropdown, exclusion_summary_table(current_session)]

    def series_options_from_settings(options_values):
        return {'fit_offset': 'offset' in (options_values or []), 'fit_prism_k': 'prism' in (options_values or []),
                'bulk_from_TIR': 'bulk' in (options_values or []),
                'match_dip_angle': 'dip' in (options_values or []),
                'angle_range_points': current_session.SPR_TIR_fitting_parameters['Fresnel_angle_range_points'],
                'TIR_ranges': [TIR_default_parameters['TIR_range_water_or_long_measurement'], TIR_default_parameters['TIR_range_air_or_few_scans']]}

    @dash.callback(
        dash.Output('series-modal', 'is_open'),
        dash.Output('series-base-sensor', 'options'),
        dash.Output('series-base-sensor', 'value'),
        dash.Input('series-button', 'n_clicks'),
        dash.Input('series-close', 'n_clicks'),
        prevent_initial_call=True)
    def open_layer_series(open_clicks, close_clicks):
        if 'series-close' == dash.ctx.triggered_id:
            return False, dash.no_update, dash.no_update
        options = [{'label': 'S{id} {name}'.format(id=sensor_id, name=sensor.name), 'value': sensor_id} for sensor_id, sensor in current_session.sensor_instances.items()]
        current = globals().get('current_sensor')
        return True, options, current.object_id if current is not None and current.object_id in current_session.sensor_instances else (options[0]['value'] if options else None)

    @dash.callback(
        dash.Output('series-base-sensor-info', 'children'),
        dash.Input('series-base-sensor', 'value'),
        prevent_initial_call=True)
    def layer_series_base_sensor_info(sensor_id):
        sensor = current_session.sensor_instances.get(sensor_id)
        if sensor is None:
            return ''
        return 'Layers {layers}. The first step fits {layer}|{parameter} (select the fitted variable in the sensor table).'.format(
            layers=', '.join(str(layer) for layer in sensor.optical_parameters['Layers']),
            layer=sensor.optical_parameters.iloc[sensor.fitted_layer_index[0], 0], parameter=sensor.optical_parameters.columns[sensor.fitted_layer_index[1]])

    @dash.callback(
        dash.Output('series-files', 'data'),
        dash.Output('series-detected', 'children'),
        dash.Output('series-recipe', 'data'),
        dash.Input('series-choose-files', 'n_clicks'),
        prevent_initial_call=True)
    def choose_layer_series_files(n_clicks):
        paths = select_files('Select the measurement files of the layer series', prompt_folder=os.path.dirname(current_data_path), file_types=[('CSV files', '*.csv')])
        if not paths:
            raise dash.exceptions.PreventUpdate  # File dialog cancelled
        series = SPRpy_series.parse_series_files(paths)
        if not series['samples']:
            return None, 'None of the selected files match the file name pattern <date>-<time>_<sample>_<step>-L<channel>_<wavelength>nm.csv.', []
        channels = ', '.join('L{n} {w} nm'.format(n=n, w=w) for n, w in series['wavelengths'].items())
        detected = [dash.html.Div('{n} files: samples {samples}; steps (in measurement order) {steps}; channels {channels}.'.format(
            n=len(paths) - len(series['skipped']), samples=', '.join(series['samples']), steps=' → '.join(series['steps']), channels=channels))]
        incomplete = ['{sample} has no {step}'.format(sample=sample, step=step) for sample, steps in series['samples'].items()
                      for step in series['steps'] if step not in steps]
        if incomplete:
            detected.append(dash.html.Div('Missing: ' + '; '.join(incomplete) + ' (the next step then builds on the previous one).'))
        if series['skipped']:
            detected.append(dash.html.Div('Skipped (file name does not match): ' + ', '.join(series['skipped'])))
        return list(paths), detected, SPRpy_series.default_recipe(series, default_sensor_values)

    @dash.callback(
        dash.Output('series-result', 'children'),
        dash.Output('series-results', 'data'),
        dash.Output('series-export', 'disabled'),
        dash.Output('fresnel-analysis-dropdown', 'children', allow_duplicate=True),
        dash.Output('exclusion-choose-background-dropdown', 'options', allow_duplicate=True),
        dash.Output('chosen-sensor-dropdown', 'children', allow_duplicate=True),
        dash.Output('result-summary-fresnel-table', 'children', allow_duplicate=True),
        dash.Output('summary-fresnel-barplot-graph', 'figure', allow_duplicate=True),
        dash.Input('series-run', 'n_clicks'),
        dash.State('series-files', 'data'),
        dash.State('series-base-sensor', 'value'),
        dash.State('series-recipe', 'data'),
        dash.State('series-options', 'value'),
        prevent_initial_call=True)
    def run_layer_series(n_clicks, paths, sensor_id, recipe, options_values):
        no_change = [dash.no_update] * 5
        if not paths:
            return [dbc.Alert('Choose the measurement files first.', color='warning'), dash.no_update, dash.no_update] + no_change
        template_sensor = current_session.sensor_instances.get(sensor_id)
        if template_sensor is None:
            return [dbc.Alert('Choose a base sensor.', color='warning'), dash.no_update, dash.no_update] + no_change
        try:
            recipe = [dict(row, n=float(row['n']), k=float(row['k']), **{'d guess [nm]': float(row['d guess [nm]'])}) for row in (recipe or [])]
        except (TypeError, ValueError, KeyError):
            return [dbc.Alert('Fill in n, k and the thickness guess of each layer.', color='warning'), dash.no_update, dash.no_update] + no_change

        series = SPRpy_series.parse_series_files(paths)
        current_session.add_log('Layer series of {n} samples ({steps}) with base sensor S{sensor} and layers {layers}'.format(
            n=len(series['samples']), steps=' -> '.join(series['steps']), sensor=sensor_id,
            layers=', '.join('{l} (n {n}, k {k})'.format(l=row['Layer'], n=row['n'], k=row['k']) for row in recipe)))
        try:
            results = SPRpy_series.run_layer_series(current_session, series, template_sensor, recipe, series_options_from_settings(options_values),
                                                    default_sensor_values, channel_spots=channel_spots)
        except (ValueError, KeyError, FileNotFoundError) as error:
            return [dbc.Alert('The layer series failed: {e}'.format(e=error), color='danger'), dash.no_update, dash.no_update] + no_change
        current_session.save_session()

        summary = SPRpy_series.series_summary(results)
        per_sample = results.pivot_table(index=['Step', 'Variable'], columns=['Sample', 'Channel'], values='Value', sort=False).round(2)
        per_sample.columns = ['{s} {c}'.format(s=sample, c=channel) for sample, channel in per_sample.columns]
        per_sample = per_sample.reset_index()
        flagged = results[results['Check'] != '']
        children = [dash.html.H5('Mean ± std over the samples', style={'margin-top': '15px'}),
                    dbc.Table.from_dataframe(summary, bordered=True, size='sm', striped=True, style={'font-size': '0.85rem'}),
                    dash.html.H5('All fits'),
                    dbc.Table.from_dataframe(per_sample, bordered=True, size='sm', striped=True, style={'font-size': '0.8rem'})]
        if not flagged.empty:
            children.append(dbc.Alert(['Flagged fits (at a bound, or rms above 0.01: the model does not describe the dip well, e.g. because of the '
                                       'optical constants for that wavelength in config.toml): ',
                                       ', '.join('{s} {st} {c} ({check})'.format(s=row['Sample'], st=row['Step'], c=row['Channel'], check=row['Check'])
                                                 for _, row in flagged.iterrows())], color='warning', style={'font-size': '0.85rem'}))
        children.append(dash.html.Div('{n} fresnel analyses (e.g. FM{first} {name}) and {sensors} sensors (one per sample and channel, with a layer added '
                                      'for each step) added to the session. Select an analysis to inspect its fit.'.format(
                                          n=len(results), first=int(results['FM'].iloc[0]), sensors=results['Sensor'].nunique(),
                                          name=current_session.fresnel_analysis_instances[int(results['FM'].iloc[0])].name),
                                      className='sprpy-hint'))

        fresnel_dropdown = [dbc.DropdownMenuItem('FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name,
                                                 id={'type': 'fresnel-analysis-list', 'index': fresnel_id}, n_clicks=0)
                            for fresnel_id in current_session.fresnel_analysis_instances]
        background_options = [{'label': 'FM' + str(fresnel_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_id].name, 'value': fresnel_id}
                              for fresnel_id in current_session.fresnel_analysis_instances]
        sensor_options = [dbc.DropdownMenuItem('S' + str(sensor_id) + ' ' + current_session.sensor_instances[sensor_id].name,
                                               id={'type': 'sensor-list', 'index': sensor_id}, n_clicks=0)
                          for sensor_id in current_session.sensor_instances]
        return [children, results.to_dict('records'), False, fresnel_dropdown, background_options, sensor_options, fresnel_summary_table(current_session),
                fresnel_summary_barplot(current_session)]

    @dash.callback(
        dash.Output('series-export-signal', 'data'),
        dash.Input('series-export', 'n_clicks'),
        dash.State('series-results', 'data'),
        prevent_initial_call=True)
    def export_layer_series(n_clicks, results):
        if not results:
            raise dash.exceptions.PreventUpdate
        save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                  file_types=[('CSV files', '*.csv')], default_extension='.csv')
        results = pd.DataFrame(results)
        results.to_csv(save_filename, sep=';', index=False)
        SPRpy_series.series_summary(results).to_csv(save_filename[:-4] + ' summary.csv', sep=';', index=False, encoding='utf-8-sig')
        raise dash.exceptions.PreventUpdate

    # Exclusion height plot of the result summary tab, updated whenever the exclusion height table is updated
    @dash.callback(
        dash.Output('summary-exclusion-barplot-graph', 'figure'),
        dash.Input('result-summary-exclusion-table', 'children'),
        prevent_initial_call=True)
    def update_exclusion_summary_barplot(table_children):
        return exclusion_summary_barplot(current_session)

    # Saving the plots of the result summary tab
    @dash.callback(
        dash.Output('summary-plots-save-signal', 'data'),
        dash.Input('barplot-save-png', 'n_clicks'),
        dash.Input('barplot-save-svg', 'n_clicks'),
        dash.Input('barplot-save-html', 'n_clicks'),
        dash.Input('barplot-save-csv', 'n_clicks'),
        dash.Input('exclusion-barplot-save-png', 'n_clicks'),
        dash.Input('exclusion-barplot-save-svg', 'n_clicks'),
        dash.Input('exclusion-barplot-save-html', 'n_clicks'),
        dash.Input('exclusion-barplot-save-csv', 'n_clicks'),
        prevent_initial_call=True)
    def save_summary_plots(*clicks):
        trigger = dash.ctx.triggered_id
        if not trigger or not dash.ctx.triggered[0]['value']:
            raise dash.exceptions.PreventUpdate
        exclusion = trigger.startswith('exclusion-')
        figure = exclusion_summary_barplot(current_session) if exclusion else fresnel_summary_barplot(current_session)
        file_format = trigger.rsplit('-', 1)[-1]
        if file_format in ('png', 'svg'):
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[(file_format.upper() + ' files', '*.' + file_format)], default_extension='.' + file_format)
            plotly.io.write_image(figure, save_filename, format=file_format)
        elif file_format == 'html':
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('HTML files', '*.html')], default_extension='.html')
            plotly.io.write_html(figure, save_filename, include_mathjax='cdn')
        else:
            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder,
                                      file_types=[('CSV files', '*.csv')], default_extension='.csv')
            if exclusion:
                rows = [{'Analysis': 'EH{id} {name}'.format(id=eh.object_id, name=eh.name),
                         'Exclusion height mean [nm]': eh.mean_exclusion_height_result[0],
                         'Exclusion height std [nm]': eh.mean_exclusion_height_result[1] if len(eh.mean_exclusion_height_result) > 1 else np.nan,
                         'Exclusion heights of the steps [nm]': ' '.join('{h:.3f}'.format(h=h) for h in np.asarray(eh.all_exclusion_results[0, :], dtype=np.float64))}
                        for eh in current_session.exclusion_height_analysis_instances.values()]
            else:
                rows = [{'Analysis': 'FM{id} {name}'.format(id=fm.object_id, name=fm.name), 'Sensor': 'S{id}'.format(id=fm.sensor_object.object_id),
                         'Fitted layer': fm.fitted_layer, 'Value': fm.fitted_result[0]} for fm in current_session.fresnel_analysis_instances.values()]
            pd.DataFrame(rows).to_csv(save_filename, sep=';', index=False)
        raise dash.exceptions.PreventUpdate

    @dash.callback(
        dash.Input('export-single-file-button', 'n_clicks'),
        prevent_initial_call=True)
    def export_result_summary(single_file):
        global current_session

        if 'export-single-file-button' == dash.ctx.triggered_id:
            fresnel_df = pd.DataFrame(
                                        {'Analysis': ['FM' + str(current_session.fresnel_analysis_instances[fresnel_inst].object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].name for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Sensor': ['S'+ str(current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.object_id) + ' ' + current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.name for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Variable': ['{layer}|{parameter}-{channel}'.format(
                                             layer=current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer,
                                             parameter=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.optical_parameters.columns[current_session.fresnel_analysis_instances[fresnel_inst].fitted_layer_index[1]],
                                             channel=current_session.fresnel_analysis_instances[fresnel_inst].sensor_object.channel) for fresnel_inst in current_session.fresnel_analysis_instances],
                                        'Value': [round(current_session.fresnel_analysis_instances[fresnel_inst].fitted_result[0], 3) for fresnel_inst in current_session.fresnel_analysis_instances]})

            exclusion_df = pd.DataFrame(
                                        {'Analysis': [current_session.exclusion_height_analysis_instances[exclusion_inst].name
                                                           for exclusion_inst in
                                                           current_session.exclusion_height_analysis_instances],
                                         'Exclusion height mean': ['{mean_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_height_result[0], 2)) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion height all': ['{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[0, :], decimals=2))) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion RI mean': ['{mean_}'.format(mean_=round(current_session.exclusion_height_analysis_instances[exclusion_inst].mean_exclusion_RI_result[0], 4)) for exclusion_inst in current_session.exclusion_height_analysis_instances],
                                         'Exclusion RI all': ['{all_}'.format(all_=str(np.round(current_session.exclusion_height_analysis_instances[exclusion_inst].all_exclusion_results[1, :], decimals=4))) for exclusion_inst in current_session.exclusion_height_analysis_instances]
                                         })

            save_filename = save_file(prompt='Choose save location and filename', prompt_folder=default_data_folder, file_types=[('CSV files', '*.csv')], default_extension='.csv')

            fresnel_df.to_csv(save_filename[:-4] + '_fresnel' + '.csv', sep=';')
            exclusion_df.to_csv(save_filename[:-4] + '_exclusion' + '.csv', sep=';')

    address = 'http://{host}:8050/'.format(host=session_host)
    if loading_splash is not None:
        loading_splash.status('Ready! Opening SPRpy at ' + address)
        time.sleep(0.6)
        loading_splash.close()
    print('\nSPRpy is running at ' + address + ' (open it in your browser). Close this window to quit SPRpy.\n')
    if open_browser:
        threading.Timer(1.5, webbrowser.open, args=[address]).start()
    app.run(debug=True, use_reloader=False, host=session_host, port=8050)

