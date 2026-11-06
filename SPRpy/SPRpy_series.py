# Layer series: fresnel modelling of repeated samples measured after each added layer (e.g. clean SiO2 -> APS -> SMCC ->
# PEG), for all channels (wavelengths) at once.
#
# The sample and the step of each measurement are taken from the file names of the converted .csv files:
#     <date>-<time>_<sample>_<step>-L<channel>_<wavelength>nm.csv   (the date-time prefix is optional)
# e.g. 20260116-135311_SiO2_1_APS-L3_670nm.csv is sample "SiO2_1", step "APS", channel L3 at 670 nm. The steps are
# ordered by measurement time. For each sample and channel, the first step is fitted with a base sensor (the fitted
# variable of a template sensor, e.g. the SiO2 thickness), and each following step adds one layer whose thickness is
# fitted with the earlier layers fixed at the results of that sample and channel.

import copy
import os
import re

import numpy as np
import pandas as pd

from fresnel_transfer_matrix import select_TIR_range, angle_range_around_minimum, fresnel_calculation
from SPRpy_classes import add_fresnel_model_object
import SPRpy_multichannel as M

SERIES_FILE_PATTERN = re.compile(r'^(?:(?P<time>\d{8}-\d{6})_)?(?P<sample>.+)_(?P<step>[^_]+)-L(?P<number>\d+)_(?P<wavelength>\d{3})nm\.csv$', re.IGNORECASE)


# ----------------------------------------------------------------------------------------------------------------------
# Files
# ----------------------------------------------------------------------------------------------------------------------

def parse_series_files(paths):
    """
    Groups measurement files by sample, step and channel from their file names.

    :param paths: .csv files of the converter (one per channel)
    :return: dict with 'samples' {sample: {step: {channel number: path}}}, 'steps' (all steps, ordered by the earliest
             measurement time, or file modification time if the names have no date-time prefix), 'wavelengths'
             {channel number: wavelength} and 'skipped' (file names that do not match the pattern)
    """
    samples, times, wavelengths, skipped = {}, {}, {}, []
    for path in paths:
        path = str(path).replace('\\', '/')
        match = SERIES_FILE_PATTERN.match(os.path.basename(path))
        if match is None:
            skipped.append(os.path.basename(path))
            continue
        number = int(match['number'])
        samples.setdefault(match['sample'], {}).setdefault(match['step'], {})[number] = path
        wavelengths[number] = int(match['wavelength'])
        time = match['time'] if match['time'] else '{t:020.6f}'.format(t=os.path.getmtime(path) if os.path.exists(path) else 0)
        times[match['step']] = min(times.get(match['step'], time), time)
    steps = sorted(times, key=lambda step: times[step])
    return {'samples': {sample: samples[sample] for sample in sorted(samples)}, 'steps': steps,
            'wavelengths': dict(sorted(wavelengths.items())), 'skipped': skipped}


def default_recipe(series, default_sensor_values):
    """
    Recipe rows for the steps after the first (base) step: the layer added in each step (named after the step), with n
    and k from config.toml if a layer of that name is listed there (else n = 1.45, k = 0) and an initial thickness guess.
    """
    rows = []
    for step in series['steps'][1:]:
        wavelength = next(iter(series['wavelengths'].values()), 670)
        constants = M.default_optical_constants(default_sensor_values, step, wavelength)
        rows.append({'Step': step, 'Layer': step, 'n': constants[0] if constants else 1.45, 'k': constants[1] if constants else 0.0, 'd guess [nm]': 2.0})
    return rows


def mean_spectrum(path):
    """
    Reflectivity spectrum averaged over all scans of a measurement (dry measurements usually contain one or a few scans).
    """
    channel = M.MeasurementChannel(path)
    return pd.DataFrame({'angles': channel.angles_df.to_numpy(dtype=np.float64),
                         'ydata': channel.ydata_df.to_numpy(dtype=np.float64).mean(axis=0)})


# ----------------------------------------------------------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------------------------------------------------------

def fit_series_step(session, sensor, path, name, options, reflectivity_df=None):
    """
    Adds a fresnel analysis of one measurement to the session and fits the fitted variable of the sensor.

    :param options: dict with 'fit_offset', 'fit_prism_k', 'bulk_from_TIR', 'angle_range_points' ([below, above]) and
                    'TIR_ranges' (TIR ranges to choose from, e.g. water and air)
    :return: the fresnel analysis object
    """
    reflectivity_df = mean_spectrum(path) if reflectivity_df is None else reflectivity_df
    parameters = dict(session.SPR_TIR_fitting_parameters)
    TIR_range = select_TIR_range(reflectivity_df['angles'], reflectivity_df['ydata'], options['TIR_ranges'])
    if TIR_range is not None:
        parameters['TIR range'] = TIR_range

    fm = add_fresnel_model_object(session, sensor, path, reflectivity_df, name)
    fm.SPR_TIR_fitting_parameters = parameters  # TIR range of this measurement (air or liquid), without changing the session
    fm.sensor_object_label = 'Sensor: ' + sensor.sensor_table_title
    fm.angle_range = angle_range_around_minimum(reflectivity_df['angles'], reflectivity_df['ydata'], *options['angle_range_points'],
                                                min_angle=parameters['TIR range'][1])
    fm.fit_offset = options['fit_offset']
    fm.fit_prism_k = options['fit_prism_k'] and options['fit_offset']
    fm.update_bulk_RI_from_TIR = options['bulk_from_TIR']
    fm.extinction_correction = 0
    fm.y_offset = M.initial_offset(fm, reflectivity_df) if fm.fit_offset else 0

    guess = float(sensor.optical_parameters.iloc[sensor.fitted_layer_index])
    if sensor.fitted_layer_index[1] == 1:
        low, high = 0.01, max(200.0, 10 * guess)  # Thickness
    else:
        low, high = guess / 4 if guess > 0 else 0.0, max(guess * 4, guess + 1)  # n or k
    guess = min(max(guess, low * 1.01), high)
    if fm.fit_offset and fm.fit_prism_k:
        fm.ini_guess = np.array([guess, fm.y_offset, 0.001])
        fm.bounds = [(low, -np.inf, 0), (high, np.inf, 0.1)]
    elif fm.fit_offset:
        fm.ini_guess = np.array([guess, fm.y_offset])
        fm.bounds = [(low, -np.inf), (high, np.inf)]
    else:
        fm.ini_guess = np.array([guess])
        fm.bounds = [low, high]

    fm.model_reflectivity_trace()
    M.update_sensor_with_fit(fm)
    return fm


def fit_rms(fm):
    """
    Root mean square deviation between the fitted model and the measured reflectivity within the angle range.
    """
    angles = fm.measurement_data['angles'].to_numpy(dtype=np.float64)
    selection = (angles >= fm.angle_range[0]) & (angles <= fm.angle_range[1])
    return float(np.sqrt(np.mean((fm.fitted_data['ydata'].to_numpy(dtype=np.float64) - fm.measurement_data['ydata'].to_numpy(dtype=np.float64)[selection]) ** 2)))


def add_layer(sensor, layer, n, k, thickness):
    """
    Adds a layer on top of the sensor (below the bulk medium) and makes its thickness the fitted variable.
    """
    table = sensor.optical_parameters.copy()
    table = pd.concat([table.iloc[:-1], pd.DataFrame({'Layers': [layer], 'd [nm]': [float(thickness)], 'n': [float(n)], 'k': [float(k)]}),
                       table.iloc[-1:]]).reset_index(drop=True)
    sensor.fitted_layer_index = (len(table) - 2, 1)
    M.set_sensor_table(sensor, table)


def register_sensor(session, sensor, name):
    session.sensor_ID_count += 1
    sensor.object_id = session.sensor_ID_count
    sensor.name = name
    M.set_sensor_table(sensor, sensor.optical_parameters)
    session.sensor_instances[sensor.object_id] = sensor
    return sensor


def run_layer_series(session, series, template_sensor, recipe, options, default_sensor_values, channel_spots=None, progress=None, log=True):
    """
    Fits all samples and channels of a layer series (see the top of this file). The first step of each sample and channel
    is fitted with a copy of the template sensor (optical constants of the channel's wavelength, fitted variable of the
    template), and each following step adds the layer of the recipe and fits its thickness, with the earlier layers fixed
    at the results of the same sample and channel. One sensor is added to the session for each sample and channel (each
    step adds its layer to it, as when adding layers manually), and a fresnel analysis for each fit. Note that the
    analyses of earlier steps then refer to the sensor with all layers: their fitted results and plots are kept, but
    running them again would use the full layer stack.

    :param series: result of parse_series_files
    :param recipe: rows (dicts) with 'Step', 'Layer', 'n', 'k' and 'd guess [nm]' for the steps after the first. Layers
                   listed in config.toml get the n and k of each wavelength from there, other layers the recipe values.
    :param options: see fit_series_step
    :param progress: optional function called with a status message
    :return: data frame with one row per fit
    """
    steps = series['steps']
    recipe_by_step = {str(row['Step']): row for row in recipe}
    missing = [step for step in steps[1:] if step not in recipe_by_step]
    if missing:
        raise ValueError('No layer given for step(s): ' + ', '.join(missing))

    rows = []
    for sample, sample_steps in series['samples'].items():
        channels = sorted({number for step_files in sample_steps.values() for number in step_files})
        for number in channels:
            sensor = None
            for step_index, step in enumerate(steps):
                path = sample_steps.get(step, {}).get(number)
                if path is None:
                    if step_index == 0:
                        break  # No base measurement for this channel
                    continue  # Step not measured for this sample/channel: the next step builds on the previous one
                label = M.channel_label(path)
                if progress is not None:
                    progress('{sample} {step} {label}'.format(sample=sample, step=step, label=label))

                if sensor is None:
                    # One sensor per sample and channel (optical constants of the channel's wavelength)
                    sensor, _ = M.adapt_sensor_to_channel(template_sensor, path, default_sensor_values)
                    sensor.fitted_layer_index = copy.deepcopy(template_sensor.fitted_layer_index)
                    register_sensor(session, sensor, '{sample} {label}'.format(sample=sample, label=label))
                else:
                    # The layer of this step is added to the same sensor, on top of the fitted earlier layers
                    sensor.data_path = path
                    row = recipe_by_step[step]
                    n, k = layer_constants(default_sensor_values, row, sensor.wavelength)
                    add_layer(sensor, row['Layer'], n, k, float(row['d guess [nm]']))

                fm = fit_series_step(session, sensor, path, '{sample} {step} {label}'.format(sample=sample, step=step, label=label), options)
                session.save_sensor(sensor.object_id)
                session.save_fresnel_analysis(fm.object_id)
                layer = str(sensor.optical_parameters.iloc[sensor.fitted_layer_index[0], 0])
                parameter = sensor.optical_parameters.columns[sensor.fitted_layer_index[1]]
                rms = fit_rms(fm)
                warnings = series_fit_warnings(fm, rms, options.get('rms_limit', 0.01))
                rows.append({'Sample': sample, 'Step': step, 'Layer': layer, 'Variable': '{layer}|{parameter}'.format(layer=layer, parameter=parameter),
                             'Channel': label, 'Spot': M.channel_spot(path, channel_spots), 'Wavelength': sensor.wavelength,
                             'Value': float(fm.fitted_result[0]),
                             'Offset': float(fm.fitted_result[1]) if fm.fit_offset and len(fm.fitted_result) > 1 else np.nan,
                             'Prism k': float(fm.fitted_result[2]) if fm.fit_prism_k and len(fm.fitted_result) > 2 else np.nan,
                             'Bulk RI': float(sensor.refractive_indices[-1]), 'n': float(sensor.refractive_indices[sensor.fitted_layer_index[0]]),
                             'RMS': rms, 'Check': ', '.join(warnings), 'FM': fm.object_id, 'Sensor': sensor.object_id,
                             'File': os.path.basename(path)})
                if log:
                    session.add_log('Layer series: FM{id} {name}: fitted {variable} = {value:.4f} (offset {offset:.4f}, rms {rms:.4f}{check})'.format(
                        id=fm.object_id, name=fm.name, variable=rows[-1]['Variable'], value=rows[-1]['Value'],
                        offset=rows[-1]['Offset'] if np.isfinite(rows[-1]['Offset']) else 0, rms=rms,
                        check=', ' + rows[-1]['Check'] if rows[-1]['Check'] else ''))
    return pd.DataFrame(rows)


def series_fit_warnings(fm, rms, rms_limit=0.01):
    """
    Warnings for a fit of the series: fitted value at a bound, or a poor fit (root mean square deviation above
    rms_limit, in reflectivity units). Intensity offsets are not flagged, since they differ between channels of dry
    measurements also for good fits.
    """
    warnings = []
    if M.fit_at_bound(fm):
        warnings.append('fitted value at a bound')
    if rms > rms_limit:
        warnings.append('poor fit (rms {r:.3f})'.format(r=rms))
    return warnings


def layer_constants(default_sensor_values, recipe_row, wavelength):
    """
    n and k of the layer of a recipe row at a wavelength: from config.toml if the layer is listed there, else the values
    of the recipe (no dispersion).
    """
    constants = M.default_optical_constants(default_sensor_values, recipe_row['Layer'], wavelength)
    return constants if constants is not None else (float(recipe_row['n']), float(recipe_row['k']))


# ----------------------------------------------------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------------------------------------------------

def series_summary(results):
    """
    Summary of a layer series: for each step (layer), the mean and standard deviation over the samples for each channel,
    and over all samples and channels.

    :return: data frame with one row per step
    """
    if results.empty:
        return pd.DataFrame()
    summary = []
    for (step, variable), group in results.groupby(['Step', 'Variable'], sort=False):
        row = {'Step': step, 'Variable': variable}
        for channel, channel_group in group.groupby('Channel', sort=False):
            values = channel_group['Value'].to_numpy(dtype=np.float64)
            row[channel] = '{m:.2f} ± {s:.2f}'.format(m=np.mean(values), s=np.std(values)) if len(values) > 1 else '{m:.2f}'.format(m=values[0])
        values = group['Value'].to_numpy(dtype=np.float64)
        row['All'] = '{m:.2f} ± {s:.2f}'.format(m=np.mean(values), s=np.std(values))
        row['Flagged fits'] = int((group['Check'] != '').sum())
        summary.append(row)
    return pd.DataFrame(summary)
