# Analysis of all channels (wavelengths) of a measurement at once: fresnel modelling and exclusion height determination
#
# A .spr2 measurement is converted to one .csv file per channel (e.g. "name-L1_670nm.csv", "name-L2_980nm.csv"). All
# channels share the same scan times, so a scan or an interval of scans selected in one channel corresponds to the same
# scans in the other channels. An analysis made for one channel (the template) can therefore be repeated for the other
# channels, using the optical constants of each wavelength.

import copy
import os
import re

import numpy as np
import pandas as pd

from SPRpy_functions import load_csv_data, calculate_sensorgram
from fresnel_transfer_matrix import fresnel_calculation, angle_range_around_minimum, SPR_minimum_index, moving_mean
from SPRpy_classes import (add_fresnel_model_object, add_exclusion_height_object, run_exclusion_height_calculation)

CHANNEL_FILE_PATTERN = re.compile(r'^(?P<base>.*)-L(?P<number>\d+)_(?P<wavelength>\d{3})nm\.csv$', re.IGNORECASE)


# ----------------------------------------------------------------------------------------------------------------------
# Channel files and data
# ----------------------------------------------------------------------------------------------------------------------

def channel_label(data_path):
    """
    Channel label like 'L2 980nm' from a converted .csv file name (same as Sensor.channel).
    """
    match = CHANNEL_FILE_PATTERN.match(os.path.basename(data_path))
    if match is None:
        return os.path.basename(data_path)
    return 'L{n} {w}nm'.format(n=match['number'], w=match['wavelength'])


def strip_channel(name):
    """
    Name without a trailing channel label, e.g. 'PEG BG L1' or 'PEG BG L1 670nm' -> 'PEG BG'.
    """
    stripped = re.sub(r'\s+L\d+(\s*\d{3}\s*nm)?\s*$', '', str(name))
    return stripped if stripped else str(name)


def channel_number(data_path):
    """
    Channel (laser) number of a converted .csv file, e.g. 2 for 'name-L2_980nm.csv', or None.
    """
    match = CHANNEL_FILE_PATTERN.match(os.path.basename(data_path))
    return int(match['number']) if match is not None else None


def channel_spot(data_path, channel_spots):
    """
    Number (1, 2, ...) of the measurement spot of a channel, from channel_spots in config.toml (e.g. [[1, 2], [3, 4]]:
    L1 and L2 measure spot 1, L3 and L4 spot 2), or None if the channel is not listed.
    """
    number = channel_number(data_path)
    for spot, numbers in enumerate(channel_spots or [], start=1):
        if number is not None and number in [int(n) for n in numbers]:
            return spot
    return None


def find_channel_files(data_path, channel_spots=None):
    """
    All converted channel files (.csv) of the same measurement as data_path, sorted by channel number.

    :param channel_spots: channels measuring the same spot (config.toml), used for the 'spot' of each channel
    :return: list of dicts with 'path', 'label' (e.g. 'L2 980nm'), 'number', 'wavelength' and 'spot'. Only data_path
             itself if the file name does not follow the 'name-L<number>_<wavelength>nm.csv' pattern of the converter.
    """
    folder, file_name = os.path.split(data_path)
    match = CHANNEL_FILE_PATTERN.match(file_name)
    if match is None:
        return [{'path': data_path, 'label': file_name, 'number': 0, 'wavelength': None, 'spot': None}]

    channels = []
    for candidate in os.listdir(folder or '.'):
        candidate_match = CHANNEL_FILE_PATTERN.match(candidate)
        if candidate_match is not None and candidate_match['base'] == match['base']:
            path = data_path if candidate == file_name else os.path.join(folder, candidate).replace('\\', '/')
            channels.append({'path': path, 'label': 'L{n} {w}nm'.format(n=candidate_match['number'], w=candidate_match['wavelength']),
                             'number': int(candidate_match['number']), 'wavelength': int(candidate_match['wavelength']),
                             'spot': channel_spot(path, channel_spots)})
    return sorted(channels, key=lambda channel: channel['number'])


class MeasurementChannel:
    """
    Measurement data of one channel (one .csv file).
    """

    def __init__(self, data_path):
        self.path, self.scanspeed, self.time_df, self.angles_df, self.ydata_df, _ = load_csv_data(path=data_path)
        self.label = channel_label(data_path)

    def reflectivity(self, scan_index):
        """
        Reflectivity spectrum of one scan (index label of ydata_df, as used when hovering the sensorgram).
        """
        return pd.DataFrame(data={'angles': self.angles_df, 'ydata': self.ydata_df.loc[scan_index]})

    def sensorgram(self, SPR_TIR_fitting_parameters):
        """
        SPR and TIR sensorgram, offset to the first scan (as the sensorgram of the exclusion height tab).
        """
        sensorgram_df = calculate_sensorgram(self.time_df, self.angles_df, self.ydata_df, SPR_TIR_fitting_parameters)
        sensorgram_df['SPR angle'] = sensorgram_df['SPR angle'] - sensorgram_df['SPR angle'].iloc[0]
        sensorgram_df['TIR angle'] = sensorgram_df['TIR angle'] - sensorgram_df['TIR angle'].iloc[0]
        return sensorgram_df


def find_scan_index(ydata_df, ydata):
    """
    Index label of the scan in ydata_df that matches the spectrum ydata (the scan used for a fresnel model).
    """
    y = np.asarray(ydata, dtype=np.float64)
    differences = np.nanmean(np.abs(ydata_df.to_numpy(dtype=np.float64) - y[np.newaxis, :]), axis=1)
    return ydata_df.index[int(np.nanargmin(differences))]


# ----------------------------------------------------------------------------------------------------------------------
# Sensors for other wavelengths
# ----------------------------------------------------------------------------------------------------------------------

def default_optical_constants(default_sensor_values, layer_name, wavelength):
    """
    n and k of a layer at a wavelength from the [default_sensor_values] of config.toml, or None if not available.
    """
    names = [str(name).strip().lower() for name in default_sensor_values.get('Layer', [])]
    key = str(layer_name).strip().lower()
    if key not in names or 'n_{w}'.format(w=wavelength) not in default_sensor_values:
        return None
    index = names.index(key)
    try:
        return (float(default_sensor_values['n_{w}'.format(w=wavelength)][index]),
                float(default_sensor_values['k_{w}'.format(w=wavelength)][index]))
    except (ValueError, IndexError, KeyError):
        return None


def adapt_sensor_to_channel(template_sensor, data_path, default_sensor_values):
    """
    Copy of a sensor for another channel of the measurement. Layer thicknesses are kept, while n and k of each layer are
    taken from config.toml for the wavelength of the channel (layers not listed there, e.g. a custom layer, keep their
    values, i.e. they are assumed to have no dispersion). The bulk refractive index is kept (normally it is updated
    from the TIR angle of the channel when fitting).

    :return: sensor object (without object id), list of notes on the values that were used
    """
    sensor = copy.deepcopy(template_sensor)
    sensor.data_path = data_path
    sensor.channel = channel_label(data_path)
    match = CHANNEL_FILE_PATTERN.match(os.path.basename(data_path))
    wavelength = int(match['wavelength']) if match is not None else template_sensor.wavelength
    sensor.wavelength = wavelength

    notes = []
    table = template_sensor.optical_parameters.copy()
    for row in range(len(table)):
        layer = str(table.iloc[row, 0])
        if row == len(table) - 1 or wavelength == template_sensor.wavelength:
            continue  # Bulk, or same wavelength
        constants = default_optical_constants(default_sensor_values, layer, wavelength)
        if constants is None:
            notes.append('{layer}: n and k from {w} nm used (no values for {wl} nm in config.toml)'.format(
                layer=layer, w=template_sensor.wavelength, wl=wavelength))
            continue
        template_constants = default_optical_constants(default_sensor_values, layer, template_sensor.wavelength)
        if template_constants is not None and not np.allclose([table.iloc[row, 2], table.iloc[row, 3]], template_constants, atol=1e-6):
            notes.append('{layer}: n and k at {w} nm differ from config.toml, config.toml values used for {wl} nm'.format(
                layer=layer, w=template_sensor.wavelength, wl=wavelength))
        table.iloc[row, 2] = constants[0]
        table.iloc[row, 3] = constants[1]

    set_sensor_table(sensor, table)
    return sensor, notes


def set_sensor_table(sensor, table):
    """
    Sets the optical parameter table of a sensor and its layer arrays (in the same way as editing the sensor table).
    """
    sensor.optical_parameters = pd.DataFrame(table).reset_index(drop=True)
    sensor.layer_thicknesses = sensor.optical_parameters['d [nm]'].to_numpy()
    sensor.refractive_indices = sensor.optical_parameters['n'].to_numpy()
    sensor.extinction_coefficients = sensor.optical_parameters['k'].to_numpy()
    sensor.fitted_var = sensor.optical_parameters.iloc[sensor.fitted_layer_index]
    sensor.sensor_table_title = 'S{sensor_number} {sensor_name} - {channel} - Fit: {fitted_layer}|{fitted_param}'.format(
        sensor_number=sensor.object_id, sensor_name=sensor.name, channel=sensor.channel,
        fitted_layer=sensor.optical_parameters.iloc[sensor.fitted_layer_index[0], 0],
        fitted_param=sensor.optical_parameters.columns[sensor.fitted_layer_index[1]])


# ----------------------------------------------------------------------------------------------------------------------
# Fresnel modelling of all channels
# ----------------------------------------------------------------------------------------------------------------------

def dip_half_widths(angles, ydata, min_angle, smooth_width_deg=0.1):
    """
    SPR angle and half-widths of the SPR dip below and above it (deg), measured at half the dip depth relative to the
    highest reflectivity above min_angle (e.g. the upper limit of the TIR range). The half-width below the minimum is
    nan if the dip does not reach half depth before min_angle (e.g. at long wavelengths where the dip lies close to
    the TIR edge).
    """
    step = abs(float(np.mean(np.diff(angles))))
    smooth = moving_mean(ydata, max(3, int(round(smooth_width_deg / step))))
    minimum = SPR_minimum_index(angles, ydata, min_angle=min_angle)
    above_TIR = angles > min_angle
    level = smooth[minimum] + 0.5 * (np.nanmax(smooth[above_TIR]) - smooth[minimum])
    widths = []
    for direction in (-1, 1):
        index = minimum
        while 0 < index < len(angles) - 1 and smooth[index] < level and angles[index] > min_angle:
            index += direction
        reached = smooth[index] >= level and angles[index] > min_angle
        widths.append(abs(angles[index] - angles[minimum]) if reached else np.nan)
    return float(angles[minimum]), widths[0], widths[1]


def channel_angle_range(template_fresnel, channel_reflectivity_df, SPR_TIR_fitting_parameters):
    """
    Angle range for a channel that covers the same part of the SPR dip as the angle range of the template fresnel model.
    The dip width differs between wavelengths, so the distances of the range ends from the SPR minimum are scaled with
    the half-widths of the dips (on the same side of the minimum, or the other side if one cannot be measured), rather
    than using a fixed number of data points. Falls back to the same number of data points around the minimum if the
    dip widths cannot be determined.
    """
    min_angle = SPR_TIR_fitting_parameters['TIR range'][1]
    angles = template_fresnel.measurement_data['angles'].to_numpy(dtype=np.float64)
    ydata = template_fresnel.measurement_data['ydata'].to_numpy(dtype=np.float64)
    channel_angles = channel_reflectivity_df['angles'].to_numpy(dtype=np.float64)
    channel_ydata = channel_reflectivity_df['ydata'].to_numpy(dtype=np.float64)
    try:
        SPR_angle, below, above = dip_half_widths(angles, ydata, min_angle)
        channel_SPR_angle, channel_below, channel_above = dip_half_widths(channel_angles, channel_ydata, min_angle)
        ratio_above = channel_above / above
        ratio_below = channel_below / below if np.isfinite(channel_below / below) else ratio_above
        if not (np.isfinite(ratio_above) and np.isfinite(ratio_below) and ratio_above > 0 and ratio_below > 0):
            raise ValueError('Dip widths not determined')
        low = max(channel_SPR_angle - (SPR_angle - template_fresnel.angle_range[0]) * ratio_below, min_angle, channel_angles[0])
        high = min(channel_SPR_angle + (template_fresnel.angle_range[1] - SPR_angle) * ratio_above, channel_angles[-1])

        # At least half as many data points as the template range (narrow dips, e.g. at 980 nm, would otherwise give
        # too few points for fitting the layer, offset and prism k)
        minimum_points = int(np.count_nonzero((angles >= template_fresnel.angle_range[0]) & (angles <= template_fresnel.angle_range[1]))) // 2
        step = abs(float(np.mean(np.diff(channel_angles))))
        while np.count_nonzero((channel_angles >= low) & (channel_angles <= high)) < minimum_points:
            if low <= max(min_angle, channel_angles[0]) and high >= channel_angles[-1]:
                break
            low = max(low - step, min_angle, channel_angles[0])
            high = min(high + step, channel_angles[-1])
        return [float(low), float(high)]
    except (ValueError, IndexError, ZeroDivisionError, FloatingPointError):
        minimum = SPR_minimum_index(angles, ydata, min_angle=min_angle)
        below = int(np.count_nonzero((angles >= template_fresnel.angle_range[0]) & (angles < angles[minimum])))
        above = int(np.count_nonzero((angles > angles[minimum]) & (angles <= template_fresnel.angle_range[1])))
        return angle_range_around_minimum(channel_angles, channel_ydata, below, above, min_angle=min_angle)


def initial_offset(fresnel_object, reflectivity_df):
    """
    Initial intensity offset: measured minimum reflectivity minus the modelled reflectivity at that angle.
    """
    minimum = int(reflectivity_df['ydata'].to_numpy(dtype=np.float64).argmin())
    modelled = fresnel_calculation(angles=reflectivity_df['angles'].iloc[[minimum]].to_numpy(dtype=np.float64),
                                   wavelength=fresnel_object.sensor_object.wavelength,
                                   layer_thicknesses=np.array(fresnel_object.sensor_object.layer_thicknesses, dtype=np.float64),
                                   n_re=np.array(fresnel_object.sensor_object.refractive_indices, dtype=np.float64),
                                   n_im=np.array(fresnel_object.sensor_object.extinction_coefficients, dtype=np.float64),
                                   ydata_type='R', polarization=fresnel_object.polarization)
    return float(reflectivity_df['ydata'].iloc[minimum] - modelled[0])


def copy_fit_settings(template_fresnel, fresnel_object, reflectivity_df):
    """
    Copies the fit settings of the template fresnel model (fitted variable, bounds, offset and prism k options).
    """
    fresnel_object.polarization = template_fresnel.polarization
    fresnel_object.fit_offset = template_fresnel.fit_offset
    fresnel_object.fit_prism_k = template_fresnel.fit_prism_k
    fresnel_object.update_bulk_RI_from_TIR = getattr(template_fresnel, 'update_bulk_RI_from_TIR', True)
    fresnel_object.match_dip_angle = getattr(template_fresnel, 'match_dip_angle', False)
    fresnel_object.extinction_correction = template_fresnel.extinction_correction
    fresnel_object.bounds = channel_bounds(template_fresnel)
    fresnel_object.y_offset = initial_offset(fresnel_object, reflectivity_df) if fresnel_object.fit_offset else 0
    ini_guess = np.array(template_fresnel.ini_guess, dtype=np.float64).ravel()
    if np.isfinite(template_fresnel.fitted_result[0]):
        ini_guess[0] = template_fresnel.fitted_result[0]  # Start from the result of the template channel
    if len(ini_guess) >= 2:
        ini_guess[1] = fresnel_object.y_offset
    fresnel_object.ini_guess = ini_guess


def channel_bounds(template_fresnel):
    """
    Fit bounds for another channel. The bounds of a fitted layer thickness are made wider than in the template (by
    default they are set tightly around its initial guess), since the fitted thickness can differ between wavelengths
    (e.g. through the optical constants): up to at least 10 times the template result.
    """
    bounds = copy.deepcopy(template_fresnel.bounds)
    if template_fresnel.fitted_layer_index[1] != 1 or not np.isfinite(template_fresnel.fitted_result[0]):
        return bounds
    lower, upper = bounds[0], bounds[1]
    minimum_upper = 10 * max(float(template_fresnel.fitted_result[0]), 1.0)
    if np.ndim(lower):
        lower, upper = list(lower), list(upper)
        lower[0], upper[0] = min(lower[0], float(template_fresnel.fitted_result[0]) / 4), max(upper[0], minimum_upper)
        return [tuple(lower), tuple(upper)]
    return [min(lower, float(template_fresnel.fitted_result[0]) / 4), max(upper, minimum_upper)]


# Limits of the fit quality checks. Matching the dip angle places the model dip at the measured SPR angle at the expense
# of the rest of the curve, which roughly doubles the rms deviation of good fits (e.g. 0.002-0.003 at 670 nm and
# 0.010-0.012 at 785 nm for dry SiO2/APS/SMCC/PEG sensors, compared to 0.001-0.002 and 0.005-0.006 without matching).
RMS_LIMIT = 0.01  # Reflectivity units, reflectivity-only fit
RMS_LIMIT_DIP_MATCHED = 0.02  # Reflectivity units, fit matching the dip angle
DIP_MISMATCH_LIMIT = 0.001  # Degrees: the fit could not place the model dip at the measured angle


def fit_rms(fresnel_object):
    """
    Root mean square deviation between the fitted model and the measured reflectivity within the angle range.
    """
    angles = fresnel_object.measurement_data['angles'].to_numpy(dtype=np.float64)
    selection = (angles >= fresnel_object.angle_range[0]) & (angles <= fresnel_object.angle_range[1])
    measured = fresnel_object.measurement_data['ydata'].to_numpy(dtype=np.float64)[selection]
    return float(np.sqrt(np.mean((fresnel_object.fitted_data['ydata'].to_numpy(dtype=np.float64) - measured) ** 2)))


def fit_quality_warnings(fresnel_object, rms=None):
    """
    Warnings for a fresnel fit (all wavelengths and layer series):
    - the fitted value is at a bound (the fit did not converge within the bounds),
    - poor fit: rms deviation above RMS_LIMIT (RMS_LIMIT_DIP_MATCHED when matching the dip angle, which raises the rms),
      i.e. the model does not describe the shape of the dip, typically because the optical constants (e.g. of gold or
      the prism) for this wavelength in config.toml do not match the sensor and instrument,
    - when matching the dip angle: the fit could not place the model dip at the measured SPR angle (more than 1 mdeg off).
      Without matching, the model dip is typically several mdeg off (shown in the fit result), which is not flagged.
    Intensity offsets are not flagged: they differ between channels also for good fits (up to about 0.1 in liquid).
    """
    warnings = []
    if fit_at_bound(fresnel_object):
        warnings.append('fitted value at a bound')
    matched = bool(getattr(fresnel_object, 'match_dip_angle', False))
    if rms is None:
        try:
            rms = fit_rms(fresnel_object)
        except (AttributeError, KeyError, TypeError, ValueError):
            rms = np.nan
    if np.isfinite(rms) and rms > (RMS_LIMIT_DIP_MATCHED if matched else RMS_LIMIT):
        warnings.append('poor fit (rms {r:.3f})'.format(r=rms))
    mismatch = getattr(fresnel_object, 'dip_angle_mismatch', np.nan)
    if matched and mismatch is not None and np.isfinite(mismatch) and abs(mismatch) > DIP_MISMATCH_LIMIT:
        warnings.append('model dip {m:+.1f} mdeg from the measured SPR angle'.format(m=mismatch * 1000))
    return warnings


def channel_fit_warnings(fresnel_object, template_fresnel=None, relative_limit=0.25, absolute_limit=0.5):
    """
    Warnings for a channel fit of "Fit all wavelengths": the checks of fit_quality_warnings, and a fitted value that
    differs from the fitted value of the template (first) channel by more than relative_limit of it (and more than
    absolute_limit, e.g. nm). A channel can describe its dip well and still give a different thickness when the optical
    constants for its wavelength (gold and prism in particular) do not match the sensor and instrument, since the fit
    then compensates with the layer thickness (and typically a larger offset).
    """
    warnings = fit_quality_warnings(fresnel_object)
    if template_fresnel is not None and template_fresnel is not fresnel_object:
        value, reference = float(fresnel_object.fitted_result[0]), float(template_fresnel.fitted_result[0])
        if np.isfinite(value) and np.isfinite(reference) and abs(value - reference) > max(relative_limit * abs(reference), absolute_limit):
            channel = getattr(template_fresnel.sensor_object, 'channel', None) or 'the first channel'
            warnings.append('{d:+.2f} from {channel}'.format(d=value - reference, channel=channel))
    return warnings


def fit_at_bound(fresnel_object):
    """
    True if the fitted value lies at its lower or upper bound (the fit did not converge to a solution within the bounds).
    """
    lower, upper = fresnel_object.bounds[0], fresnel_object.bounds[1]
    lower, upper = (lower[0], upper[0]) if np.ndim(lower) else (lower, upper)
    value = float(fresnel_object.fitted_result[0])
    tolerance = 1e-4 * max(abs(upper - lower), 1e-9)
    return abs(value - lower) <= tolerance or abs(value - upper) <= tolerance


def fit_fresnel_all_channels(session, template_fresnel, default_sensor_values, channel_paths=None, channels=None, log=True):
    """
    Repeats a fitted fresnel model (template) for the other channels of the measurement, using the same scan and the
    optical constants of each wavelength. A sensor and a fresnel analysis are added to the session for each channel
    (and reused when the analysis is repeated).

    :param session: session object
    :param template_fresnel: fitted fresnel analysis object of one channel
    :param default_sensor_values: [default_sensor_values] of config.toml
    :param channel_paths: data files of the channels to include (default: all channels of the measurement)
    :param channels: optional dict {path: MeasurementChannel} with already loaded data
    :return: list of fresnel analysis objects (one per channel, the template included) and list of notes
    """
    channels = {} if channels is None else channels
    if channel_paths is None:
        channel_paths = [channel['path'] for channel in find_channel_files(template_fresnel.initial_data_path)]

    template_path = template_fresnel.initial_data_path
    if template_path not in channels:
        channels[template_path] = MeasurementChannel(template_path)
    scan_index = find_scan_index(channels[template_path].ydata_df, template_fresnel.measurement_data['ydata'])
    template_fresnel.scan_index = scan_index

    if not hasattr(template_fresnel, 'channel_fresnel_ids'):
        template_fresnel.channel_fresnel_ids = {}
    template_fresnel.channel_fresnel_ids[template_path] = template_fresnel.object_id

    fresnel_objects, notes = [], []
    for path in channel_paths:
        if path == template_path:
            fresnel_objects.append(template_fresnel)
            continue
        if path not in channels:
            channels[path] = MeasurementChannel(path)
        channel = channels[path]
        reflectivity_df = channel.reflectivity(scan_index)

        # Sensor with the optical constants of this wavelength (a new one each time, the template may have changed)
        sensor, sensor_notes = adapt_sensor_to_channel(template_fresnel.sensor_object, path, default_sensor_values)
        notes += ['{channel}: {note}'.format(channel=channel.label, note=note) for note in sensor_notes]
        existing_id = template_fresnel.channel_fresnel_ids.get(path)
        fresnel_object = session.fresnel_analysis_instances.get(existing_id) if existing_id is not None else None
        if fresnel_object is not None and fresnel_object.sensor_object.object_id in session.sensor_instances:
            sensor.object_id = fresnel_object.sensor_object.object_id
        else:
            session.sensor_ID_count += 1
            sensor.object_id = session.sensor_ID_count
        sensor.name = '{name} {channel}'.format(name=strip_channel(template_fresnel.sensor_object.name), channel=channel.label)
        set_sensor_table(sensor, sensor.optical_parameters)
        session.sensor_instances[sensor.object_id] = sensor

        name = '{name} {channel}'.format(name=strip_channel(template_fresnel.name), channel=channel.label)
        if fresnel_object is None:
            fresnel_object = add_fresnel_model_object(session, sensor, path, reflectivity_df, name)
        fresnel_object.name = name
        fresnel_object.sensor_object = sensor
        fresnel_object.initial_data_path = path
        fresnel_object.measurement_data = reflectivity_df
        fresnel_object.scan_index = scan_index
        fresnel_object.fitted_layer_index = copy.deepcopy(template_fresnel.fitted_layer_index)
        fresnel_object.fitted_layer = copy.deepcopy(template_fresnel.fitted_layer)
        fresnel_object.sensor_object_label = 'Sensor: ' + sensor.sensor_table_title
        fresnel_object.angle_range = channel_angle_range(template_fresnel, reflectivity_df, session.SPR_TIR_fitting_parameters)
        copy_fit_settings(template_fresnel, fresnel_object, reflectivity_df)
        template_fresnel.channel_fresnel_ids[path] = fresnel_object.object_id

        fresnel_object.model_reflectivity_trace()
        update_sensor_with_fit(fresnel_object)
        fresnel_objects.append(fresnel_object)

        session.save_sensor(sensor.object_id)
        session.save_fresnel_analysis(fresnel_object.object_id)
        if log:
            session.add_log('FM{id} {name}: fitted {layer}|{param} = {value:.4f}{extra} (all wavelengths from FM{template}, scan {scan}, angle range {low:.3f}-{high:.3f} deg)'.format(
                id=fresnel_object.object_id, name=fresnel_object.name, layer=fresnel_object.fitted_layer,
                param=sensor.optical_parameters.columns[fresnel_object.fitted_layer_index[1]], value=fresnel_object.fitted_result[0],
                extra=fit_extra_text(fresnel_object), template=template_fresnel.object_id, scan=scan_index,
                low=fresnel_object.angle_range[0], high=fresnel_object.angle_range[1]))

    session.save_fresnel_analysis(template_fresnel.object_id)
    return fresnel_objects, notes


def fit_extra_text(fresnel_object):
    extra = ''
    if fresnel_object.fit_offset and len(fresnel_object.fitted_result) > 1:
        extra += ', offset {o:.4f}'.format(o=fresnel_object.fitted_result[1])
    if fresnel_object.fit_prism_k and len(fresnel_object.fitted_result) > 2:
        extra += ', prism k {k:.5f}'.format(k=fresnel_object.fitted_result[2])
    return extra


def update_sensor_with_fit(fresnel_object):
    """
    Writes the fitted value (and fitted prism k) to the sensor table, as after fitting in the GUI.
    """
    sensor = fresnel_object.sensor_object
    sensor.optical_parameters.iloc[fresnel_object.fitted_layer_index] = round(fresnel_object.fitted_result[0], 4)
    if fresnel_object.fit_prism_k and len(fresnel_object.fitted_result) > 2:
        sensor.optical_parameters.iloc[(0, 3)] = round(fresnel_object.fitted_result[2], 5)
    else:
        sensor.optical_parameters.iloc[(0, 3)] = sensor.extinction_coefficients[0]
    # Keep the layer arrays in step with the table (they are views of it with pandas 2, but do not rely on that)
    sensor.layer_thicknesses = sensor.optical_parameters['d [nm]'].to_numpy()
    sensor.refractive_indices = sensor.optical_parameters['n'].to_numpy()
    sensor.extinction_coefficients = sensor.optical_parameters['k'].to_numpy()


# ----------------------------------------------------------------------------------------------------------------------
# Exclusion height of all channels
# ----------------------------------------------------------------------------------------------------------------------

def exclusion_height_all_channels(session, template_exclusion, logical_cores, channel_paths=None, channels=None,
                                  bulk_RI_method=None, progress=None, log=True):
    """
    Runs the exclusion height analysis (template) for all channels of the measurement, with the same injection, buffer
    and probe points and settings. The fresnel background of each channel is the corresponding channel analysis of the
    template's fresnel background (see fit_fresnel_all_channels, which must be run first). An exclusion height analysis
    is added to the session for each channel (and reused when the analysis is repeated).

    :param progress: optional function called with a status message
    :return: list of exclusion height analysis objects (one per channel, the template included)
    """
    channels = {} if channels is None else channels
    template_fresnel = template_exclusion.fresnel_object
    channel_fresnel_ids = getattr(template_fresnel, 'channel_fresnel_ids', {})
    if channel_paths is None:
        channel_paths = [channel['path'] for channel in find_channel_files(template_exclusion.initial_data_path)]
    template_path = template_exclusion.initial_data_path
    missing = [channel_label(path) for path in channel_paths
               if path != template_path and channel_fresnel_ids.get(path) not in session.fresnel_analysis_instances]
    if missing:
        raise ValueError('First fit the fresnel background FM{id} for all wavelengths (missing: {missing}).'.format(
            id=template_fresnel.object_id, missing=', '.join(missing)))
    if bulk_RI_method is None:
        bulk_RI_method = getattr(template_exclusion, 'bulk_RI_method', 'edge_overlay')

    if not hasattr(template_exclusion, 'channel_exclusion_ids'):
        template_exclusion.channel_exclusion_ids = {}
    template_exclusion.channel_exclusion_ids[template_path] = template_exclusion.object_id

    exclusion_objects = []
    for path in channel_paths:
        if path not in channels:
            channels[path] = MeasurementChannel(path)
        channel = channels[path]
        if progress is not None:
            progress('Exclusion height {channel}...'.format(channel=channel.label))

        if path == template_path:
            exclusion_object = template_exclusion
        else:
            fresnel_object = session.fresnel_analysis_instances[channel_fresnel_ids[path]]
            sensorgram_df = channel.sensorgram(session.SPR_TIR_fitting_parameters)
            existing_id = template_exclusion.channel_exclusion_ids.get(path)
            exclusion_object = session.exclusion_height_analysis_instances.get(existing_id) if existing_id is not None else None
            name = '{name} {channel}'.format(name=strip_channel(template_exclusion.name), channel=channel.label)
            if exclusion_object is None:
                exclusion_object = add_exclusion_height_object(session, fresnel_object, sensorgram_df, path, name)
            exclusion_object.name = name
            exclusion_object.fresnel_object = fresnel_object
            exclusion_object.fresnel_object_label = 'Fresnel background: FM{id} {name}'.format(id=fresnel_object.object_id, name=fresnel_object.name)
            exclusion_object.sensor_object = fresnel_object.sensor_object
            exclusion_object.polarization = fresnel_object.polarization
            exclusion_object.initial_data_path = path
            exclusion_object.sensorgram_data = sensorgram_df.iloc[:, :3]
            for attribute in ('sensorgram_offset_ind', 'd_n_pair_resolution', 'RI_bounds', 'RI_initial_guess',
                              'fit_offset', 'fit_prism', 'height_search', 'match_dip_angle', 'injection_points', 'buffer_points', 'probe_points'):
                setattr(exclusion_object, attribute, copy.deepcopy(getattr(template_exclusion, attribute, {'height_search': 'grid', 'match_dip_angle': False}.get(attribute))))
            exclusion_object.height_bounds = list(template_exclusion.height_bounds)
            exclusion_object.height_steps = np.array(template_exclusion.height_steps, dtype=np.float64)
            template_exclusion.channel_exclusion_ids[path] = exclusion_object.object_id

        exclusion_object.initialize_model(channel.ydata_df, bulk_RI_method=bulk_RI_method)
        run_exclusion_height_calculation(exclusion_object, logical_cores)
        session.save_exclusion_height_analysis(exclusion_object.object_id)
        exclusion_objects.append(exclusion_object)
        if log:
            session.add_log('EH{id} {name}: {height:.2f} nm (std {std:.2f}), heights {heights} nm, TIR shifts {shifts} mdeg (all wavelengths from EH{template})'.format(
                id=exclusion_object.object_id, name=exclusion_object.name, height=exclusion_object.mean_exclusion_height_result[0],
                std=exclusion_object.mean_exclusion_height_result[1], heights=np.round(exclusion_object.all_exclusion_results[0, :], 2),
                shifts=np.round(np.asarray(exclusion_object.TIR_shifts, dtype=np.float64) * 1000, 1), template=template_exclusion.object_id))

    return exclusion_objects


def results_by_spot(analyses, channel_spots, value):
    """
    Results of the channel analyses grouped by measurement spot. Channels of the same spot measure the same part of the
    sample (and should agree, apart from the different sensitivity of each wavelength), while differences between spots
    can indicate an inhomogeneous sample. Results are therefore compared per spot rather than averaged over all channels.

    :param analyses: fresnel or exclusion height analysis objects (one per channel)
    :param channel_spots: channels measuring the same spot (config.toml)
    :param value: function giving the result of an analysis
    :return: list of (spot number or None, [(channel label, wavelength, result), ...]), sorted by spot
    """
    groups = {}
    for analysis in analyses:
        spot = channel_spot(analysis.initial_data_path, channel_spots)
        groups.setdefault(spot, []).append((analysis.sensor_object.channel, analysis.sensor_object.wavelength, value(analysis)))
    return sorted(groups.items(), key=lambda item: (item[0] is None, item[0] or 0))


def spot_summary_lines(analyses, channel_spots, value, unit='', decimals=2):
    """
    One line per measurement spot comparing the results of its channels, e.g. 'Spot 1 (L1 670nm, L2 980nm): 6.84, 7.10 nm'.
    """
    lines = []
    for spot, channels in results_by_spot(analyses, channel_spots, value):
        name = 'Spot {spot}'.format(spot=spot) if spot is not None else 'Channels not in channel_spots'
        values = ', '.join('nan' if not np.isfinite(result) else '{v:.{d}f}'.format(v=result, d=decimals) for _, _, result in channels)
        lines.append('{name} ({channels}): {values}{unit}'.format(name=name, channels=', '.join(label for label, _, _ in channels),
                                                                  values=values, unit=' ' + unit if unit else ''))
    return lines
