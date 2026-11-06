# Analysis of all channels (wavelengths) of a measurement: simulated two-channel measurement with a known layer
import copy
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from fresnel_transfer_matrix import fresnel_calculation
from SPRpy_functions import load_config
from SPRpy_classes import Sensor, FresnelModel
import SPRpy_multichannel as M

LAYER_D = 4.0  # nm, true thickness of the extra layer
BUFFER_N, PROBE_N = 1.3350, 1.3362
BRUSH = (8.0, 1.38)
WAVELENGTHS = {1: 670, 4: 785}
TIR_PARAMETERS = {'TIR range': [60.8, 63.0], 'TIR method': 'savitzky_golay', 'TIR smooth width deg': 0.15,
                  'TIR fit halfwidth deg': 0.045, 'SPR method': 'dip_level', 'SPR dip level': 0.2, 'TIR fit points': 200,
                  'SPR fit points': 200, 'Fresnel_angle_range_points': [40, 60], 'sensorgram_angle_range_points': [70, 70]}


@pytest.fixture(scope='module')
def defaults():
    return load_config()['default_sensor_values']


def config_stack(defaults, wavelength, layer=None, bulk_n=BUFFER_N):
    """
    d, n, k of a Cr/Au sensor (optical constants from config.toml) with an optional extra layer (thickness, n).
    """
    n = [float(defaults['n_{w}'.format(w=wavelength)][i]) for i in (0, 1, 2)]
    k = [float(defaults['k_{w}'.format(w=wavelength)][i]) for i in (0, 1, 2)]
    d = [np.nan, float(defaults['d_nm'][1]), float(defaults['d_nm'][2])]
    if layer is not None:
        d.append(layer[0]); n.append(layer[1]); k.append(0.0)
    return np.array(d + [np.nan]), np.array(n + [bulk_n]), np.array(k + [0.0])


def write_measurement(folder, defaults, scans):
    """
    Writes one .csv file per channel (converter format). scans: list of (layer, bulk_n) for each scan.
    """
    angles = np.linspace(58.0, 78.0, 2731)  # Same angle step as a measurement (~0.0073 deg)
    paths = {}
    for number, wavelength in WAVELENGTHS.items():
        rows = []
        for index, (layer, bulk_n) in enumerate(scans):
            d, n, k = config_stack(defaults, wavelength, layer=layer, bulk_n=bulk_n)
            rows.append([0.1 * index] + list(fresnel_calculation(angles=angles, wavelength=wavelength, layer_thicknesses=d, n_re=n, n_im=k)))
        path = os.path.join(folder, 'Sample-L{n}_{w}nm.csv'.format(n=number, w=wavelength)).replace('\\', '/')
        with open(path, 'w') as file:
            file.write('# Left most column is Time (min), First row is Angles (deg), Scanspeed=5\n')
            file.write(';'.join(['0'] + ['{a:.6f}'.format(a=a) for a in angles]) + '\n')
            for row in rows:
                file.write(';'.join('{v:.7f}'.format(v=v) for v in row) + '\n')
        paths[wavelength] = path
    return paths


def make_session():
    session = SimpleNamespace(sensor_instances={}, fresnel_analysis_instances={}, exclusion_height_analysis_instances={},
                              sensor_ID_count=0, fresnel_analysis_ID_count=0, exclusion_height_analysis_ID_count=0,
                              SPR_TIR_fitting_parameters=dict(TIR_PARAMETERS), log=[])
    session.save_sensor = session.save_fresnel_analysis = session.save_exclusion_height_analysis = lambda *args: None
    session.add_log = lambda message, save=False: session.log.append(message)
    return session


def template_fresnel(session, defaults, path, layer=('PEG', 2.0, 1.456), scan_index=1, bulk_RI_from_TIR=True):
    """
    Fitted fresnel model of the 670 nm channel, fitting the thickness of an extra layer (as made in the GUI).
    """
    sensor = Sensor(path, 1, defaults)
    table = sensor.optical_parameters
    table = pd.concat([table.iloc[:-1], pd.DataFrame({'Layers': [layer[0]], 'd [nm]': [layer[1]], 'n': [layer[2]], 'k': [0.0]}), table.iloc[-1:]]).reset_index(drop=True)
    table.iloc[-1, 2] = BUFFER_N
    sensor.fitted_layer_index = (3, 1)
    M.set_sensor_table(sensor, table)
    session.sensor_ID_count = 1
    session.sensor_instances[1] = sensor

    channel = M.MeasurementChannel(path)
    reflectivity_df = channel.reflectivity(scan_index)
    fm = FresnelModel(session, sensor, path, reflectivity_df, 1, 'Layer')
    session.fresnel_analysis_ID_count = 1
    session.fresnel_analysis_instances[1] = fm
    fm.angle_range = M.angle_range_around_minimum(reflectivity_df['angles'], reflectivity_df['ydata'], 60, 60, min_angle=63.0)
    fm.fit_offset, fm.fit_prism_k = True, True
    fm.update_bulk_RI_from_TIR = bulk_RI_from_TIR
    fm.y_offset = 0.0
    fm.ini_guess = np.array([layer[1], 0.0, 0.001])
    fm.bounds = [(0.1, -np.inf, 0), (20, np.inf, 0.1)]
    fm.model_reflectivity_trace()
    M.update_sensor_with_fit(fm)
    return fm


def test_channel_files_and_labels(tmp_path):
    for name in ('Sample-L1_670nm.csv', 'Sample-L2_980nm.csv', 'Sample-L10_785nm.csv', 'Other-L1_670nm.csv', 'Sample.spr2'):
        (tmp_path / name).write_text('')
    path = str(tmp_path / 'Sample-L1_670nm.csv').replace('\\', '/')
    channels = M.find_channel_files(path)
    assert [c['label'] for c in channels] == ['L1 670nm', 'L2 980nm', 'L10 785nm']
    assert channels[0]['path'] == path
    assert [c['wavelength'] for c in channels] == [670, 980, 785]
    single = str(tmp_path / 'Sample.spr2')
    assert M.find_channel_files(single)[0]['path'] == single
    spots = [[1, 2], [3, 4]]
    assert [c['spot'] for c in M.find_channel_files(path, spots)] == [1, 1, None]  # L10 is not in channel_spots
    assert M.channel_spot(str(tmp_path / 'Sample-L4_785nm.csv'), spots) == 2
    assert M.channel_spot(path, None) is None
    assert M.strip_channel('PEG BG L1') == 'PEG BG'
    assert M.strip_channel('PEG BG L1 670nm') == 'PEG BG'
    assert M.strip_channel('Lipid layer') == 'Lipid layer'


def test_sensor_for_other_wavelength_uses_config_values(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [(None, BUFFER_N)])
    sensor = Sensor(paths[670], 1, defaults)
    table = sensor.optical_parameters
    table = pd.concat([table.iloc[:-1], pd.DataFrame({'Layers': ['PEG', 'Custom'], 'd [nm]': [2.5, 3.0], 'n': [1.456, 1.40], 'k': [0.0, 0.0]}), table.iloc[-1:]]).reset_index(drop=True)
    table.iloc[2, 1] = 48.0  # Fitted gold thickness is kept
    M.set_sensor_table(sensor, table)

    adapted, notes = M.adapt_sensor_to_channel(sensor, paths[785], defaults)
    assert adapted.wavelength == 785 and adapted.channel == 'L4 785nm'
    assert adapted.refractive_indices[2] == pytest.approx(float(defaults['n_785'][2]))
    assert adapted.extinction_coefficients[2] == pytest.approx(float(defaults['k_785'][2]))
    assert adapted.refractive_indices[0] == pytest.approx(float(defaults['n_785'][0]))
    assert adapted.layer_thicknesses[2] == 48.0 and adapted.layer_thicknesses[4] == 3.0
    assert adapted.refractive_indices[4] == 1.40  # Custom layer kept
    assert any('Custom' in note for note in notes)
    assert sensor.wavelength == 670 and sensor.refractive_indices[2] == pytest.approx(float(defaults['n_670'][2]))  # Template unchanged


def test_fresnel_all_channels_finds_layer_at_each_wavelength(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [(('PEG', LAYER_D, 1.456)[1:], BUFFER_N)] * 3)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670], bulk_RI_from_TIR=False)  # Exact bulk RI, so that the true thickness is found
    assert fm.fitted_result[0] == pytest.approx(LAYER_D, abs=0.05)

    fresnel_objects, notes = M.fit_fresnel_all_channels(session, fm, defaults)
    assert len(fresnel_objects) == 2 and fresnel_objects[0] is fm
    other = fresnel_objects[1]
    assert other.sensor_object.wavelength == 785
    assert other.fitted_result[0] == pytest.approx(LAYER_D, abs=0.05)
    assert other.sensor_object.optical_parameters.iloc[3, 1] == pytest.approx(LAYER_D, abs=0.05)
    assert other.sensor_object.layer_thicknesses[3] == pytest.approx(LAYER_D, abs=0.05)
    assert other.update_bulk_RI_from_TIR is False and other.sensor_object.refractive_indices[-1] == BUFFER_N
    assert other.scan_index == fm.scan_index == 1
    assert len(session.fresnel_analysis_instances) == 2 and len(session.sensor_instances) == 2

    # Fitting again reuses the analysis and sensor of each channel
    M.fit_fresnel_all_channels(session, fm, defaults)
    assert len(session.fresnel_analysis_instances) == 2 and len(session.sensor_instances) == 2


def test_fresnel_all_channels_with_bulk_RI_from_TIR(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [(('PEG', LAYER_D, 1.456)[1:], BUFFER_N)] * 3)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670])
    fresnel_objects, _ = M.fit_fresnel_all_channels(session, fm, defaults)
    other = fresnel_objects[1]
    assert other.sensor_object.refractive_indices[-1] == pytest.approx(BUFFER_N, abs=5e-4)  # From the TIR angle of the channel
    assert other.fitted_result[0] == pytest.approx(fm.fitted_result[0], abs=0.5)  # Same small bias from the TIR-based bulk RI


def test_exclusion_height_needs_fresnel_background_of_all_channels(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [(BRUSH, BUFFER_N)] * 3)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670], layer=('Brush', BRUSH[0], BRUSH[1]))
    eh = M.add_exclusion_height_object(session, fm, M.MeasurementChannel(paths[670]).sensorgram(TIR_PARAMETERS), paths[670], 'Brush')
    with pytest.raises(ValueError, match='fresnel background'):
        M.exclusion_height_all_channels(session, eh, 1)


@pytest.mark.slow
def test_exclusion_height_all_channels_finds_brush_height(tmp_path, defaults):
    scans = [(BRUSH, BUFFER_N)] * 6 + [(BRUSH, PROBE_N)] * 6 + [(BRUSH, BUFFER_N)] * 6
    paths = write_measurement(str(tmp_path), defaults, scans)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670], layer=('Brush', BRUSH[0], BRUSH[1]))
    M.fit_fresnel_all_channels(session, fm, defaults)

    eh = M.add_exclusion_height_object(session, fm, M.MeasurementChannel(paths[670]).sensorgram(TIR_PARAMETERS), paths[670], 'Brush')
    eh.injection_points = [(3, 0.3, 0.0), (15, 1.5, 0.0)]
    eh.buffer_points = [(0, 0.0, 0.0), (6, 0.6, 0.0), (12, 1.2, 0.0), (18, 1.8, 0.0)]
    eh.probe_points = [(6, 0.6, 0.0), (12, 1.2, 0.0)]
    eh.height_bounds = [3, 15]
    eh.height_steps = np.linspace(3, 15, 25)
    eh.d_n_pair_resolution = 25

    exclusion_objects = M.exclusion_height_all_channels(session, eh, 2)
    assert [e.sensor_object.wavelength for e in exclusion_objects] == [670, 785]
    for exclusion_object in exclusion_objects:
        assert exclusion_object.mean_exclusion_height_result[0] == pytest.approx(BRUSH[0], abs=1.0)
        assert exclusion_object.injection_points == eh.injection_points
    lines = M.spot_summary_lines(exclusion_objects, [[1, 2], [3, 4]], lambda e: float(e.mean_exclusion_height_result[0]), unit='nm')
    assert lines[0].startswith('Spot 1 (L1 670nm): ') and lines[1].startswith('Spot 2 (L4 785nm): ') and lines[1].endswith(' nm')
    assert len(session.exclusion_height_analysis_instances) == 2


def test_results_are_grouped_by_spot():
    def analysis(number, wavelength, result):
        return SimpleNamespace(initial_data_path='C:/data/Sample-L{n}_{w}nm.csv'.format(n=number, w=wavelength), result=result,
                               sensor_object=SimpleNamespace(channel='L{n} {w}nm'.format(n=number, w=wavelength), wavelength=wavelength))
    analyses = [analysis(4, 785, 9.94), analysis(1, 670, 6.84), analysis(3, 670, 12.28), analysis(2, 980, np.nan), analysis(5, 850, 1.0)]
    groups = M.results_by_spot(analyses, [[1, 2], [3, 4]], lambda a: a.result)
    assert [spot for spot, _ in groups] == [1, 2, None]
    assert [label for label, _, _ in groups[1][1]] == ['L4 785nm', 'L3 670nm']
    lines = M.spot_summary_lines(analyses, [[1, 2], [3, 4]], lambda a: a.result, unit='nm')
    assert lines[0] == 'Spot 1 (L1 670nm, L2 980nm): 6.84, nan nm'
    assert lines[1] == 'Spot 2 (L4 785nm, L3 670nm): 9.94, 12.28 nm'
    assert lines[2].startswith('Channels not in channel_spots (L5 850nm)')


def test_channel_angle_range_scales_with_dip_width(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [((LAYER_D, 1.456), BUFFER_N)] * 3)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670], bulk_RI_from_TIR=False)
    channel = M.MeasurementChannel(paths[785])
    reflectivity = channel.reflectivity(1)
    angles = reflectivity['angles'].to_numpy(float)
    low, high = M.channel_angle_range(fm, reflectivity, TIR_PARAMETERS)
    SPR_670, below_670, above_670 = M.dip_half_widths(fm.measurement_data['angles'].to_numpy(float), fm.measurement_data['ydata'].to_numpy(float), 63.0)
    SPR_785, below_785, above_785 = M.dip_half_widths(angles, reflectivity['ydata'].to_numpy(float), 63.0)
    assert above_785 < above_670  # The dip is narrower at the longer wavelength
    assert low < SPR_785 < high
    assert (high - SPR_785) == pytest.approx((fm.angle_range[1] - SPR_670) * above_785 / above_670, rel=0.1)
    template_points = np.count_nonzero((fm.measurement_data['angles'] >= fm.angle_range[0]) & (fm.measurement_data['angles'] <= fm.angle_range[1]))
    assert np.count_nonzero((angles >= low) & (angles <= high)) >= template_points // 2


def test_channel_bounds_and_fit_warnings(tmp_path, defaults):
    paths = write_measurement(str(tmp_path), defaults, [((LAYER_D, 1.456), BUFFER_N)] * 3)
    session = make_session()
    fm = template_fresnel(session, defaults, paths[670], bulk_RI_from_TIR=False)
    fm.bounds = [(LAYER_D / 4, -np.inf, 0), (LAYER_D * 1.2, np.inf, 0.1)]  # Tight bounds around the template value
    bounds = M.channel_bounds(fm)
    assert bounds[1][0] >= 10 * fm.fitted_result[0] and bounds[0][0] <= fm.fitted_result[0] / 4
    fresnel_objects, _ = M.fit_fresnel_all_channels(session, fm, defaults)
    assert M.channel_fit_warnings(fresnel_objects[1]) == []
    other = fresnel_objects[1]
    other.fitted_result = np.array([other.bounds[1][0], 0.05, 0.0])
    assert M.channel_fit_warnings(other) == ['fitted value at a bound', 'large offset']
