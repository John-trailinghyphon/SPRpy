# Layer series: simulated dry measurements of two samples (SiO2 -> APS -> PEG) at two wavelengths
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from fresnel_transfer_matrix import fresnel_calculation
from SPRpy_functions import load_config
from SPRpy_classes import Sensor
import SPRpy_series as S

WAVELENGTHS = {1: 670, 4: 785}
TRUE = {'SiO2_1': {'SiO2': 21.0, 'APS': 1.0, 'PEG': 2.0}, 'SiO2_2': {'SiO2': 23.0, 'APS': 0.6, 'PEG': 1.5}}
STEPS = [('clean', None), ('APS', ('APS', 1.42)), ('2kDaPEG', ('PEG', 1.456))]
AIR = 1.0003


@pytest.fixture(scope='module')
def config():
    return load_config()


def stack(defaults, wavelength, layers):
    """
    Prism / Ti / Au / SiO2 / added layers / air with the optical constants of config.toml.
    """
    column = lambda name: [str(layer).strip().lower() for layer in defaults['Layer']].index(name.lower())
    n_key, k_key = 'n_{w}'.format(w=wavelength), 'k_{w}'.format(w=wavelength)
    d = [np.nan, float(defaults['d_nm'][column('Ti')]), float(defaults['d_nm'][column('Au')])]
    n = [float(defaults[n_key][column(name)]) for name in ('Prism', 'Ti', 'Au')]
    k = [float(defaults[k_key][column(name)]) for name in ('Prism', 'Ti', 'Au')]
    for name, thickness, layer_n in layers:
        d.append(thickness); n.append(layer_n); k.append(0.0)
    return np.array(d + [np.nan]), np.array(n + [AIR]), np.array(k + [0.0])


def write_series(folder, defaults):
    angles = np.linspace(39.2, 49.7, 1440)
    paths = []
    for sample_index, (sample, true) in enumerate(TRUE.items()):
        layers = []
        for step_index, (step, added) in enumerate(STEPS):
            silica_n = float(defaults['n_670'][[str(l).lower() for l in defaults['Layer']].index('sio2')])
            if added is None:
                layers = [('SiO2', true['SiO2'], None)]
            else:
                layers = layers + [(added[0], true[added[0]], added[1])]
            for number, wavelength in WAVELENGTHS.items():
                silica = float(defaults['n_{w}'.format(w=wavelength)][[str(l).lower() for l in defaults['Layer']].index('sio2')])
                peg = float(defaults['n_{w}'.format(w=wavelength)][[str(l).lower() for l in defaults['Layer']].index('peg')])
                d, n, k = stack(defaults, wavelength, [(name, thickness, silica if name == 'SiO2' else (peg if name == 'PEG' else layer_n))
                                                       for name, thickness, layer_n in layers])
                reflectivity = fresnel_calculation(angles=angles, wavelength=wavelength, layer_thicknesses=d, n_re=n, n_im=k)
                name = '20260116-1{s}{t}000_{sample}_{step}-L{n}_{w}nm.csv'.format(s=sample_index, t=step_index, sample=sample, step=step, n=number, w=wavelength)
                path = os.path.join(folder, name).replace('\\', '/')
                with open(path, 'w') as file:
                    file.write('# Left most column is Time (min), First row is Angles (deg), Scanspeed=5\n')
                    file.write(';'.join(['0'] + ['{a:.6f}'.format(a=a) for a in angles]) + '\n')
                    file.write(';'.join(['0.1'] + ['{v:.7f}'.format(v=v) for v in reflectivity]) + '\n')
                paths.append(path)
    return paths


def make_session(config):
    session = SimpleNamespace(sensor_instances={}, fresnel_analysis_instances={}, sensor_ID_count=0, fresnel_analysis_ID_count=0,
                              SPR_TIR_fitting_parameters={'TIR range': [40.2, 41.8], 'TIR method': 'savitzky_golay', 'TIR smooth width deg': 0.15,
                                                          'TIR fit halfwidth deg': 0.045, 'TIR fit points': 200, 'Fresnel_angle_range_points': [40, 60]},
                              log=[])
    session.save_sensor = session.save_fresnel_analysis = lambda *args: None
    session.add_log = lambda message, save=False: session.log.append(message)
    return session


def test_parse_series_files(tmp_path):
    names = ['20260116-135311_SiO2_1_APS-L1_670nm.csv', '20260116-134607_SiO2_1_clean-L1_670nm.csv', '20260116-134607_SiO2_1_clean-L2_980nm.csv',
             '20260116-134835_SiO2_2_clean-L1_670nm.csv', '20260116-154527_SiO2_2_2kDaPEG-L1_670nm.csv', 'notes.csv']
    paths = [str(tmp_path / name) for name in names]
    series = S.parse_series_files(paths)
    assert list(series['samples']) == ['SiO2_1', 'SiO2_2']
    assert series['steps'] == ['clean', 'APS', '2kDaPEG']  # Ordered by measurement time
    assert series['wavelengths'] == {1: 670, 2: 980}
    assert series['samples']['SiO2_1']['clean'][2].endswith('SiO2_1_clean-L2_980nm.csv')
    assert series['skipped'] == ['notes.csv']


def test_default_recipe_uses_config_values(config):
    series = {'steps': ['clean', 'APS', 'PEG'], 'wavelengths': {1: 670}}
    recipe = S.default_recipe(series, config['default_sensor_values'])
    assert [row['Step'] for row in recipe] == ['APS', 'PEG']
    assert recipe[0]['n'] == 1.45  # Not in config.toml
    assert recipe[1]['n'] == pytest.approx(1.456)  # PEG from config.toml


def test_layer_series_finds_layer_thicknesses(tmp_path, config):
    defaults = config['default_sensor_values']
    paths = write_series(str(tmp_path), defaults)
    series = S.parse_series_files(paths)
    assert series['steps'] == ['clean', 'APS', '2kDaPEG']
    session = make_session(config)
    template = Sensor(series['samples']['SiO2_1']['clean'][1], 1, defaults, sensor_metal='SiO2')
    recipe = [{'Step': 'APS', 'Layer': 'APS', 'n': 1.42, 'k': 0.0, 'd guess [nm]': 2.0},
              {'Step': '2kDaPEG', 'Layer': 'PEG', 'n': 1.456, 'k': 0.0, 'd guess [nm]': 2.0}]
    options = {'fit_offset': True, 'fit_prism_k': True, 'bulk_from_TIR': False, 'angle_range_points': [80, 120],
               'TIR_ranges': [[60.8, 63.0], [40.2, 41.8]]}
    template.refractive_indices[-1] = AIR
    results = S.run_layer_series(session, series, template, recipe, options, defaults, channel_spots=[[1, 2], [3, 4]])

    assert len(results) == 2 * 3 * 2
    assert len(session.fresnel_analysis_instances) == len(results)
    assert len(session.sensor_instances) == 2 * 2  # One sensor per sample and channel, layers added per step
    assert results.groupby(['Sample', 'Channel'])['Sensor'].nunique().eq(1).all()
    assert sorted(sensor.name for sensor in session.sensor_instances.values())[0] == 'SiO2_1 L1 670nm'
    for _, row in results.iterrows():
        assert row['Value'] == pytest.approx(TRUE[row['Sample']][row['Layer']], abs=0.05), (row['Sample'], row['Step'], row['Channel'])
        assert row['Check'] == ''
    assert set(results['Spot']) == {1, 2}
    # The last step of each sample and channel contains all layers, with the earlier ones at the fitted values
    last = session.fresnel_analysis_instances[int(results['FM'].iloc[-1])]
    assert list(last.sensor_object.optical_parameters['Layers']) == ['Prism', 'Ti', 'Au', 'SiO2', 'APS', 'PEG', 'Bulk']
    assert last.sensor_object.wavelength == 785
    summary = S.series_summary(results)
    assert list(summary['Step']) == ['clean', 'APS', '2kDaPEG']
    assert 'L1 670nm' in summary.columns and 'All' in summary.columns


def test_recipe_must_cover_all_steps(tmp_path, config):
    paths = write_series(str(tmp_path), config['default_sensor_values'])
    series = S.parse_series_files(paths)
    template = Sensor(paths[0], 1, config['default_sensor_values'], sensor_metal='SiO2')
    with pytest.raises(ValueError, match='2kDaPEG'):
        S.run_layer_series(make_session(config), series, template, [{'Step': 'APS', 'Layer': 'APS', 'n': 1.42, 'k': 0, 'd guess [nm]': 2}],
                           {'fit_offset': True, 'fit_prism_k': True, 'bulk_from_TIR': True, 'angle_range_points': [40, 60],
                            'TIR_ranges': [[40.2, 41.8]]}, config['default_sensor_values'])
