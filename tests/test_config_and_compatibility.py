# Configuration file handling and compatibility with sessions from older SPRpy versions.
import tomllib

import numpy as np

from conftest import SPRPY_DIR
from SPRpy_functions import load_config, default_SPR_TIR_fitting_parameters, ensure_SPR_TIR_fitting_parameters
from SPRpy_classes import point_pair_indices


def test_shipped_config_has_all_default_entries():
    with open(SPRPY_DIR + '/config.toml', 'rb') as f:
        config = tomllib.load(f)
    with open(SPRPY_DIR + '/config.toml.defaults', 'rb') as f:
        defaults = tomllib.load(f)

    def keys(d, prefix=''):
        result = set()
        for key, value in d.items():
            result.add(prefix + key)
            if isinstance(value, dict):
                result |= keys(value, prefix + key + '.')
        return result

    assert keys(defaults) <= keys(config)


def test_default_sensor_value_columns_have_equal_length():
    config = load_config()
    lengths = {key: len(value) for key, value in config['default_sensor_values'].items()}
    assert len(set(lengths.values())) == 1, lengths
    # Extinction coefficients of the dielectric layers (PEG, protein, DNA) are zero at every wavelength
    for key in ['k_670', 'k_785', 'k_850', 'k_980']:
        assert [float(v) for v in config['default_sensor_values'][key][-3:]] == [0, 0, 0]


def test_missing_config_entries_are_filled_from_defaults(tmp_path):
    (tmp_path / 'config.toml').write_text("session_host = '127.0.0.2'\n[TIR_fitting_parameters]\nTIR_fit_points = 100\n")
    (tmp_path / 'config.toml.defaults').write_text("session_host = '127.0.0.1'\nmax_logical_cores = 0\n"
                                                   "[TIR_fitting_parameters]\nTIR_fit_points = 2000\nTIR_method = 'savitzky_golay'\n")
    config = load_config(str(tmp_path / 'config.toml'), str(tmp_path / 'config.toml.defaults'))
    assert config['session_host'] == '127.0.0.2'  # Own values are kept
    assert config['max_logical_cores'] == 0
    assert config['TIR_fitting_parameters'] == {'TIR_fit_points': 100, 'TIR_method': 'savitzky_golay'}


def test_new_sessions_use_new_methods_and_old_sessions_keep_legacy():
    config = load_config()
    new = default_SPR_TIR_fitting_parameters(config['TIR_fitting_parameters'], config['SPR_fitting_parameters'])
    assert new['TIR method'] == 'savitzky_golay' and new['SPR method'] == 'dip_level'

    old = {'TIR range': [60.8, 63.0], 'TIR window count': 7}
    ensure_SPR_TIR_fitting_parameters(old, config['TIR_fitting_parameters'])
    assert old['TIR method'] == 'legacy' and old['SPR method'] == 'legacy'
    assert old['SPR dip level'] == 0.2


def test_point_pairs_can_be_selected_in_any_order():
    points = [(40, 1.0, 0.1), (10, 0.2, 0.1), (50, 2.0, 0.1), (70, 3.0, 0.1)]
    assert point_pair_indices(points, 0) == (10, 40)
    assert point_pair_indices(points, 2) == (50, 70)
