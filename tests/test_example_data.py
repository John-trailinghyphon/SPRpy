# End-to-end checks on the example measurement shipped with SPRpy: .spr2 conversion, loading and sensorgram.
# These take ~20 s. Skip them with:  python -m pytest -m "not slow"
import glob
import os
import shutil

import numpy as np
import pytest

from conftest import SPRPY_DIR
from SPRpy_functions import (load_config, load_csv_data, default_SPR_TIR_fitting_parameters,
                             update_TIR_settings_for_measurement, calculate_sensorgram)
from SPRpy_classes import Session, add_sensor_backend

EXAMPLE_SPR2 = os.path.join(SPRPY_DIR, 'example data', 'Example_height_probing_data.spr2')

pytestmark = pytest.mark.slow


@pytest.fixture(scope='module')
def example_csv_files(tmp_path_factory):
    import SPRpy_spr2_to_csv
    folder = tmp_path_factory.mktemp('example')
    spr2_file = str(folder / os.path.basename(EXAMPLE_SPR2))
    shutil.copy(EXAMPLE_SPR2, spr2_file)
    SPRpy_spr2_to_csv.convert_spr2_file(spr2_file)
    return sorted(glob.glob(str(folder / '*.csv')))


def test_conversion_creates_one_file_per_laser(example_csv_files):
    names = [os.path.basename(f) for f in example_csv_files]
    assert names == ['Example_height_probing_data-L1_670nm.csv', 'Example_height_probing_data-L2_980nm.csv',
                     'Example_height_probing_data-L3_670nm.csv', 'Example_height_probing_data-L4_785nm.csv']


@pytest.mark.parametrize('channel', range(4))
def test_example_sensorgram(example_csv_files, channel):
    config = load_config()
    path, scanspeed, time_df, angles_df, ydata_df, reflectivity_df = load_csv_data(path=example_csv_files[channel])
    assert scanspeed == 5
    assert ydata_df.shape == (1012, len(angles_df))
    assert np.all(np.diff(angles_df.to_numpy(dtype=float)) > 0)

    parameters = default_SPR_TIR_fitting_parameters(config['TIR_fitting_parameters'], config['SPR_fitting_parameters'])
    update_TIR_settings_for_measurement(parameters, scanspeed, angles_df, ydata_df, config['TIR_fitting_parameters'])
    assert parameters['TIR range'] == [float(i) for i in config['TIR_fitting_parameters']['TIR_range_water_or_long_measurement']]

    sensorgram = calculate_sensorgram(time_df, angles_df, ydata_df, parameters)
    assert (sensorgram['SPR angle'] > 0).all()
    assert (sensorgram['TIR angle'] > 0).all()
    # Liquid measurement: TIR around 61-62 deg and the SPR angle well above it
    assert sensorgram['TIR angle'].between(60.8, 63.0).all()
    assert (sensorgram['SPR angle'] > sensorgram['TIR angle'] + 2).all()
    # Successive spectra should not jump (the probe injections give gradual changes). The first scans are skipped:
    # in the 980 nm channel the TIR edge only appears after ~6 scans (start of the measurement)
    assert np.abs(np.diff(sensorgram['TIR angle'][10:])).max() < 0.05
    assert np.abs(np.diff(sensorgram['SPR angle'][10:])).max() < 0.05


def test_session_and_sensor_creation(example_csv_files, tmp_path):
    config = load_config()
    parameters = default_SPR_TIR_fitting_parameters(config['TIR_fitting_parameters'], config['SPR_fitting_parameters'])
    session = Session('test', parameters, directory=str(tmp_path).replace('\\', '/'), current_data_path=example_csv_files[0])
    sensor = add_sensor_backend(session, example_csv_files[0].replace('\\', '/'), config['default_sensor_values'])
    assert sensor.wavelength == 670
    assert list(sensor.optical_parameters['Layers']) == ['Prism', 'Cr', 'Au', 'Bulk']
    session.save_all()
    assert os.path.exists(session.location + '/Sensors')
