# Shared test setup: the SPRpy modules import each other by name and read config.toml from the working directory,
# so the SPRpy folder is put on the path and used as working directory.
import os
import sys

import numpy as np
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
SPRPY_DIR = os.path.join(os.path.dirname(TESTS_DIR), 'SPRpy')
sys.path.insert(0, SPRPY_DIR)
sys.path.insert(0, TESTS_DIR)
os.chdir(SPRPY_DIR)

# Sensor values used in the tests (gold sensor at 670 nm, from config.toml.defaults)
N_PRISM = 1.5202
WAVELENGTH = 670


def gold_sensor(bulk_n=1.333, layer=None):
    """
    Layer thicknesses, n and k of a Cr/Au sensor with an optional extra layer (thickness, n) and the bulk medium.
    """
    d = [np.nan, 2.0, 50.0]
    n = [N_PRISM, 3.5295, 0.2238]
    k = [0.0, 4.2683, 3.88]
    if layer is not None:
        d.append(layer[0])
        n.append(layer[1])
        k.append(0.0)
    d.append(np.nan)
    n.append(bulk_n)
    k.append(0.0)
    return np.array(d), np.array(n), np.array(k)


@pytest.fixture
def liquid_angles():
    return np.linspace(58.0, 78.0, 2731)  # Same step (~0.0073 deg) as a medium speed liquid scan


@pytest.fixture
def air_angles():
    return np.linspace(39.0, 50.0, 1500)
