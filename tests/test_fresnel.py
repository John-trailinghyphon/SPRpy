# The vectorised Fresnel model must reproduce the angle-by-angle implementation of SPRpy <= 1.2 exactly, and fitting
# with it must recover known layer parameters.
import numpy as np
import pytest
import scipy

import legacy_reference
from conftest import gold_sensor, WAVELENGTH
from fresnel_transfer_matrix import fresnel_calculation

STACKS = {
    'gold in air': gold_sensor(bulk_n=1.0003),
    'gold + layer in water': gold_sensor(bulk_n=1.333, layer=(30.0, 1.45)),
    'single interface': (np.array([np.nan, np.nan]), np.array([1.52, 1.333]), np.array([0.0, 0.0])),
    'palladium in air': (np.array([np.nan, 2.0, 20.0, np.nan]), np.array([1.517, 2.78, 3.03, 1.0003]), np.array([0.0, 2.68, 6.10, 0.0])),
    'absorbing prism, many layers': (np.array([np.nan, 2.0, 50.0, 20.0, 100.0, np.nan]), np.array([1.52, 2.15, 0.2238, 1.46, 1.38, 1.333]), np.array([0.002, 2.36, 3.88, 0.0, 0.0, 0.0])),
}


@pytest.mark.parametrize('stack', STACKS.keys())
@pytest.mark.parametrize('ydata_type', ['R', 'T', 'A'])
@pytest.mark.parametrize('polarization', [1.0, 0.0, 0.5])
def test_matches_legacy_implementation(stack, ydata_type, polarization):
    d, n, k = STACKS[stack]
    angles = np.linspace(39, 78, 800)
    expected = legacy_reference.fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy(),
                                                    ydata_type=ydata_type, polarization=polarization)
    result = fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy(),
                                 ydata_type=ydata_type, polarization=polarization)
    np.testing.assert_allclose(result, expected, rtol=1e-9, atol=1e-12)


def test_residuals_offset_and_weights_match_legacy():
    d, n, k = gold_sensor(bulk_n=1.0003, layer=(4.0, 1.5))
    angles = np.linspace(42, 48, 200)
    ydata = np.linspace(0.1, 0.3, 200)
    weights = np.linspace(1, 2, 200)
    fitted_var = [5.0, 0.01, 0.001]
    kwargs = dict(fitted_layer_index=(3, 1), angles=angles, wavelength=WAVELENGTH, ydata=ydata, weights=weights)
    expected = legacy_reference.fresnel_calculation(list(fitted_var), layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy(), **kwargs)
    result = fresnel_calculation(list(fitted_var), layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy(), **kwargs)
    np.testing.assert_allclose(result, expected, rtol=1e-9, atol=1e-12)


def test_fitted_variable_is_written_to_the_layer_arrays():
    # The GUI relies on the fitted value being written into the sensor arrays (as in SPRpy <= 1.2)
    d, n, k = gold_sensor(bulk_n=1.0003, layer=(4.0, 1.5))
    fresnel_calculation([7.5], fitted_layer_index=(3, 1), angles=np.linspace(42, 48, 10), layer_thicknesses=d, n_re=n, n_im=k)
    assert d[3] == 7.5
    fresnel_calculation([3.9], fitted_layer_index=(2, 3), angles=np.linspace(42, 48, 10), layer_thicknesses=d, n_re=n, n_im=k)
    assert k[2] == 3.9


def test_single_angle():
    d, n, k = gold_sensor(bulk_n=1.0003)
    single = fresnel_calculation(angles=[44.0], layer_thicknesses=d, n_re=n, n_im=k)
    several = fresnel_calculation(angles=[43.0, 44.0, 45.0], layer_thicknesses=d, n_re=n, n_im=k)
    assert single.shape == (1,)
    assert single[0] == pytest.approx(several[1])


def test_energy_conservation_for_lossless_stack():
    d, n, k = np.array([np.nan, 100.0, np.nan]), np.array([1.52, 1.45, 1.333]), np.zeros(3)
    angles = np.linspace(20, 50, 100)  # Below the critical angle
    R = fresnel_calculation(angles=angles, layer_thicknesses=d, n_re=n, n_im=k, ydata_type='R')
    T = fresnel_calculation(angles=angles, layer_thicknesses=d, n_re=n, n_im=k, ydata_type='T')
    np.testing.assert_allclose(R + T, 1.0, atol=1e-10)


def test_fit_recovers_layer_thickness_and_offset():
    true_thickness, true_offset = 6.3, 0.012
    d, n, k = gold_sensor(bulk_n=1.0003, layer=(true_thickness, 1.45))
    angles = np.linspace(42.0, 47.0, 400)
    measured = fresnel_calculation(angles=angles, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy()) - true_offset
    measured += np.random.default_rng(1).normal(0, 2e-4, len(angles))

    d[3] = 2.0  # Initial guess
    result = scipy.optimize.least_squares(fresnel_calculation, [2.0, 0.0], bounds=[(0, -np.inf), (20, np.inf)],
                                          kwargs={'fitted_layer_index': (3, 1), 'angles': angles, 'wavelength': WAVELENGTH,
                                                  'layer_thicknesses': d, 'n_re': n, 'n_im': k, 'ydata': measured})
    assert result.x[0] == pytest.approx(true_thickness, abs=0.1)
    assert result.x[1] == pytest.approx(true_offset, abs=1e-3)
