# Fitting with the SPR dip angle (fresnel modelling and exclusion height)
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from conftest import gold_sensor, WAVELENGTH
from fresnel_transfer_matrix import (fresnel_calculation, model_dip_angle, measured_dip_angle, angle_range_around_minimum,
                                     fit_spectrum, DIP_ANGLE_WEIGHT)
from SPRpy_classes import FresnelModel, ExclusionHeight, fit_surface_layer

LAYER = (3.0, 1.456)
PARAMETERS = {'TIR range': [60.8, 63.0], 'SPR dip level': 0.2, 'TIR method': 'savitzky_golay', 'TIR smooth width deg': 0.15,
              'TIR fit halfwidth deg': 0.045, 'TIR fit points': 200, 'Fresnel_angle_range_points': [40, 60]}


def exact_minimum(angles, d, n, k, wavelength=WAVELENGTH):
    fine = np.linspace(angles[0], angles[-1], 400001)
    return fine[np.argmin(fresnel_calculation(angles=fine, wavelength=wavelength, layer_thicknesses=d, n_re=n, n_im=k))]


@pytest.mark.parametrize('bulk_n, span', [(1.333, (58.0, 78.0)), (1.0003, (39.0, 50.0))])
def test_model_dip_angle_is_accurate(bulk_n, span):
    d, n, k = gold_sensor(bulk_n=bulk_n, layer=LAYER)
    exact = exact_minimum(np.array(span), d, n, k)
    found = model_dip_angle(None, exact + 0.1, wavelength=WAVELENGTH, layer_thicknesses=d, n_re=n, n_im=k)
    assert found == pytest.approx(exact, abs=5e-5)  # 0.05 mdeg
    # The minimum is also found when the first guess is outside the default search width
    assert model_dip_angle(None, exact + 1.0, wavelength=WAVELENGTH, layer_thicknesses=d, n_re=n, n_im=k) == pytest.approx(exact, abs=5e-5)


def distorted_spectrum(angles):
    """
    Spectrum of a 3 nm layer with a reflectivity tilt across the angles (like imperfect optics or calibration), which the
    model cannot describe.
    """
    d, n, k = gold_sensor(bulk_n=1.333, layer=LAYER)
    reflectivity = fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d, n_re=n, n_im=k)
    return reflectivity + 0.03 * (angles - 69.0)


def make_fresnel_model(angles, ydata, match_dip_angle):
    d, n, k = gold_sensor(bulk_n=1.333, layer=(2.0, LAYER[1]))
    sensor = SimpleNamespace(layer_thicknesses=d, refractive_indices=n, extinction_coefficients=k, wavelength=WAVELENGTH,
                             data_type='R', fitted_layer_index=(3, 1),
                             optical_parameters=pd.DataFrame({'Layers': ['Prism', 'Cr', 'Au', 'Layer', 'Bulk'], 'd [nm]': d, 'n': n, 'k': k}))
    session = SimpleNamespace(SPR_TIR_fitting_parameters=dict(PARAMETERS))
    fm = FresnelModel(session, sensor, 'test.csv', pd.DataFrame({'angles': angles, 'ydata': ydata}), 1, 'test')
    fm.match_dip_angle = match_dip_angle
    fm.update_bulk_RI_from_TIR = False
    fm.fit_offset, fm.fit_prism_k = True, True
    fm.ini_guess = np.array([2.0, 0.0, 0.001])
    fm.bounds = [(0.1, -np.inf, 0), (20, np.inf, 0.1)]
    return fm


def fit_with_range(angles, ydata, below, above, match_dip_angle):
    fm = make_fresnel_model(angles, ydata, match_dip_angle)
    fm.angle_range = angle_range_around_minimum(angles, ydata, below, above, min_angle=63.0)
    fm.model_reflectivity_trace()
    return fm


def test_matching_the_dip_angle_makes_the_result_independent_of_the_angle_range():
    angles = np.linspace(58.0, 78.0, 2731)
    ydata = distorted_spectrum(angles)
    results = {}
    for match in (False, True):
        fits = [fit_with_range(angles, ydata, below, above, match) for below, above in ((40, 60), (60, 40), (80, 120), (30, 30))]
        thicknesses = np.array([fm.fitted_result[0] for fm in fits])
        mismatches = np.array([fm.dip_angle_mismatch for fm in fits]) * 1000  # mdeg
        results[match] = (np.ptp(thicknesses), np.max(np.abs(mismatches)))
    spread_without, mismatch_without = results[False]
    spread_with, mismatch_with = results[True]
    assert mismatch_with < 0.3  # Model dip at the measured SPR angle for every range
    assert mismatch_without > 1.0  # The distortion shifts the reflectivity-only fit
    assert spread_with < spread_without / 3


def test_fit_spectrum_without_dip_weight_is_a_reflectivity_fit():
    angles = np.linspace(66.0, 72.0, 800)
    d, n, k = gold_sensor(bulk_n=1.333, layer=LAYER)
    ydata = fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy())
    kwargs = dict(fitted_layer_index=(3, 1), wavelength=WAVELENGTH, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy(), angles=angles, ydata=ydata)
    plain = fit_spectrum([2.0], [0.1, 20], kwargs)
    matched = fit_spectrum([2.0], [0.1, 20], kwargs, measured_dip_angle=measured_dip_angle(angles, ydata, PARAMETERS), dip_weight=DIP_ANGLE_WEIGHT)
    assert plain['x'][0] == pytest.approx(LAYER[0], abs=1e-3)
    assert matched['x'][0] == pytest.approx(LAYER[0], abs=0.01)  # Undistorted data: the same layer either way


def test_analyses_from_older_sessions_keep_the_reflectivity_fit():
    angles = np.linspace(58.0, 78.0, 2731)
    ydata = distorted_spectrum(angles)
    old = make_fresnel_model(angles, ydata, False)
    del old.match_dip_angle  # Fresnel analyses saved with SPRpy <= 2.0.0
    old.angle_range = angle_range_around_minimum(angles, ydata, 40, 60, min_angle=63.0)
    old.model_reflectivity_trace()
    reference = fit_with_range(angles, ydata, 40, 60, False)
    assert old.fitted_result[0] == pytest.approx(reference.fitted_result[0], abs=1e-6)
    assert np.isfinite(old.dip_angle_mismatch)  # The mismatch is reported also without matching


def test_exclusion_height_fit_places_the_model_dip_at_the_measured_angle():
    angles = np.linspace(66.0, 73.0, 960)
    d, n, k = gold_sensor(bulk_n=1.335, layer=(8.0, 1.38))
    ydata = fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d.copy(), n_re=n.copy(), n_im=k.copy()) + 0.02 * (angles - 69.5)
    sensor = SimpleNamespace(layer_thicknesses=d, refractive_indices=n, extinction_coefficients=k, wavelength=WAVELENGTH, data_type='R')
    analysis = SimpleNamespace(sensor_object=sensor, fresnel_object=SimpleNamespace(y_offset=0.0, ini_guess=np.array([8.0])), polarization=1,
                               SPR_TIR_fitting_parameters=dict(PARAMETERS), fit_offset=False, fit_prism=False, RI_initial_guess=1.35,
                               RI_bounds=[1.2, 1.6], match_dip_angle=True)
    spectrum = pd.DataFrame({'angles': angles, 'reflectivity': ydata})
    RI = fit_surface_layer(analysis, spectrum, 1.335, 8.0)[0]
    d_fit = d.copy(); n_fit = n.copy(); n_fit[-2] = RI
    model_dip = model_dip_angle(None, 69.5, wavelength=WAVELENGTH, layer_thicknesses=d_fit, n_re=n_fit, n_im=k.copy())
    assert model_dip == pytest.approx(measured_dip_angle(angles, ydata, PARAMETERS), abs=3e-4)
