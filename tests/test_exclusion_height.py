# Exclusion height (non-interacting probe method) on simulated measurements with a known brush height
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from conftest import gold_sensor, N_PRISM, WAVELENGTH
from fresnel_transfer_matrix import fresnel_calculation, TIR_determination, TIR_shift_overlay
from SPRpy_classes import (ExclusionHeight, calculate_exclusion_height, exclusion_height_crossing,
                           exclusion_height_sensitivity, adaptive_exclusion_height, run_exclusion_height_calculation)

BRUSH = (8.0, 1.38)  # True brush height (nm) and RI
BUFFER_N = 1.3350
PROBE_N = 1.3362  # The probe raises the bulk RI but does not enter the brush
TIR_PARAMETERS = {'TIR range': [60.8, 63.0], 'TIR method': 'savitzky_golay', 'TIR smooth width deg': 0.15,
                  'TIR fit halfwidth deg': 0.045, 'TIR window count': 7, 'points_below_TIR_peak': 4,
                  'points_above_TIR_peak': 5, 'TIR fit points': 2000}


def spectrum(angles, bulk_n, layer=BRUSH):
    d, n, k = gold_sensor(bulk_n=bulk_n, layer=layer)
    return fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d, n_re=n, n_im=k)


def critical_angle(bulk_n):
    return np.degrees(np.arcsin(bulk_n / N_PRISM))


def make_analysis(angles, TIR_method='savitzky_golay', noise=0.0, brush=BRUSH):
    """
    Exclusion height analysis of a simulated measurement: 10 buffer spectra, 10 probe spectra and 10 buffer spectra.
    """
    rng = np.random.default_rng(1)
    buffer, probe = spectrum(angles, BUFFER_N, brush), spectrum(angles, PROBE_N, brush)
    rows = [buffer] * 10 + [probe] * 10 + [buffer] * 10
    ydata_df = pd.DataFrame(np.array(rows) + noise * rng.standard_normal((30, len(angles))))

    d, n, k = gold_sensor(bulk_n=BUFFER_N, layer=brush)
    sensor = SimpleNamespace(refractive_indices=n, extinction_coefficients=k, layer_thicknesses=d,
                             wavelength=WAVELENGTH, data_type='R')
    SPR_angle = angles[np.argmin(buffer)]
    fresnel_object = SimpleNamespace(object_id=1, name='Background', sensor_object=sensor, polarization=1, y_offset=0.0,
                                     ini_guess=np.array([brush[0], 0.0, 0.0]), angle_range=[SPR_angle - 0.6, SPR_angle + 0.6],
                                     measurement_data=pd.DataFrame({'angles': angles, 'ydata': buffer}))
    session = SimpleNamespace(SPR_TIR_fitting_parameters=dict(TIR_PARAMETERS, **{'TIR method': TIR_method}))
    sensorgram = pd.DataFrame({'time': np.arange(30.0), 'SPR angle': np.zeros(30), 'TIR angle': np.zeros(30)})
    analysis = ExclusionHeight(session, fresnel_object, sensorgram, 'simulated.csv', 1, 'Test')
    analysis.injection_points = [(5, 5.0, 0.0), (25, 25.0, 0.0)]
    analysis.buffer_points = [(0, 0.0, 0.0), (10, 10.0, 0.0), (20, 20.0, 0.0), (30, 30.0, 0.0)]
    analysis.probe_points = [(20, 20.0, 0.0), (10, 10.0, 0.0)]  # Points of a pair may be selected in any order
    analysis.height_steps = np.linspace(2, 20, 37)
    return analysis, ydata_df


def exclusion_height(analysis, buffer_index=0):
    buffer_RI = calculate_exclusion_height(analysis, 'buffer', buffer_index)[0]
    probe_RI = calculate_exclusion_height(analysis, 'probe', buffer_index // 2)[0]
    return exclusion_height_crossing(analysis.height_steps, buffer_RI, probe_RI), buffer_RI, probe_RI


def test_TIR_shift_overlay_finds_bulk_step(liquid_angles):
    buffer, probe = spectrum(liquid_angles, BUFFER_N), spectrum(liquid_angles, PROBE_N)
    TIR_buffer = TIR_determination(liquid_angles, buffer, TIR_PARAMETERS)[0]
    TIR_probe = TIR_determination(liquid_angles, probe, TIR_PARAMETERS)[0]
    shift = TIR_shift_overlay(liquid_angles, buffer, probe, TIR_buffer, TIR_probe)
    assert shift == pytest.approx(critical_angle(PROBE_N) - critical_angle(BUFFER_N), abs=3e-4)


def test_TIR_shift_overlay_ignores_scale_and_offset(liquid_angles):
    buffer = spectrum(liquid_angles, BUFFER_N)
    shifted = 0.97 * np.interp(liquid_angles - 0.0937, liquid_angles, buffer) + 0.004
    TIR_buffer = TIR_determination(liquid_angles, buffer, TIR_PARAMETERS)[0]
    shift = TIR_shift_overlay(liquid_angles, buffer, shifted, TIR_buffer, TIR_buffer + 0.09)
    assert shift == pytest.approx(0.0937, abs=2e-4)


def test_exclusion_height_crossing():
    heights = np.array([1.0, 2.0, 3.0, 4.0])
    buffer = np.array([1.30, 1.32, 1.34, 1.36])
    probe = np.array([1.29, 1.32, 1.35, 1.38])  # Crossing at 2 nm (equal), first lower at 3 nm
    height, RI = exclusion_height_crossing(heights, buffer, probe)
    assert 2.0 <= height <= 2.01
    assert RI == pytest.approx(1.32, abs=1e-3)
    assert np.isnan(exclusion_height_crossing(heights, buffer, buffer - 0.01)[0])


@pytest.mark.parametrize('bulk_RI_method', ['edge_overlay', 'TIR_angle'])
@pytest.mark.parametrize('TIR_method', ['savitzky_golay', 'legacy'])
def test_exclusion_height_finds_true_brush_height(liquid_angles, bulk_RI_method, TIR_method):
    analysis, ydata_df = make_analysis(liquid_angles, TIR_method)
    analysis.initialize_model(ydata_df, bulk_RI_method=bulk_RI_method)

    assert analysis.bulk_RI_method == bulk_RI_method
    assert len(analysis.buffer_bulk_RIs) == 2 and len(analysis.probe_bulk_RIs) == 1
    assert np.array(analysis.TIR_shifts) == pytest.approx(critical_angle(PROBE_N) - critical_angle(BUFFER_N), abs=1e-3)
    (height, RI), _, _ = exclusion_height(analysis)
    assert height == pytest.approx(BRUSH[0], abs=1.0)
    assert RI == pytest.approx(BRUSH[1], abs=0.01)


def test_edge_overlay_is_independent_of_TIR_settings(liquid_angles):
    shifts = []
    for smooth_width in (0.05, 0.15, 0.25):
        analysis, ydata_df = make_analysis(liquid_angles, noise=1e-4)
        analysis.SPR_TIR_fitting_parameters['TIR smooth width deg'] = smooth_width
        analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
        shifts.append(analysis.TIR_shifts[0])
    assert np.ptp(shifts) < 2e-4  # deg


def test_sensitivity_matches_recalculated_height(liquid_angles):
    analysis, ydata_df = make_analysis(liquid_angles)
    analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
    (height, _), buffer_RI, probe_RI = exclusion_height(analysis)
    sensitivity = exclusion_height_sensitivity(analysis, 0, height, buffer_RI, probe_RI)

    # Recalculate with a probe bulk RI corresponding to a 2 mdeg larger TIR shift
    TIR_probe = np.arcsin(analysis.probe_bulk_RIs[0] / N_PRISM)
    analysis.probe_bulk_RIs = [N_PRISM * np.sin(TIR_probe + np.radians(2e-3))]
    (shifted_height, _), _, _ = exclusion_height(analysis)
    assert sensitivity == pytest.approx((shifted_height - height) / 2, rel=0.2, abs=0.05)
    assert sensitivity > 0  # A larger bulk RI change gives a larger exclusion height


def test_old_analysis_objects_keep_separate_TIR_angles(liquid_angles):
    analysis, ydata_df = make_analysis(liquid_angles)
    del analysis.bulk_RI_method  # Analysis objects from sessions made with SPRpy <= 1.2
    analysis.initialize_model(ydata_df)
    assert analysis.bulk_RI_method == 'TIR_angle'


@pytest.mark.parametrize('brush', [BRUSH, (120.0, 1.36)])
def test_adaptive_search_finds_thin_and_thick_brushes(liquid_angles, brush):
    analysis, ydata_df = make_analysis(liquid_angles, brush=brush)
    analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
    analysis.height_steps = np.linspace(0, 300, 31)  # Coarse 10 nm sweep over a wide range
    result = adaptive_exclusion_height(analysis, 0, 0)
    assert result['height'] == pytest.approx(brush[0], rel=0.1)
    assert result['RI'] == pytest.approx(brush[1], abs=0.01)
    assert result['sensitivity'] > 0
    assert len(result['heights']) < 40  # The sweep stops at the crossing

    # Same crossing as a fine grid around it
    analysis.height_steps = np.linspace(brush[0] * 0.8, brush[0] * 1.2, 201)
    (grid_height, _), _, _ = exclusion_height(analysis)
    assert result['height'] == pytest.approx(grid_height, abs=0.01 * brush[0])


def test_adaptive_search_ignores_zero_height(liquid_angles):
    analysis, ydata_df = make_analysis(liquid_angles)
    analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
    analysis.height_steps = np.linspace(0, 300, 31)  # First usable step (10 nm) is already past the crossing at 8 nm
    assert adaptive_exclusion_height(analysis, 0, 0)['height'] == pytest.approx(BRUSH[0], abs=1.0)
    analysis.height_steps = np.linspace(50, 300, 11)  # Crossing below the range
    assert np.isnan(adaptive_exclusion_height(analysis, 0, 0)['height'])


@pytest.mark.parametrize('height_search', ['adaptive', 'grid'])
def test_run_exclusion_height_calculation(liquid_angles, height_search):
    analysis, ydata_df = make_analysis(liquid_angles)
    analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
    analysis.height_search = height_search
    analysis.height_steps = np.linspace(2, 20, 37)
    run_exclusion_height_calculation(analysis)
    assert analysis.all_exclusion_results.shape == (2, 2)
    assert analysis.all_exclusion_results[0] == pytest.approx([BRUSH[0]] * 2, abs=1.0)
    assert analysis.mean_exclusion_height_result[0] == pytest.approx(BRUSH[0], abs=1.0)
    assert np.all(np.isfinite(analysis.height_sensitivity))
    assert list(analysis.d_n_pair_dfs[0].columns) == ['height', 'buffer RI', 'probe RI']


def test_two_pass_search_resolves_curves_around_the_crossing(liquid_angles):
    analysis, ydata_df = make_analysis(liquid_angles, brush=(120.0, 1.36))
    analysis.initialize_model(ydata_df, bulk_RI_method='edge_overlay')
    analysis.height_steps = np.linspace(2, 300, 50)  # Resolution 50 = points of the fine pass
    result = adaptive_exclusion_height(analysis, 0, 0)
    heights = result['heights']
    assert len(heights) >= 50
    assert heights.min() < result['height'] < heights.max()
    assert heights.max() - heights.min() < 25  # Fine window of about 1.5 coarse steps (~10 nm), not the whole range
    assert np.max(np.diff(heights)) < 0.5
