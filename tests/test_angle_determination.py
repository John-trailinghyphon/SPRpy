# TIR and SPR angle determination on modelled spectra with known answers.
import numpy as np
import pandas as pd
import pytest

from conftest import gold_sensor, N_PRISM, WAVELENGTH
from fresnel_transfer_matrix import (fresnel_calculation, TIR_determination, SPR_determination, select_TIR_range,
                                     moving_mean, SPR_minimum_index, angle_range_around_minimum)
from SPRpy_functions import calculate_sensorgram, sensorgram_fit_details, linear_drift, sensorgram_traces, surface_coverage, SENSORGRAM_TRACES

LIQUID_RANGE = [60.8, 63.0]
AIR_RANGE = [40.2, 41.8]


def parameters(TIR_method='savitzky_golay', SPR_method='dip_level', TIR_range=LIQUID_RANGE):
    return {'TIR range': TIR_range, 'TIR method': TIR_method, 'SPR method': SPR_method,
            'TIR smooth width deg': 0.15, 'TIR fit halfwidth deg': 0.045, 'SPR dip level': 0.2,
            'TIR window count': 7, 'points_below_TIR_peak': 4, 'points_above_TIR_peak': 5, 'TIR fit points': 2000,
            'sensorgram_angle_range_points': [70, 70], 'SPR fit points': 4000}


def spectrum(angles, bulk_n=1.333, layer=None):
    d, n, k = gold_sensor(bulk_n=bulk_n, layer=layer)
    return fresnel_calculation(angles=angles, wavelength=WAVELENGTH, layer_thicknesses=d, n_re=n, n_im=k)


def model_SPR_angle(bulk_n=1.333, layer=None):
    fine = np.linspace(64, 76, 120001)
    return fine[np.argmin(spectrum(fine, bulk_n, layer))]


@pytest.mark.parametrize('TIR_method', ['savitzky_golay', 'legacy'])
def test_TIR_angle_follows_bulk_refractive_index(liquid_angles, TIR_method):
    # The TIR angle shift between two bulk refractive indices must match Snell's law
    bulk = [1.3330, 1.3380]
    TIR = [TIR_determination(liquid_angles, spectrum(liquid_angles, n), parameters(TIR_method))[0] for n in bulk]
    expected_shift = np.degrees(np.arcsin(bulk[1] / N_PRISM) - np.arcsin(bulk[0] / N_PRISM))
    assert TIR[1] - TIR[0] == pytest.approx(expected_shift, rel=0.03)
    assert TIR[0] == pytest.approx(np.degrees(np.arcsin(bulk[0] / N_PRISM)), abs=0.05)


def test_savitzky_golay_TIR_is_less_noisy_than_legacy(liquid_angles):
    clean = spectrum(liquid_angles)
    rng = np.random.default_rng(0)
    results = {'savitzky_golay': [], 'legacy': []}
    for _ in range(40):
        noisy = clean + rng.normal(0, 1e-3, len(clean))
        for method in results:
            results[method].append(TIR_determination(liquid_angles, noisy, parameters(method))[0])
    assert np.std(results['savitzky_golay']) < np.std(results['legacy'])


def test_TIR_determination_accepts_pandas(liquid_angles):
    y = spectrum(liquid_angles)
    a = TIR_determination(pd.Series(liquid_angles), pd.Series(y), parameters())[0]
    b = TIR_determination(liquid_angles, y, parameters())[0]
    assert a == b


def test_TIR_range_without_data_raises(liquid_angles):
    with pytest.raises(ValueError):
        TIR_determination(liquid_angles, spectrum(liquid_angles), parameters(TIR_range=[30.0, 31.0]))


@pytest.mark.parametrize('layer', [None, (20.0, 1.45), (100.0, 1.40)])
def test_dip_level_SPR_angle_matches_model_minimum(liquid_angles, layer):
    SPR_theta = SPR_determination(liquid_angles, spectrum(liquid_angles, layer=layer), parameters())[0]
    assert SPR_theta == pytest.approx(model_SPR_angle(layer=layer), abs=0.005)


def test_legacy_SPR_angle_matches_model_minimum(liquid_angles):
    SPR_theta = SPR_determination(liquid_angles, spectrum(liquid_angles), parameters(SPR_method='legacy'))[0]
    assert SPR_theta == pytest.approx(model_SPR_angle(), abs=0.01)


def test_SPR_without_dip_falls_back_to_lowest_point():
    # Spectrum without a dip that ends at the lowest point (like Pd films in air)
    angles = np.linspace(39, 78, 3000)
    y = np.linspace(0.9, 0.3, 3000)
    SPR_theta = SPR_determination(angles, y, parameters(TIR_range=AIR_RANGE))[0]
    assert np.isfinite(SPR_theta)
    assert SPR_theta == pytest.approx(angles[-1], abs=0.1)


def test_select_TIR_range_finds_liquid_and_air():
    angles = np.linspace(39, 78, 5000)
    assert select_TIR_range(angles, spectrum(angles, bulk_n=1.333), [LIQUID_RANGE, AIR_RANGE]) == LIQUID_RANGE
    assert select_TIR_range(angles, spectrum(angles, bulk_n=1.0003), [LIQUID_RANGE, AIR_RANGE]) == AIR_RANGE


def test_moving_mean_ignores_nan():
    y = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    np.testing.assert_allclose(moving_mean(y, 3), [1.5, 1.5, 3.0, 4.5, 4.5])


def test_angle_range_is_clamped_at_the_end_of_the_scan():
    angles = np.linspace(39, 78, 3000)
    y = np.linspace(0.9, 0.3, 3000)  # Minimum at the last point
    low, high = angle_range_around_minimum(angles, y, 40, 60, min_angle=41.8)
    assert high == angles[-1]
    assert low == pytest.approx(angles[-41])
    assert SPR_minimum_index(angles, y, min_angle=41.8) == len(angles) - 1


def test_angle_range_ignores_minimum_below_TIR(air_angles):
    # Low reflectivity below the TIR angle must not be taken for the SPR dip
    y = spectrum(air_angles, bulk_n=1.0003)
    y[:20] = 0.01
    low, high = angle_range_around_minimum(air_angles, y, 40, 60, min_angle=AIR_RANGE[1])
    assert low > AIR_RANGE[1]


def test_sensorgram_tracks_bulk_steps_and_never_returns_nan(liquid_angles):
    bulk = np.concatenate((np.full(10, 1.333), np.full(10, 1.336), np.full(10, 1.333)))
    spectra = np.array([spectrum(liquid_angles, n) for n in bulk])
    spectra[5, :] = np.nan  # A broken spectrum
    time = pd.Series(np.arange(len(bulk)) * 0.1)
    sensorgram = calculate_sensorgram(time, pd.Series(liquid_angles), pd.DataFrame(spectra), parameters())

    assert list(sensorgram.columns) == ['time', 'SPR angle', 'TIR angle']
    assert not sensorgram.isna().any().any()
    assert sensorgram.loc[5, 'SPR angle'] == 0 and sensorgram.loc[5, 'TIR angle'] == 0
    TIR_step = sensorgram.loc[15, 'TIR angle'] - sensorgram.loc[0, 'TIR angle']
    assert TIR_step == pytest.approx(np.degrees(np.arcsin(1.336 / N_PRISM) - np.arcsin(1.333 / N_PRISM)), rel=0.05)
    assert sensorgram.loc[15, 'SPR angle'] > sensorgram.loc[0, 'SPR angle']
    assert sensorgram.loc[25, 'SPR angle'] == pytest.approx(sensorgram.loc[0, 'SPR angle'], abs=1e-6)


@pytest.mark.parametrize('TIR_method', ['savitzky_golay', 'legacy'])
@pytest.mark.parametrize('SPR_method', ['dip_level', 'legacy'])
def test_sensorgram_fit_details(liquid_angles, TIR_method, SPR_method):
    details = sensorgram_fit_details(liquid_angles, spectrum(liquid_angles), parameters(TIR_method, SPR_method))
    for key in ['TIR angle', 'TIR deriv x', 'TIR deriv y', 'TIR deriv fit x', 'TIR deriv fit y', 'SPR angle', 'SPR fit x', 'SPR fit y', 'SPR data x', 'SPR data y']:
        assert key in details
    assert np.isfinite(details['TIR angle']) and np.isfinite(details['SPR angle'])
    assert len(details['TIR deriv x']) == len(details['TIR deriv y'])
    assert details['SPR data x'][0] <= details['SPR angle'] <= details['SPR data x'][-1]


def test_linear_drift():
    time = np.linspace(0, 60, 601)
    signal = 70.0 + 0.002 * time + np.random.default_rng(2).normal(0, 1e-4, len(time))
    signal[300:] += 0.1  # An injection after the baseline
    assert linear_drift(time, signal, [0, 25]) == pytest.approx(0.002, rel=0.05)
    assert linear_drift(time, signal, [25, 0]) == pytest.approx(0.002, rel=0.05)
    assert linear_drift(time, signal, [100, 200]) is None


def test_sensorgram_traces_bulk_correction_and_drift():
    sensorgram = pd.DataFrame({'time': [0.0, 1.0, 2.0], 'SPR angle': [69.0, 69.1, 69.3], 'TIR angle': [61.0, 61.05, 61.1]})
    traces = sensorgram_traces(sensorgram, S_SPR=116, S_TIR=74, layer_thickness=0, decay_length=226)
    assert list(traces) == SENSORGRAM_TRACES
    assert traces['Bulk corrected'] == pytest.approx(sensorgram['SPR angle'] - sensorgram['TIR angle'] * 116 / 74)
    thick = sensorgram_traces(sensorgram, 116, 74, 50, 226)['Bulk corrected']
    assert thick == pytest.approx(sensorgram['SPR angle'] - sensorgram['TIR angle'] * 116 / 74 * np.exp(-2 * 50 / 226))
    drifted = sensorgram_traces(sensorgram, 116, 74, 0, 226, drift={'slope': 0.1})
    assert drifted['SPR angle'] == pytest.approx([69.0, 69.0, 69.1])
    assert sensorgram_traces(sensorgram, 116, None, 0, 226)['Bulk corrected'] == pytest.approx(sensorgram['SPR angle'])  # Missing S_TIR


def test_surface_coverage_units():
    # 0.1 deg with delta = 226 nm, dn/dC = 0.18 cm3/g, S = 116 deg/RIU: 226e-7 cm * 0.1 / (2 * 0.18 * 116) g/cm2 = 54.1 ng/cm2
    assert surface_coverage(0.1, 226, 0.18, 116) == pytest.approx(226e-7 * 0.1 / (2 * 0.18 * 116) * 1e9)
    assert surface_coverage(0.1, 226, 0.18, 116) == pytest.approx(54.12, abs=0.01)
    assert surface_coverage(-0.05, 226, 0.18, 116) < 0


def test_surface_coverage_with_layer_thickness():
    # m_S = d * dtheta / (dn/dC * S * (1 - exp(-2 d / delta))), d in nm -> factor 100 for ng/cm2
    expected = 20e-7 * 0.1 / (0.18 * 116 * (1 - np.exp(-2 * 20 / 226))) * 1e9
    assert surface_coverage(0.1, 226, 0.18, 116, thickness=20) == pytest.approx(expected)
    # A thicker layer extends further out in the evanescent field, so the same shift means more mass
    assert surface_coverage(0.1, 226, 0.18, 116, thickness=20) > surface_coverage(0.1, 226, 0.18, 116)
    # The thin film approximation is the limit d -> 0, and is used for d = 0 or None
    assert surface_coverage(0.1, 226, 0.18, 116, thickness=1e-6) == pytest.approx(surface_coverage(0.1, 226, 0.18, 116), rel=1e-6)
    assert surface_coverage(0.1, 226, 0.18, 116, thickness=0) == surface_coverage(0.1, 226, 0.18, 116, thickness=None) == surface_coverage(0.1, 226, 0.18, 116)
