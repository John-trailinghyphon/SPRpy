# Contains fresnel functions for the transfer-matrix method and curve fitting, and the TIR and SPR angle determination

import numpy as np
import pandas as pd
import bottleneck
import scipy.interpolate
import scipy.optimize

# Available angle determination methods (keys stored in the session SPR_TIR_fitting_parameters)
TIR_METHODS = {'savitzky_golay': 'Savitzky-Golay derivative (recommended)',
               'legacy': 'Moving average derivative (SPRpy <= 1.2)'}
SPR_METHODS = {'dip_level': 'Adaptive dip window (recommended)',
               'legacy': 'Fixed number of points (SPRpy <= 1.2)'}


def fresnel_calculation(fitted_var=None,
                        fitted_layer_index=(2, 3),
                        angles=np.linspace(39, 50, 1567),
                        wavelength=670,
                        layer_thicknesses=np.array([np.nan, 2, 50, 4, np.nan]),
                        n_re=np.array([1.5202, 3.3105, 0.2238, 1.5, 1.0003]),
                        n_im=np.array([0, 3.4556, 3.9259, 0, 0]),
                        ydata=None,
                        ydata_type='R',
                        weights=None,
                        polarization=1.0,
                        ydata_offset=0,
                        ):

    """
    Function for calculating fresnel coefficients or for fitting angular reflectivity traces based on the residuals of
    a measurement. By default, the function provides the thickness of a monolayer of BSA on gold in air.

    The transfer-matrix product is evaluated for all angles at once (vectorised with numpy), which gives identical
    results to the angle-by-angle loop of SPRpy <= 1.2 but is typically 50-200 times faster. This speeds up every
    fitting procedure in SPRpy (fresnel modelling, batch analysis and exclusion height determination).

    :param fitted_var: variable to be fitted
    :param angles: ndarray
    :param fitted_layer_index: tuple
    :param wavelength: int
    :param layer_thicknesses: ndarray
    :param n_re: ndarray
    :param n_im: ndarray
    :param ydata: ndarray (default None), if provided the function will instead return residuals between the modelled intensity and measurement
    :param ydata_type: string, specify if reflectivity ('R'), transmission ('T') or absorption ('A') is fitted against
    :param polarization: int, 1 (default) or 0
    :return: ndarray(s), either the fresnel coefficients or the residuals between modelled intensity and measurement
    """

    # Check first if fitting is performed or not
    if fitted_var is not None:

        # Selecting main layer to fit
        match fitted_layer_index[1]:
            case 0:
                print('Invalid fitting variable!')
                return 0
            case 1:
                layer_thicknesses[fitted_layer_index[0]] = fitted_var[0]
            case 2:
                n_re[fitted_layer_index[0]] = fitted_var[0]
            case 3:
                n_im[fitted_layer_index[0]] = fitted_var[0]

        # Include fitting intensity offset
        if len(fitted_var) >= 2:
            ydata_offset = fitted_var[1]

        # Include fitting prism extinction value
        if len(fitted_var) == 3:
            n_im[0] = fitted_var[2]

    # Merge real and imaginary refractive indices
    n = np.asarray(n_re, dtype=np.float64) + 1j * np.asarray(n_im, dtype=np.float64)
    d = np.asarray(layer_thicknesses, dtype=np.float64)
    angles = np.atleast_1d(np.asarray(angles, dtype=np.float64))
    layers = len(n)

    # Snell's law, layer by layer for all angles at once (rows: layers, columns: angles)
    theta = np.empty((layers, len(angles)), dtype=np.complex128)
    theta[0] = angles * np.pi / 180
    for a in range(layers - 1):
        arcsin_value = np.arcsin(n[a] / n[a + 1] * np.sin(theta[a]))
        theta[a + 1] = np.real(arcsin_value) - 1j * np.abs(np.imag(arcsin_value))
    cos_theta = np.cos(theta)
    n_col = n[:, np.newaxis]

    # Fresnel coefficients of each interface (rows: interfaces, columns: angles)
    n1, n2 = n_col[:-1], n_col[1:]
    cos1, cos2 = cos_theta[:-1], cos_theta[1:]
    s_reflection = (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    s_transmission = 2 * n1 * cos1 / (n1 * cos1 + n2 * cos2)
    p_reflection = (n1 * cos2 - n2 * cos1) / (n1 * cos2 + n2 * cos1)
    p_transmission = 2 * n1 * cos1 / (n1 * cos2 + n2 * cos1)
    fresnel_reflection = s_reflection * (polarization - 1) + p_reflection * polarization
    fresnel_transmission = s_transmission * (polarization - 1) + p_transmission * polarization

    if layers == 2:
        # Special case of single interface
        fr_tot = fresnel_reflection[0]
        ft_tot = fresnel_transmission[0]
    else:
        # Transfer matrix M = prod_a (1/t_a) [[1, r_a], [r_a, 1]] [[exp(-i delta_a), 0], [0, exp(i delta_a)]], built
        # up element-wise for all angles at once
        m00 = np.ones(len(angles), dtype=np.complex128)
        m01 = np.zeros(len(angles), dtype=np.complex128)
        m10 = np.zeros(len(angles), dtype=np.complex128)
        m11 = np.ones(len(angles), dtype=np.complex128)
        for a in range(layers - 1):
            r = fresnel_reflection[a]
            inv_t = 1 / fresnel_transmission[a]
            if a < layers - 2:
                # Phase shift factors of the layer after the interface
                delta = 2 * np.pi * d[a + 1] / wavelength * n[a + 1] * cos_theta[a + 1]
                e_minus = np.exp(-1j * delta)
                e_plus = np.exp(1j * delta)
            else:
                e_minus = e_plus = 1
            b00 = inv_t * e_minus
            b01 = inv_t * r * e_plus
            b10 = inv_t * r * e_minus
            b11 = inv_t * e_plus
            m00, m01, m10, m11 = (m00 * b00 + m01 * b10, m00 * b01 + m01 * b11,
                                  m10 * b00 + m11 * b10, m10 * b01 + m11 * b11)

        # Total fresnel coefficients
        fr_tot = m10 / m00
        ft_tot = 1 / m00

    # Total fresnel coefficients in intensity
    fresnel_coefficients_reflection = np.absolute(fr_tot) ** 2
    fresnel_coefficients_transmission = np.absolute(ft_tot) ** 2 * np.real(n[-1] * cos_theta[-1]) / np.real(n[0] * cos_theta[0])
    fresnel_coefficients_absorption = 1 - fresnel_coefficients_reflection - fresnel_coefficients_transmission

    match ydata_type:
        case 'T':
            fresnel_coefficients = fresnel_coefficients_transmission
        case 'A':
            fresnel_coefficients = fresnel_coefficients_absorption
        case _:
            fresnel_coefficients = fresnel_coefficients_reflection

    # Return fresnel coefficients or residuals depending on if fitting is performed against ydata
    if ydata is None:
        return fresnel_coefficients - ydata_offset

    fresnel_residuals = (fresnel_coefficients - ydata_offset) - np.asarray(ydata, dtype=np.float64)
    if weights is not None:
        return fresnel_residuals * weights
    return fresnel_residuals


# ----------------------------------------------------------------------------------------------------------------------
# TIR angle determination
# ----------------------------------------------------------------------------------------------------------------------

def _as_numpy(data):
    if isinstance(data, (pd.Series, pd.DataFrame)):
        return data.to_numpy(dtype=np.float64)
    return np.asarray(data, dtype=np.float64)


def TIR_determination(xdata, ydata, SPR_TIR_fitting_parameters):
    """
    Determines the TIR angle of a reflectivity spectrum as the angle of steepest increase in reflectivity within
    SPR_TIR_fitting_parameters['TIR range'].

    The method is selected with SPR_TIR_fitting_parameters['TIR method']:
     - 'savitzky_golay' (default for new sessions): derivative from a Savitzky-Golay (local straight line) filter with
       a width set in degrees, followed by a parabola fitted to the top of the derivative peak. Adapts automatically to
       the angular resolution (scan speed) of the measurement.
     - 'legacy' (default for sessions created with SPRpy <= 1.2): moving-average smoothing, first difference and a
       cubic fitted around the derivative maximum.

    :return: TIR_theta, derivative angles, derivative values, fit angles, fit values
    """

    xdata = _as_numpy(xdata)
    ydata = _as_numpy(ydata)

    if SPR_TIR_fitting_parameters.get('TIR method', 'legacy') == 'savitzky_golay':
        return TIR_determination_savitzky_golay(xdata, ydata, SPR_TIR_fitting_parameters['TIR range'],
                                                smooth_width_deg=SPR_TIR_fitting_parameters.get('TIR smooth width deg', 0.15),
                                                fit_halfwidth_deg=SPR_TIR_fitting_parameters.get('TIR fit halfwidth deg', 0.045),
                                                fit_points=SPR_TIR_fitting_parameters.get('TIR fit points', 200))
    return TIR_determination_legacy(xdata, ydata, SPR_TIR_fitting_parameters)


def TIR_determination_legacy(xdata, ydata, SPR_TIR_fitting_parameters):

    TIR_ydata = ydata[(xdata >= SPR_TIR_fitting_parameters['TIR range'][0]) & (xdata <= SPR_TIR_fitting_parameters['TIR range'][1])]
    TIR_xdata = xdata[(xdata >= SPR_TIR_fitting_parameters['TIR range'][0]) & (xdata <= SPR_TIR_fitting_parameters['TIR range'][1])]

    # Filter the data with a moving-average filter to smoothen the signal
    TIR_ydata_filtered = bottleneck.move_mean(TIR_ydata, window=SPR_TIR_fitting_parameters['TIR window count'], min_count=1)
    TIR_xdata_filtered = bottleneck.move_mean(TIR_xdata, window=SPR_TIR_fitting_parameters['TIR window count'], min_count=1)

    # Find maximum derivative
    deriv_ydata = np.concatenate(([np.diff(TIR_ydata_filtered)[0]], np.diff(TIR_ydata_filtered)))  # Add extra value for dimensions
    dTIR_i = np.argmax(deriv_ydata)

    # Fit against the derivative spike where the derivative is max, considering also nearest neighbors
    below = SPR_TIR_fitting_parameters['points_below_TIR_peak']
    above = SPR_TIR_fitting_parameters['points_above_TIR_peak']
    if dTIR_i - below < 0 or dTIR_i + above >= len(TIR_xdata_filtered):
        raise ValueError('TIR derivative peak too close to the edge of the TIR range.')
    poly_c = np.polyfit(TIR_xdata_filtered[dTIR_i - below:dTIR_i + above + 1], deriv_ydata[dTIR_i - below:dTIR_i + above + 1], 3)

    # Recreate the curve with a lot more points
    deriv_TIR_fit_x = np.linspace(TIR_xdata_filtered[dTIR_i - below], TIR_xdata_filtered[dTIR_i + above], SPR_TIR_fitting_parameters['TIR fit points'])

    # Find TIR from max of deriv fit
    deriv_TIR_fit_y = np.polyval(poly_c, deriv_TIR_fit_x)
    dTIR_final = np.argmax(deriv_TIR_fit_y)
    TIR_theta = deriv_TIR_fit_x[dTIR_final]

    return TIR_theta, TIR_xdata_filtered, deriv_ydata, deriv_TIR_fit_x, deriv_TIR_fit_y


def TIR_determination_savitzky_golay(xdata, ydata, TIR_range, smooth_width_deg=0.15, fit_halfwidth_deg=0.045, fit_points=200):
    """
    TIR angle as the angle of steepest increase in reflectivity within TIR_range (port of TIR_det.m by J. Andersson
    and J. Medin, Chalmers). The derivative is computed with a Savitzky-Golay (local straight-line) filter and a
    parabola is fitted to the top of the derivative peak. The filter widths are given in degrees, so the method adapts
    to the angle step of the measurement.

    :param xdata: incident angles (deg), uniformly spaced
    :param ydata: reflectivity
    :param TIR_range: [min, max] angle (deg) in which the TIR angle is searched
    :param smooth_width_deg: width of the derivative filter (deg). Wider = less noise, but distorts the edge
    :param fit_halfwidth_deg: half-width of the parabola fit around the derivative peak (deg)
    :param fit_points: number of points of the returned fit curve (only used for plotting)
    """

    ydata = np.where(ydata > 0, ydata, np.nan)  # Ignore non-positive values
    dx = np.mean(np.diff(xdata))
    half = max(2, int(round(smooth_width_deg / dx / 2)))
    w = max(2, int(round(fit_halfwidth_deg / dx)))

    # Data in TIR_range, padded by half a filter width on each side so that the derivative is valid over the whole range
    in_range = np.flatnonzero((xdata >= TIR_range[0]) & (xdata <= TIR_range[1]))
    if len(in_range) < 2 * w + 1:
        raise ValueError('Too few data points in TIR range.')
    idx = np.arange(max(in_range[0] - half, 0), min(in_range[-1] + half, len(xdata) - 1) + 1)
    x = xdata[idx]
    y = ydata[idx]
    if len(x) < 2 * half + 1:
        raise ValueError('TIR range too narrow for the derivative filter.')

    # Savitzky-Golay first derivative: slope of a straight line fitted to the 2*half+1 points around each angle
    k = np.arange(-half, half + 1)
    kernel = k / np.sum(k ** 2) / dx
    deriv = np.convolve(y, kernel[::-1], mode='same')
    valid = np.zeros(len(x), dtype=bool)
    valid[half:len(x) - half] = True
    valid &= (x >= TIR_range[0]) & (x <= TIR_range[1])
    deriv[~valid] = np.nan
    if np.all(np.isnan(deriv)):
        raise ValueError('No valid derivative in TIR range (range too narrow for the filter, or no valid data).')

    # Fit a parabola to the derivative peak, centred on the highest point for good numerical conditioning
    ip = int(np.nanargmax(deriv))
    fit_idx = np.arange(max(ip - w, 0), min(ip + w, len(x) - 1) + 1)
    fit_idx = fit_idx[~np.isnan(deriv[fit_idx])]
    x0 = x[ip]
    TIR_theta = x0
    if len(fit_idx) >= 3:
        p = np.polyfit(x[fit_idx] - x0, deriv[fit_idx], 2)
        if p[0] < 0:
            TIR_theta = x0 - p[1] / (2 * p[0])
        else:
            p = None
        # Fall back to the highest point if the parabola has no clear maximum
        if p is None or abs(TIR_theta - x0) > w * dx:
            TIR_theta = x0
    else:
        p = None

    fit_x = np.linspace(x[fit_idx[0]], x[fit_idx[-1]], int(fit_points))
    fit_y = np.polyval(p, fit_x - x0) if p is not None else np.full(len(fit_x), deriv[ip])

    return TIR_theta, x[valid], deriv[valid], fit_x, fit_y


def select_TIR_range(xdata, ydata, TIR_ranges):
    """
    Selects the TIR range (e.g. for liquid or air) that contains the TIR edge of a spectrum: the range with the steepest
    increase in reflectivity. This works for single scans and long measurements alike, unlike guessing the medium from
    the number of scans.

    :param TIR_ranges: list of [min, max] ranges
    :return: the selected range, or None if none of the ranges could be evaluated
    """
    xdata = _as_numpy(xdata)
    ydata = _as_numpy(ydata)
    best_range, best_slope = None, -np.inf
    for TIR_range in TIR_ranges:
        try:
            _, _, deriv, _, _ = TIR_determination_savitzky_golay(xdata, ydata, TIR_range)
        except ValueError:
            continue
        # Normalise with the reflectivity level so that the comparison does not depend on intensity calibration
        level = np.nanmean(ydata[(xdata >= TIR_range[0]) & (xdata <= TIR_range[1])])
        slope = np.nanmax(deriv) / level if level > 0 else np.nanmax(deriv)
        if slope > best_slope:
            best_range, best_slope = [float(TIR_range[0]), float(TIR_range[1])], slope
    return best_range


def TIR_shift_overlay(xdata, y_reference, y, TIR_reference, TIR_guess, window_deg=(0.2, 0.5), search_deg=0.03):
    """
    Shift of the TIR edge of spectrum y relative to spectrum y_reference (deg), found by overlaying the two edges: within
    a window around the TIR edge, y is fitted with a * y_reference(angle - shift) + b (scale and offset absorb small
    changes in reflectivity level).

    The difference between two separately determined TIR angles depends somewhat on how the TIR angle is defined
    (filter widths etc.) when the edge changes shape slightly, as it does between buffer and probe injections. The
    overlay shift does not, which matters for the exclusion height since it depends strongly on the bulk refractive
    index change between buffer and probe.

    :param xdata: incident angles (deg), shared by both spectra
    :param y_reference: reference reflectivity spectrum (e.g. buffer)
    :param y: shifted reflectivity spectrum (e.g. probe)
    :param TIR_reference: TIR angle of y_reference (deg)
    :param TIR_guess: TIR angle of y (deg), used as a starting point and to place the window
    :param window_deg: (below, above) the TIR edge of y (deg) included in the overlay
    :param search_deg: the shift is searched within +-search_deg around TIR_guess - TIR_reference
    :return: shift (deg)
    """
    xdata = _as_numpy(xdata)
    y_reference = _as_numpy(y_reference)
    y = _as_numpy(y)

    valid_reference = np.isfinite(y_reference) & np.isfinite(xdata)
    reference = scipy.interpolate.CubicSpline(xdata[valid_reference], y_reference[valid_reference], extrapolate=False)

    shift_guess = TIR_guess - TIR_reference
    in_window = (xdata >= TIR_guess - window_deg[0]) & (xdata <= TIR_guess + window_deg[1]) & np.isfinite(y)
    x_window, y_window = xdata[in_window], y[in_window]
    if len(x_window) < 10:
        raise ValueError('Too few data points around the TIR edge for overlaying spectra.')

    def residual_sum(shift):
        model = reference(x_window - shift)
        ok = np.isfinite(model)
        design = np.column_stack([model[ok], np.ones(np.count_nonzero(ok))])
        coefficients, *_ = np.linalg.lstsq(design, y_window[ok], rcond=None)
        return np.mean((y_window[ok] - design @ coefficients) ** 2)

    bounds = (shift_guess - search_deg, shift_guess + search_deg)
    result = scipy.optimize.minimize_scalar(residual_sum, bounds=bounds, method='bounded', options={'xatol': 1e-7})
    shift = float(result.x)

    # A minimum at the edge of the search range means the overlay failed, fall back to the TIR angle difference
    if min(abs(shift - bounds[0]), abs(shift - bounds[1])) < 1e-4:
        print('WARNING: overlaying the TIR edges failed, using the difference between TIR angles instead.')
        return float(shift_guess)
    return shift


# ----------------------------------------------------------------------------------------------------------------------
# SPR angle determination
# ----------------------------------------------------------------------------------------------------------------------

def moving_mean(y, window):
    """
    Centred moving mean ignoring NaN values (like MATLAB movmean(y, window, 'omitnan')).
    """
    y = np.asarray(y, dtype=np.float64)
    window = max(1, int(window))
    before = window // 2
    after = window - 1 - before
    valid = ~np.isnan(y)
    values = np.concatenate(([0.0], np.cumsum(np.where(valid, y, 0.0))))
    counts = np.concatenate(([0], np.cumsum(valid)))
    index = np.arange(len(y))
    lo = np.clip(index - before, 0, len(y))
    hi = np.clip(index + after + 1, 0, len(y))
    n_valid = counts[hi] - counts[lo]
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(n_valid > 0, (values[hi] - values[lo]) / np.maximum(n_valid, 1), np.nan)


def SPR_minimum_index(xdata, ydata, min_angle=None, smooth_width_deg=0.1):
    """
    Index of the SPR minimum: the lowest point of the smoothed spectrum above min_angle (e.g. the upper limit of the TIR
    range). Falls back to the lowest point of the whole spectrum. Always returns a valid index, also for spectra without
    a clear dip (e.g. Pd or Pt films in air).
    """
    xdata = _as_numpy(xdata)
    ydata = _as_numpy(ydata)
    dx = abs(np.mean(np.diff(xdata))) if len(xdata) > 1 else 1
    y_search = np.where(ydata > 0, ydata, np.nan)
    if min_angle is not None:
        y_search = np.where(xdata > min_angle, y_search, np.nan)
    if np.all(np.isnan(y_search)):
        y_search = ydata
    smoothed = moving_mean(y_search, max(3, int(round(smooth_width_deg / dx))))
    smoothed[np.isnan(y_search)] = np.nan
    if np.all(np.isnan(smoothed)):
        return int(np.nanargmin(ydata))
    return int(np.nanargmin(smoothed))


def angle_range_around_minimum(xdata, ydata, points_below, points_above, min_angle=None):
    """
    Angle range [low, high] spanning points_below and points_above data points around the SPR minimum, clamped to the
    measured angles (so that it also works when the minimum lies close to the end of the scan, e.g. for Pd films).
    """
    xdata = _as_numpy(xdata)
    min_index = SPR_minimum_index(xdata, ydata, min_angle=min_angle)
    low = min(max(min_index - int(points_below), 0), len(xdata) - 1)
    high = min(max(min_index + int(points_above), 0), len(xdata) - 1)
    return [float(xdata[low]), float(xdata[high])]


def SPR_determination(xdata, ydata, SPR_TIR_fitting_parameters):
    """
    Determines the SPR angle (reflectivity minimum) of a spectrum.

    The method is selected with SPR_TIR_fitting_parameters['SPR method']:
     - 'dip_level' (default for new sessions): cubic fitted to the points below SPR dip level (fraction of the dip
       depth) around the minimum, so that the fit window adapts to the width of the dip (broad at 670 nm, narrow at
       980 nm). The minimum is found analytically from the roots of the derivative of the cubic.
     - 'legacy' (default for sessions created with SPRpy <= 1.2): cubic fitted to a fixed number of points below and
       above the lowest point.

    :return: SPR_theta, fit angles, fit values, angles of the fit window, reflectivity of the fit window
    """
    xdata = _as_numpy(xdata)
    ydata = _as_numpy(ydata)

    if SPR_TIR_fitting_parameters.get('SPR method', 'legacy') == 'dip_level':
        return SPR_determination_dip_level(xdata, ydata,
                                           dip_level=SPR_TIR_fitting_parameters.get('SPR dip level', 0.2),
                                           min_angle=SPR_TIR_fitting_parameters['TIR range'][1],
                                           fit_points=SPR_TIR_fitting_parameters.get('SPR fit points', 200))
    return SPR_determination_legacy(xdata, ydata, SPR_TIR_fitting_parameters)


def SPR_determination_legacy(xdata, ydata, SPR_TIR_fitting_parameters):
    below, above = SPR_TIR_fitting_parameters['sensorgram_angle_range_points']
    min_index = int(np.argmin(ydata))
    if min_index - below < 0 or min_index + above >= len(xdata):
        raise ValueError('SPR minimum too close to the end of the scan for the selected number of fit points.')
    x_window = xdata[min_index - below:min_index + above]
    y_window = ydata[min_index - below:min_index + above]
    polynomial = np.polyfit(x_window, y_window, 3)
    x_selection = np.linspace(xdata[min_index - below], xdata[min_index + above], SPR_TIR_fitting_parameters['SPR fit points'])
    y_polyfit = np.polyval(polynomial, x_selection)
    SPR_theta = x_selection[np.argmin(y_polyfit)]
    return SPR_theta, x_selection, y_polyfit, xdata[min_index - below:min_index + above + 1], ydata[min_index - below:min_index + above + 1]


def SPR_determination_dip_level(xdata, ydata, dip_level=0.2, min_angle=None, fit_points=200, smooth_width_deg=0.1, min_points=10):
    """
    SPR angle as the minimum of a cubic fitted to the bottom of the SPR dip (port of SPR_det.m by J. Andersson and
    J. Medin, Chalmers). The fit window contains the contiguous points around the minimum with a reflectivity below
    Rmin + dip_level*(Rtop - Rmin) in the smoothed spectrum.

    If no proper fit window is found (e.g. no complete dip, like for Pd or Pt films in air), the lowest point of the
    smoothed spectrum is returned instead of raising an error, so that the spectrum stays accessible in the GUI.
    """
    dx = np.mean(np.diff(xdata))

    # Smoothed curve above min_angle, ignoring non-positive values
    y_search = np.where(ydata > 0, ydata, np.nan)
    if min_angle is not None:
        y_search = np.where(xdata > min_angle, y_search, np.nan)
    smoothed = moving_mean(y_search, max(3, int(round(smooth_width_deg / dx))))
    smoothed[np.isnan(y_search)] = np.nan
    if np.all(np.isnan(smoothed)):
        raise ValueError('No data above the TIR range for SPR angle determination.')

    # Locate the dip and the reflectivity level that bounds the fit window
    im = int(np.nanargmin(smoothed))
    R_min = smoothed[im]
    R_top = np.nanmax(smoothed)
    threshold = R_min + dip_level * (R_top - R_min)

    # Contiguous points around the minimum that lie below the threshold (NaN counts as above the threshold)
    above_threshold = ~(smoothed < threshold)
    before = np.flatnonzero(above_threshold[:im + 1])
    after = np.flatnonzero(above_threshold[im:])
    lo = before[-1] + 1 if len(before) else 0
    hi = im + after[0] - 1 if len(after) else len(xdata) - 1
    idx = np.arange(lo, hi + 1)
    idx = idx[~np.isnan(y_search[idx])]

    x0 = xdata[im]
    if len(idx) < min_points:
        # Fall back to the minimum of the smoothed curve
        window = np.arange(max(im - min_points, 0), min(im + min_points, len(xdata) - 1) + 1)
        return x0, np.array([x0]), np.array([smoothed[im]]), xdata[window], ydata[window]

    # Cubic fit with centred and scaled angles for good numerical conditioning
    s = (xdata[idx[-1]] - xdata[idx[0]]) / 2
    u = (xdata[idx] - x0) / s
    p = np.polyfit(u, ydata[idx], 3)

    # Minimum from the roots of the derivative: real, inside the window, and a true minimum (positive curvature)
    roots = np.roots(np.polyder(p))
    roots = np.real(roots[np.abs(np.imag(roots)) < 1e-9])
    roots = roots[(roots >= u.min()) & (roots <= u.max()) & (np.polyval(np.polyder(p, 2), roots) > 0)]
    if len(roots) == 0:
        SPR_theta = x0
    else:
        SPR_theta = x0 + roots[np.argmin(np.polyval(p, roots))] * s

    fit_u = np.linspace(u.min(), u.max(), int(fit_points))
    return SPR_theta, x0 + fit_u * s, np.polyval(p, fit_u), xdata[idx], ydata[idx]
