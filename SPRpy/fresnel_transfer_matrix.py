# Contains fresnel functions for the transfer-matrix method and curve fitting and TIR angle determination

import numpy as np
import pandas as pd
import bottleneck


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
    :return: ndarray(s), either the fresnel coefficients or the residuals between modelled intensity and measured intensity
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
    n = n_re + 1j * n_im

    # Calculate fresnel coefficients for every angle
    fresnel_coefficients_reflection = np.zeros(len(angles))
    fresnel_coefficients_transmission = np.zeros(len(angles))
    fresnel_coefficients_absorption = np.zeros(len(angles))

    for angle_ind, angle_val in enumerate(angles):

        # Snell's law
        theta = np.zeros(len(n), dtype=np.complex128)
        theta[0] = angle_val * np.pi / 180
        for a in range(len(n) - 1):
            theta[a + 1] = np.real(np.arcsin(n[a] / n[a + 1] * np.sin(theta[a]))) - 1j * np.abs(np.imag(np.arcsin(n[a] / n[a + 1] * np.sin(theta[a]))))

        # Calculating fresnel coefficients:
        fresnel_reflection = np.zeros(len(n) - 1, dtype=np.complex128)
        fresnel_transmission = np.zeros(len(n) - 1, dtype=np.complex128)

        for a in range(len(n) - 1):
            # formulas for s polarization
            s_reflection = (n[a] * np.cos(theta[a]) - n[a + 1] * np.cos(theta[a + 1])) / (n[a] * np.cos(theta[a]) + n[a + 1] * np.cos(theta[a + 1]))
            s_transmission = 2 * n[a] * np.cos(theta[a]) / (n[a] * np.cos(theta[a]) + n[a + 1] * np.cos(theta[a + 1]))

            # formulas for p polarization
            p_reflection = (n[a] * np.cos(theta[a + 1]) - n[a + 1] * np.cos(theta[a])) / (n[a] * np.cos(theta[a + 1]) + n[a + 1] * np.cos(theta[a]))
            p_transmission = 2 * n[a] * np.cos(theta[a]) / (n[a] * np.cos(theta[a + 1]) + n[a + 1] * np.cos(theta[a]))

            fresnel_reflection[a] = s_reflection*(polarization-1) + p_reflection*polarization
            fresnel_transmission[a] = s_transmission*(polarization-1) + p_transmission*polarization

        # Phase shift factors:
        delta = np.zeros(len(n) - 2, dtype=np.complex128)
        for a in range(len(n) - 2):
            delta[a] = 2 * np.pi * layer_thicknesses[a + 1] / wavelength * n[a + 1] * np.cos(theta[a + 1])

        # Build up transfer matrix:
        transfer_matrix = np.array([[1, 0], [0, 1]], dtype=np.complex128)
        for a in range(len(n) - 2):
            transfer_matrix = np.dot(transfer_matrix, np.dot(1 / fresnel_transmission[a], np.dot(np.array([[1, fresnel_reflection[a]], [fresnel_reflection[a], 1]]), np.array([[np.exp(-1j * delta[a]), 0], [0, np.exp(1j * delta[a])]]))))
        transfer_matrix = np.dot(transfer_matrix, np.dot(1 / fresnel_transmission[len(n) - 2], np.array([[1, fresnel_reflection[len(n) - 2]], [fresnel_reflection[len(n) - 2], 1]])))

        # Total fresnel coefficients:
        fr_tot = transfer_matrix[1, 0] / transfer_matrix[0, 0]
        ft_tot = 1 / transfer_matrix[0, 0]

        # Special case of single interface:
        if len(n) == 2:
            fr_tot = fresnel_reflection[0]
            ft_tot = fresnel_transmission[0]

        # Total fresnel coefficients in intensity:
        fresnel_coefficients_reflection[angle_ind] = np.absolute(fr_tot)**2
        fresnel_coefficients_transmission[angle_ind] = np.absolute(ft_tot)**2 * np.real(n[-1] * np.cos(theta[-1])) / np.real(n[0] * np.cos(theta[0]))
        fresnel_coefficients_absorption[angle_ind] = 1 - fresnel_coefficients_reflection[angle_ind] - fresnel_coefficients_transmission[angle_ind]

    # Return fresnel coefficients or residuals depending on if fitting is performed against ydata
    if ydata is None:
        match ydata_type:
            case 'R':
                return fresnel_coefficients_reflection - ydata_offset
            case 'T':
                return fresnel_coefficients_transmission - ydata_offset
            case 'A':
                return fresnel_coefficients_absorption - ydata_offset

    else:
        fresnel_residuals = np.array([]*len(ydata))
        match ydata_type:
            case 'R':
                fresnel_residuals = (fresnel_coefficients_reflection - ydata_offset) - ydata
            case 'T':
                fresnel_residuals = (fresnel_coefficients_transmission - ydata_offset) - ydata
            case 'A':
                fresnel_residuals = (fresnel_coefficients_absorption - ydata_offset) - ydata
        if weights is not None:
            return fresnel_residuals*weights
        else:
            return fresnel_residuals


def TIR_determination(xdata, ydata, TIR_range, scanspeed, _window_count=None, _fit_lower_ind=None, _fit_higher_ind=None):

    # Convert to numpy array first if necessary
    if isinstance(xdata, pd.Series):
        xdata = xdata.to_numpy()

    if isinstance(ydata, pd.Series):
        ydata = ydata.to_numpy()

    TIR_ydata = ydata[(xdata >= TIR_range[0]) & (xdata <= TIR_range[1])]
    TIR_xdata = xdata[(xdata >= TIR_range[0]) & (xdata <= TIR_range[1])]

    # Scanspeed dependent filtering and fitting
    if scanspeed <= 5:
        window_count = (_window_count or 7)
        fit_lower_ind = (_fit_lower_ind or 4)
        fit_higher_ind = (_fit_higher_ind or 5)
    elif scanspeed > 5:
        window_count = (_window_count or 3)
        fit_lower_ind = (_fit_lower_ind or 3)
        fit_higher_ind = (_fit_higher_ind or 3)

    # Filter the data with a moving-average filter to smoothen the signal
    TIR_ydata_filtered = bottleneck.move_mean(TIR_ydata, window=window_count, min_count=1)
    TIR_xdata_filtered = bottleneck.move_mean(TIR_xdata, window=window_count, min_count=1)

    # Find maximum derivative
    deriv_ydata = np.concatenate(([np.diff(TIR_ydata_filtered)[0]], np.diff(TIR_ydata_filtered)))  # Add extra value for dimensions
    dTIR_i = np.argmax(deriv_ydata)

    # Fit against the derivative spike where the derivative is max, considering also nearest neighbors
    poly_c = np.polyfit(TIR_xdata_filtered[dTIR_i-fit_lower_ind:dTIR_i+fit_higher_ind],
                    deriv_ydata[dTIR_i-fit_lower_ind:dTIR_i+fit_higher_ind], 3)

    # Recreate the curve with a lot more points
    deriv_TIR_fit_x = np.linspace(TIR_xdata_filtered[dTIR_i-fit_lower_ind], TIR_xdata_filtered[dTIR_i+fit_lower_ind], 2000)

    # Find TIR from max of deriv fit
    deriv_TIR_fit_y = np.polyval(poly_c, deriv_TIR_fit_x)
    dTIR_final = np.argmax(deriv_TIR_fit_y)
    TIR_theta = deriv_TIR_fit_x[dTIR_final]

    return TIR_theta, deriv_TIR_fit_x, deriv_TIR_fit_y
