# Contains fresnel function for the transfer-matrix method and curve fitting

import numpy as np
import pandas as pd
import bottleneck


def fresnel_calculation(fitted_var=None,
                        fitted_layer_index=(2, 3),
                        angles=np.linspace(39, 50, 1567),
                        wavelength=670,
                        layer_thicknesses=np.array([np.NaN, 2, 50, 4, np.NaN]),
                        n_re=np.array([1.5202, 3.3105, 0.2238, 1.5, 1.0003]),
                        n_im=np.array([0, 3.4556, 3.9259, 0, 0]),
                        ydata=None,
                        ydata_type='R',
                        polarization=1
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

        # Selecting layer to fit
        match fitted_layer_index[1]:
            case 0:
                print('Invalid fitting variable!')
                return 0
            case 1:
                layer_thicknesses[fitted_layer_index[0]] = fitted_var
            case 2:
                n_re[fitted_layer_index[0]] = fitted_var
            case 3:
                n_im[fitted_layer_index[0]] = fitted_var

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

        if polarization == 0:  # formulas for s polarization
            for a in range(len(n) - 1):
                fresnel_reflection[a] = (n[a] * np.cos(theta[a]) - n[a + 1] * np.cos(theta[a + 1])) / (n[a] * np.cos(theta[a]) + n[a + 1] * np.cos(theta[a + 1]))
                fresnel_transmission[a] = 2 * n[a] * np.cos(theta[a]) / (n[a] * np.cos(theta[a]) + n[a + 1] * np.cos(theta[a + 1]))
        elif polarization == 1:  # formulas for p polarization
            for a in range(len(n) - 1):
                fresnel_reflection[a] = (n[a] * np.cos(theta[a + 1]) - n[a + 1] * np.cos(theta[a])) / (n[a] * np.cos(theta[a + 1]) + n[a + 1] * np.cos(theta[a]))
                fresnel_transmission[a] = 2 * n[a] * np.cos(theta[a]) / (n[a] * np.cos(theta[a + 1]) + n[a + 1] * np.cos(theta[a]))

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
                return fresnel_coefficients_reflection
            case 'T':
                return fresnel_coefficients_transmission
            case 'A':
                return fresnel_coefficients_absorption

    else:
        fresnel_residuals = np.array([]*len(ydata))
        match ydata_type:
            case 'R':
                fresnel_residuals = fresnel_coefficients_reflection - ydata
            case 'T':
                fresnel_residuals = fresnel_coefficients_transmission - ydata
            case 'A':
                fresnel_residuals = fresnel_coefficients_absorption - ydata

        return fresnel_residuals


def TIR_determination(xdata, ydata, TIR_range, scanspeed):

    # Convert to numpy array first if necessary
    if isinstance(xdata, pd.Series):
        xdata = xdata.to_numpy()

    if isinstance(ydata, pd.Series):
        ydata = ydata.to_numpy()

    TIR_ydata = ydata[(xdata >= TIR_range[0]) & (xdata <= TIR_range[1])]
    TIR_xdata = xdata[(xdata >= TIR_range[0]) & (xdata <= TIR_range[1])]

    if scanspeed == 5 or scanspeed == 1:
        # Filter the data with a moving-average filter to smoothen the signal
        TIR_ydata_filtered = bottleneck.move_mean(TIR_ydata, window=7, min_count=1)
        TIR_xdata_filtered = bottleneck.move_mean(TIR_xdata, window=7, min_count=1)

        # Find maximum derivative
        deriv_ydata = np.concatenate(([0], np.diff(TIR_ydata_filtered)))  # Add extra 0 for dimensions
        dTIR_i = np.argmax(deriv_ydata)

        # Fit against the derivative spike where the derivative is max, considering also -3 +2 nearest neighbors
        pd = np.polyfit(TIR_xdata_filtered[dTIR_i-4:dTIR_i+5],
                        deriv_ydata[dTIR_i-4:dTIR_i+5], 3)

        # Recreate the curve with a lot more points
        deriv_TIR_fit_x = np.linspace(TIR_xdata_filtered[dTIR_i-4], TIR_xdata_filtered[dTIR_i+4], 2000)

    elif scanspeed == 10:
        # Filter the data with a moving-average filter to smoothen the signal
        TIR_ydata_filtered = bottleneck.move_mean(TIR_ydata, window=3, min_count=1)
        TIR_xdata_filtered = bottleneck.move_mean(TIR_xdata, window=3, min_count=1)

        # Find maximum derivative
        deriv_ydata = np.concatenate(([0], np.diff(TIR_ydata_filtered)))  # Add extra 0 for dimensions
        dTIR_i = np.argmax(deriv_ydata)

        # Fit against the derivative spike where the derivative is max, considering also -3 +2 nearest neighbors
        pd = np.polyfit(TIR_xdata_filtered[dTIR_i-3:dTIR_i+3],
                        deriv_ydata[dTIR_i-3:dTIR_i+3], 3)

        # Recreate the curve with a lot more points
        deriv_TIR_fit_x = np.linspace(TIR_xdata_filtered[dTIR_i-3], TIR_xdata_filtered[dTIR_i+3], 2000)

    else:
        raise ValueError('Invalid scanspeed value')

    # Find TIR from max of deriv fit
    deriv_TIR_fit_y = np.polyval(pd, deriv_TIR_fit_x)
    dTIR_final = np.argmax(deriv_TIR_fit_y)
    TIR_theta = deriv_TIR_fit_x[dTIR_final]

    return TIR_theta, deriv_TIR_fit_x, deriv_TIR_fit_y
