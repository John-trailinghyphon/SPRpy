# Contains fresnel function for the transfer-matrix method and curve fitting

import numpy as np


def fresnel_calculation(fitted_var,
                        angles=np.linspace(39, 50, 1567),
                        layer='h_surf', wavelength=670,
                        layer_thicknesses=np.array([np.NaN, 2, 50, 4, np.NaN]),
                        n_re=np.array([1.5202, 3.3105, 0.2238, 1.35, 1.0003]),
                        n_im=np.array([0, 3.4556, 3.9259, 0, 0]),
                        ydata=None,
                        ydata_type='R',
                        polarization=1
                        ):

    """
    Function for calculating fresnel coefficients or for fitting angular reflectivity traces based on the residuals of
    a measurement. By default the function provides the thickness of a monolayer of BSA on gold in air.

    :param angles: ndarray
    :param fitted_var: float
    :param layer: string, determines which layer is fitted, see cases below for valid options
    :param wavelength: int
    :param layer_thicknesses: ndarray
    :param n_re: ndarray
    :param n_im: ndarray
    :param ydata: ndarray (default None), if provided the function will instead return residuals between the modelled intensity and measurement
    :param ydata_type: string, specify if reflectivity ('R'), transmission ('T') or absorption ('A') is fitted against
    :param polarization: int, 1 (default) or 0
    :return: ndarray(s), either the fresnel coefficients or the residuals between modelled intensity and measured intensity
    """

    # Selecting layer to fit
    match layer:
        case 'h_Cr':
            layer_thicknesses[1] = fitted_var
        case 'h_Au':
            layer_thicknesses[2] = fitted_var
        case 'h_SiO2':
            layer_thicknesses[3] = fitted_var
        case 'h_surf':
            layer_thicknesses[-2] = fitted_var
        case 'n_re_prism':
            n_re[0] = fitted_var
        case 'n_re_Cr':
            n_re[1] = fitted_var
        case 'n_im_Cr':
            n_im[1] = fitted_var
        case 'n_re_metal':
            n_re[2] = fitted_var
        case 'n_im_metal':
            n_im[2] = fitted_var
        case 'n_re_surf':
            n_re[-2] = fitted_var
        case 'n_im_surf':
            n_im[-2] = fitted_var
        case 'n_re_sol':
            n_re[-1] = fitted_var
        case 'n_im_sol':
            n_im[-1] = fitted_var
        case 'N_surf':
            n_im[-2] = fitted_var[1]
            n_re[-2] = fitted_var[0]
        case 'H_surf':
            n_im[-2] = fitted_var[1]
            layer_thicknesses[-2] = fitted_var[0]
        case 'N_metal':
            n_im[2] = fitted_var[1]
            n_re[2] = fitted_var[0]

    # # Create complex refractive index (this is unnecessary)
    # n = np.array([0] * len(n_re))
    # for value in range(len(n_re)):
    #     n[value] = complex(n_re[value], n_im[value])

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
            transfer_matrix = np.dot(transfer_matrix, 1 / fresnel_transmission[a] * np.array([[1, fresnel_reflection[a]], [fresnel_reflection[a], 1]]) * np.array([[np.exp(-1j * delta[a]), 0], [0, np.exp(1j * delta[a])]]))
        transfer_matrix = np.dot(transfer_matrix, 1 / fresnel_transmission[len(n) - 1] * np.array([[1, fresnel_reflection[len(n) - 1]], [fresnel_reflection[len(n) - 1], 1]]))

        # Total fresnel coefficients:
        fr_tot = transfer_matrix[2, 1] / transfer_matrix[1, 1]
        ft_tot = 1 / transfer_matrix[1, 1]

        # Special case of single interface:
        if len(n) == 2:
            fr_tot = fresnel_reflection[1]
            ft_tot = fresnel_transmission[1]

        # Total fresnel coefficients in intensity:
        fresnel_coefficients_reflection[angle_ind] = (abs(fr_tot)) ^ 2
        fresnel_coefficients_transmission[angle_ind] = (abs(ft_tot)) ^ 2 * np.real(n[len(n)] * np.cos(theta[len(n)])) / np.real(n[1] * np.cos(theta[1]))
        fresnel_coefficients_absorption[angle_ind] = 1 - fresnel_coefficients_reflection[angle_ind] - fresnel_coefficients_transmission[angle_ind]

    # Return fresnel coefficients or residuals depending on if fitting is performed against ydata
    if ydata is None:
        return fresnel_coefficients_reflection, fresnel_coefficients_transmission, fresnel_coefficients_absorption

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
