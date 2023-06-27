# Contains fresnel function for the transfer-matrix method and curve fitting

import math
import cmath
import numpy as np


def calculate_fresnel(fitted_var, layer, wavelength, layer_thicknesses, n_re, n_im, angles, ydata=None, ydata_type='R', polarization=1):

    # Selecting layer to fit
    match layer:
        case 'h_Cr':
            layer_thicknesses[2] = fitted_var
        case 'h_Au':
            layer_thicknesses[3] = fitted_var
        case 'h_SiO2':
            layer_thicknesses[4] = fitted_var
        case 'h_surf':
            layer_thicknesses[-2] = fitted_var
        case 'n_re_prism':
            n_re[1] = fitted_var
        case 'n_re_Cr':
            n_re[2] = fitted_var
        case 'n_im_Cr':
            n_im[2] = fitted_var
        case 'n_re_Au':
            n_re[3] = fitted_var
        case 'n_im_Au':
            n_im[3] = fitted_var
        case 'n_re_Cr2':
            n_re[-3] = fitted_var
        case 'n_im_Cr2':
            n_im[-3] = fitted_var
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
        case 'N_Au':
            n_im[3] = fitted_var[1]
            n_re[3] = fitted_var[0]

    # Create complex refractive index
    n = np.array([0]*len(n_re))
    for value in range(len(n_re)):
        n[value] = complex(n_re[value], n_im[value])

    # Calculate fresnel coefficients for every angle
    fresnel_coefficients_reflection = np.array([0]*len(angles))
    fresnel_coefficients_transmission = np.array([0] * len(angles))
    fresnel_coefficients_absorption = np.array([0] * len(angles))

    for angle in angles:

        # Snell's law
        angles_rad = np.array([0] * len(angles))
        angles_rad[0] = angles[angle] * math.pi / 180
        for a in range(len(n)-1):
            angles_rad[a + 1] = (cmath.asin(n[a] / n[a + 1] * cmath.sin(angles_rad[a]))).real - 1j * abs((cmath.asin(n[a] / n[a + 1] * cmath.sin(angles_rad[a]))).imag)

        # Calculating fresnel coefficients:
        fresnel_reflection = np.array([0] * (len(n)-1))
        fresnel_transmission = np.array([0] * (len(n)-1))

        if polarization.real == 0:  # Formulas for s polarization
            for a in range(len(n)-1):
                fresnel_reflection[a] = (n[a] * cmath.cos(angles_rad[a]) - n[a + 1] * cmath.cos(angles_rad[a + 1])) / (n[a] * cmath.cos(angles_rad[a]) + n[a + 1] * cmath.cos(angles_rad[a + 1]))
                fresnel_transmission[a] = 2 * n[a] * cmath.cos(angles_rad[a]) / (n[a] * cmath.cos(angles_rad[a]) + n[a + 1] * cmath.cos(angles_rad[a + 1]))

        elif polarization.real == 1:  # Formulas for p polarization
            for a in range(len(n)-1):
                fresnel_reflection[a] = (n[a] * cmath.cos(angles_rad[a + 1]) - n[a + 1] * cmath.cos(angles_rad[a])) / (n[a] * cmath.cos(angles_rad[a + 1]) + n[a + 1] * cmath.cos(angles_rad[a]))
                fresnel_transmission[a] = 2 * n[a] * cmath.cos(angles_rad[a]) / (n[a] * cmath.cos(angles_rad[a + 1]) + n[a + 1] * cmath.cos(angles_rad[a]))

        # Phase shift factors:
        delta = np.array([0] * (len(n)-2))
        for a in range(len(n)-2):
            delta[a] = 2 * math.pi * layer_thicknesses[a + 1] / wavelength * n[a + 1] * math.cos(angles_rad[a + 1])

        # Build up transfer matrix:
        transfer_matrix = np.identity(2)  # Start with unity matrix
        for a in range(len(n)-2):
            transfer_matrix = transfer_matrix * 1 / fresnel_transmission[a] * np.matrix([[1, fresnel_reflection[a]], [fresnel_reflection[a], 1]]) * np.matrix([[cmath.exp(-1j * delta[a]), 0], [0, cmath.exp(1j * delta[a])]])

        transfer_matrix = transfer_matrix * 1 / fresnel_transmission[len(n) - 1] * np.matrix([[1, fresnel_reflection[len(n) - 1]], [fresnel_reflection[len(n) - 1], 1]])

        # total Fresnel coefficients:
        fr_tot = transfer_matrix[2, 1] / transfer_matrix[1, 1]
        ft_tot = 1 / transfer_matrix[1, 1]

        # special case of single interface:
        if len(n) == 2:
            fr_tot = fresnel_reflection[1]
            ft_tot = fresnel_transmission[1]

        # Total fresnel coefficients in intensity:
        fresnel_coefficients_reflection[angle] = (abs(fr_tot)) ^ 2
        fresnel_coefficients_transmission[angle] = (abs(ft_tot)) ^ 2 * (n[len(n)] * cmath.cos(angles_rad[len(n)])).real / (n[1] * cmath.cos(angles_rad[1])).real
        fresnel_coefficients_absorption[angle] = 1 - fresnel_coefficients_reflection[angle] - fresnel_coefficients_transmission[angle]

    if ydata is None:
        return fresnel_coefficients_reflection, fresnel_coefficients_transmission, fresnel_coefficients_absorption

    else:
        match ydata_type:
            case 'R':
                fresnel_residuals = fresnel_coefficients_reflection - ydata
            case 'T':
                fresnel_residuals = fresnel_coefficients_transmission - ydata
            case 'A':
                fresnel_residuals = fresnel_coefficients_absorption - ydata

        return fresnel_residuals
