import numpy as np
import itertools
def K_lambda(line='Ha'):
    """
    Calculate the dust attenuation value k(lambda) from the Calzetti et al. (2000) attenuation curve
    for either H-alpha (656.3 nm) or H-beta (486.1 nm).
    
    Parameters:
        line (str): 'Ha' for H-alpha (656.3 nm) or 'Hb' for H-beta (486.1 nm).
        
    Returns:
        float: Attenuation value k(λ)
    """
    if line == 'Ha':
        wavelength_um = 0.6563  # Hα (656.3 nm)
    elif line == 'Hb':
        wavelength_um = 0.4861  # Hβ (486.1 nm)
    else:
        raise ValueError("Invalid input! Use 'Ha' for H-alpha or 'Hb' for H-beta.")

    if 0.12 <= wavelength_um <= 0.63:
        # UV to optical range (Calzetti 2000)
        k_lambda = 2.659 * (-2.156 + 1.509 / wavelength_um - 0.198 / (wavelength_um ** 2) + 0.011 / (wavelength_um ** 3)) + 4.05
    elif 0.63 < wavelength_um <= 2.2:
        # Near-Infrared range (Calzetti 2000)
        k_lambda = 2.659 * (-1.857 + 1.040 / wavelength_um) + 4.05
    else:
        raise ValueError("Wavelength out of range (0.12 - 2.2 μm).")

    return k_lambda


R0_BALMER = 2.86  # intrinsic Ha/Hb (Case B); keep consistent everywhere

def calculate_attenuation_from_ratio(R, R0=R0_BALMER):
    """
    A_Ha from Balmer ratio R = Ha/Hb.
    """
    R = np.asarray(R, dtype=float)
    kHa = K_lambda('Ha')
    kHb = K_lambda('Hb')
    coef = 2.5 * kHa / (kHb - kHa)
    return coef * np.log10(R / R0)

def calculate_attenuation_err_from_ratio(R, R_err, R0=R0_BALMER):
    """
    1-sigma error propagation for A_Ha given R and sigma_R.
    Uses: sigma_A = coef * (1/ln10) * (sigma_R / R)
    """
    R = np.asarray(R, dtype=float)
    R_err = np.asarray(R_err, dtype=float)

    kHa = K_lambda('Ha')
    kHb = K_lambda('Hb')
    coef = 2.5 * kHa / (kHb - kHa)

    with np.errstate(divide='ignore', invalid='ignore'):
        sigma_A = coef * (1.0 / np.log(10.0)) * (R_err / R)

    # mask invalid cases
    sigma_A = np.where((R > 0) & np.isfinite(R) & np.isfinite(R_err) & (R_err >= 0), sigma_A, np.nan)
    return sigma_A

def log10_flux_err(F, F_err):
    """
    sigma(log10 F) = (1/ln10) * sigma_F / F
    """
    F = np.asarray(F, dtype=float)
    F_err = np.asarray(F_err, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        s = (1.0 / np.log(10.0)) * (F_err / F)
    return np.where((F > 0) & (F_err >= 0) & np.isfinite(F) & np.isfinite(F_err), s, np.nan)

def list_flatten(lis):
    """
    Flatten a nested list from individual profiles.
    """
    return list(itertools.chain.from_iterable(lis))