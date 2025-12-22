def Calzetti(line='Ha'):
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