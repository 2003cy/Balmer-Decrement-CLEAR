import numpy as np

#here we apply a correction of Nii for Ha +[NII] linemap
# #using Andreas L. Faisst 2018
def nii_ha_ratio_faisst(stellar_mass, redshift):
    """
    Calculate the [NII]/Hb flux ratio based on the empirical model by Faisst et al. (2018).
    
    Parameters:
    stellar_mass (float): Stellar mass of the galaxy (in solar masses)
    redshift (float): Redshift (z)
    
    Returns:
    float: Computed [NII]/Hb flux ratio
    """
    # Compute the scaled variables
    x = np.log10(stellar_mass) - 10.5
    y = redshift - 2.3

    # Apply the empirical relation
    log_nii_ha = -0.773 + 1.235 * x - 0.281 * x**2 - 0.020 * y + 0.019 * y**2

    # Convert log value to actual ratio
    return 10**log_nii_ha  


# Zahid+13 Table 1 (KK04): sample, z, Z0, log(M0), gamma
ZAHID13_TABLE1 = (
    ("SDSS", 0.08, 9.121, 8.999, 0.85),
    ("SHELS", 0.29, 9.130, 9.304, 0.77),
    ("DEEP2", 0.78, 9.161, 9.661, 0.65),
    ("Y12", 1.40, 9.06, 9.6, 0.7),
    ("E06", 2.26, 9.06, 9.7, 0.6),
)
# Kewley & Ellison 2008 Table 3: KK04 -> PP04 N2, valid for 8.2 <= KK04 <= 9.2
KE08_KK04_TO_PP04N2 = (916.7484, -309.54480, 35.051680, -1.3188000)
# median redshift of the sample: the Zahid+13 MZR is evaluated at this z for all galaxies (no redshift dependence)
Z_MED = 1.06


def zahid13_metallicity_kk04(logm, redshift):
    """Zahid+13 Eq. 4 MZR (KK04), interpolated linearly in redshift between the epochs."""
    logm, redshift = np.broadcast_arrays(np.asarray(logm, float), np.asarray(redshift, float))
    zs = np.array([row[1] for row in ZAHID13_TABLE1])
    vals = np.array([z0 - np.log10(1.0 + 10 ** (-g * (logm - m0))) for _, _, z0, m0, g in ZAHID13_TABLE1])
    zc = np.clip(redshift, zs[0], zs[-1])
    j = np.clip(np.searchsorted(zs, zc) - 1, 0, len(zs) - 2)
    w = (zc - zs[j]) / (zs[j + 1] - zs[j])
    v0 = np.take_along_axis(vals, j[None], 0)[0]
    v1 = np.take_along_axis(vals, (j + 1)[None], 0)[0]
    return (1 - w) * v0 + w * v1


def nii_ha_ratio_zahid(stellar_mass, redshift):
    """
    Calculate the ([NII]6548+6584)/Ha flux ratio from the Zahid et al. (2013) MZR
    (KK04 -> PP04 N2 with Kewley & Ellison 2008; N2 from Pettini & Pagel 2004).

    Parameters:
    stellar_mass (float): Stellar mass of the galaxy (in solar masses)
    redshift (float): Redshift (z); not used, the MZR is taken at Z_MED

    Returns:
    float: Computed ([NII]6548+6584)/Ha flux ratio
    """
    x = np.clip(zahid13_metallicity_kk04(np.log10(stellar_mass), Z_MED), 8.2, 9.2)
    a, b, c, d = KE08_KK04_TO_PP04N2
    n2 = (a + b * x + c * x**2 + d * x**3 - 8.90) / 0.57
    return 4.0 / 3.0 * 10 ** n2
