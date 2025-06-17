import numpy as np

### the following is the SFMS model from eq3 of Whitaker 2014
def whitaker_params(z):
    a = 0.7 - 0.13 * z  # Slope varies with redshift
    b = 1.13 + 0.12 * z  # Intercept varies with redshift
    return a, b

# SFMS for 0.5 < z < 1.0 (Whitaker+2014)
def z_05_10(m_array):
    z = 0.75  # Representative of 0.5 < z < 1.0
    a, b = whitaker_params(z)
    sfr = a * (m_array - 10.5) + b
    sfr_err = np.sqrt((0.09 * (m_array - 10.5))**2 + 0.04**2)  # Error calculation
    return sfr, sfr_err

# SFMS for 1.0 < z < 1.5 (Whitaker+2014)
def z_10_15(m_array):
    z = 1.25  # Representative of 1.0 < z < 1.5
    a, b = whitaker_params(z)
    sfr = a * (m_array - 10.5) + b
    sfr_err = np.sqrt((0.08 * (m_array - 10.5))**2 + 0.03**2)  # Error calculation
    return sfr, sfr_err

# Generate SFMS model
def SFMS_model(m_min, m_max):
    m = np.linspace(m_min, m_max, 100)
    sfr_1, sfr_1_err = z_05_10(m)
    sfr_2, sfr_2_err = z_10_15(m)

    data_points = np.array([
    [9.00985010706638, 0.3750833400895641],
    [9.198715203426124, 0.4684361631605811],
    [9.310492505353318, 0.5224858476843619],
    [9.406852248394005, 0.5863141027169592],
    [9.510920770877945, 0.6501528497779989],
    [9.603426124197002, 0.7678734089077328],
    [9.672805139186295, 0.8316649418407791],
    [9.807708779443255, 0.9200445434298445],
    [9.923340471092077, 1.0083979149478028],
    [10.058244111349037, 1.1457752893653756],
    [10.220128479657387, 1.1900936175083348],
    [10.408993576017131, 1.2834464405793513],
    [10.759743040685226, 1.4505162554904314]])
    m,sfr = data_points.reshape(-1,2).transpose()
    err = 0.6539 - sfr[0]
    sfr_1_l, sfr_1_h = sfr - err, sfr + err

    sfr_1_l, sfr_1_h = sfr_1 - sfr_1_err, sfr_1 + sfr_1_err
    sfr_2_l, sfr_2_h = sfr_2 - sfr_2_err, sfr_2 + sfr_2_err



    return m, sfr_1_l, sfr_1_h, sfr_2_l, sfr_2_h