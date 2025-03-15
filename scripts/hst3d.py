import numpy as np

lowmass_r = np.array([0.00879765395894433, 0.24633431085043983, 0.49706744868035185, 
              0.7434017595307918, 1.0029325513196483, 1.4956011730205279, 
              1.9970674486803517, 2.4985337243401755])

lowmass_balmer = np.array([3.024000000000001, 2.896000000000001, 2.4960000000000004, 
              1.8719999999999999, 1.4880000000000013, 1.6800000000000015, 
              1.6959999999999997, 2.0000000000000018])

lowmass_err = np.array([3.7119999999999997,3.4880000000000013, 3.087999999999999, 2.464000000000002, 
              1.984, 2.208000000000002, 2.3360000000000003, 2.7200000000000006])

lowmass_err -= lowmass_balmer



midmass_r = np.array([
    0.004132231404958664, 0.2479338842975206, 0.4999999999999999, 
    0.7520661157024793, 0.9958677685950413, 1.495867768595041, 
    1.9917355371900825, 2.495867768595041, 2.995867768595041
])


midmass_balmer = np.array([
    4.526315789473684, 3.6992481203007515, 2.406015037593985, 
    2.1203007518796984, 2.3157894736842106, 2.1052631578947363, 
    2.390977443609021, 2.6466165413533833, 2.3157894736842106
])

midmass_err = np.array([4.96240602, 4.01503759, 2.63157895, 2.33082707, 2.61654135,
        2.40601504, 2.91729323, 3.29323308, 2.88721805])

midmass_err -= midmass_balmer




highmass_r = np.array([
    0.01490683, 0.22732919, 0.45838509, 0.68571429, 0.91677019,
    1.36397516, 1.81118012, 2.26583851, 2.73167702
])

highmass_balmer = np.array([
    9.48074534, 7.99006211, 5.05341615, 3.48819876, 2.8173913,
    2.37018634, 2.11677019, 2.92173913, 1.66956522
])

highmass_err = np.array([
    7.349068322981368, 10.27080745341615, 6.350310559006212, 4.337888198757764, 
    3.5329192546583865, 2.9962732919254655, 2.6683229813664617, 
    3.7416149068323, 1.9677018633540388
])

highmass_err -= highmass_balmer


def hst3d_result(mass_type):
    if mass_type == 'low':
        return lowmass_r, lowmass_balmer+1, lowmass_err
    elif mass_type == 'middle':
        return midmass_r, midmass_balmer+1, midmass_err
    elif mass_type == 'high':
        return highmass_r, highmass_balmer+1, np.abs(highmass_err)
    else:
        raise ValueError("Invalid mass type. Choose from 'low', 'middle', or 'high'.")

# this is the integrated balmer decrement as function of mass from
#A J Battisti, M B Bagley, I Baronchelli, Y S Dai, A L Henry, M A Malkan, A Alavi, D Calzetti, J Colbert, P J McCarthy, V Mehta, M Rafelski, C Scarlata, I Shivaei, E Wisnioski, The average dust attenuation curve at z ∼ 1.3 based on HST grism surveys, Monthly Notices of the Royal Astronomical Society, Volume 513, Issue 3, July 2022, Pages 4431–4450, https://doi.org/10.1093/mnras/stac1052
def integrated_hst3d():
    m = np.array([ 8.41,  8.81,  9.06,  9.27,  9.43,  9.64,  9.88, 10.15, 10.69])
    balmer = np.array([2.93, 4.18, 3.69, 3.6 , 4.38, 4.44, 4.72, 4.74, 4.1 ])
    balmer_err = np.array([0.26, 0.31, 0.32, 0.31, 0.28, 0.31, 0.43, 0.37, 1.29])

    m = np.array([ 9.58009709, 10.03640777, 10.60436893])
    balmer =  np.array([ 0.6047619 ,  1.20952381,  2.34761905])
    balmer_err = np.array([0.61428571,  1.2,2.35714286])
    balmer_err = np.array([1.0952381 ,  1.86190476,  1.3952381 ])


    mask = m>9.0
    return m[mask], balmer[mask], np.abs(balmer_err-balmer)[mask]


def integrated_hst3d_var():
    m = np.array([ 9.25157412,  9.72443936, 10.00539173, 10.43086408])  

    balmer = np.array([ 0.39589151,  0.49876488,  0.7082964 ,  1.1757607])   
    balmer_err = np.array([0.2986401 ,0.33833216,  0.50895897,  1.0396253 ])
    return m, balmer, np.abs(balmer_err-balmer)
