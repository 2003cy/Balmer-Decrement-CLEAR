import numpy as np

lowmass_r = x_more = np.array([
0.43856921, 1.29237947, 2.15085537
])

lowmass_balmer = np.array([4.46441948, 2.0, 3.64044944])

lowmass_err = y_more = np.array([4.65917603, 1.80524345, 5.33333333])
    

lowmass_err -= lowmass_balmer



midmass_r = np.array([
3.82115086, 4.7029549, 5.55209953, 6.41057543
]) - 3.4012441679626746


midmass_balmer = np.array([
4.2247191, 3.0411985, 3.55805243, 3.11610487
])

midmass_err = np.array([
4.32209738, 3.16104869, 3.72284644, 3.38576779
])

midmass_err -= midmass_balmer




highmass_r = np.array([
    7.22239502332815, 8.090202177293934, 8.948678071539657, 9.825816485225506
]) - 6.7838258164852245

highmass_balmer = np.array([
    2.426966292134832, 2.5617977528089897, 4.134831460674158, 2.9063670411985028
])

highmass_err = np.array([
    2.5617977528089897, 2.711610486891386, 4.719101123595506, 3.3108614232209743

])

highmass_err -= highmass_balmer


def jwst_result(mass_type):
    if mass_type == 'low':
        return lowmass_r, lowmass_balmer, lowmass_err
    elif mass_type == 'middle':
        return midmass_r, midmass_balmer, midmass_err
    elif mass_type == 'high':
        return highmass_r, highmass_balmer, np.abs(highmass_err)
    else:
        return np.array([]), np.array([]), np.array([])
    

def integrated_jwst():
    m = np.array([ 8.41101624,  9.41611003, 10.03068485])
    balmer = np.array([2.76057784, 3.12781194, 3.34633152])
    balmer_err =np.array([3.3888942 , 3.43757367,
        3.93923561])
    balmer_err -= balmer
    return m, balmer, balmer_err
