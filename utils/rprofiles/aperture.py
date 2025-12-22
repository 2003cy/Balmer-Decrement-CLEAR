from    astropy.io              import fits
import  numpy                   as     np
from    photutils.aperture      import EllipticalAnnulus, EllipticalAperture
import  astropy.units           as     u

#generate a seiries of ellliptical apertures based on object properties
def ellip_aperture_series(obj, image, annuli_width=1):
    #effective radius 
    re = obj['re']/image.header['PIXASEC'] #in pixels
    # The axis ratio a*b
    q = obj['q']
    # The position angle
    pa = obj['pa']
    # The center of the object
    center = (image.data.shape[0]/2 - 0.5, image.data.shape[1]/2 - 0.5)
    
    # We now generate a series of equi width elliptical annuli
    if 3*re <5:
        radius_limit = 5 #at least 5 pixels
    elif 3*re >=50:
        radius_limit = 50 #cap at 50 pixels to save computation time
    else:
        radius_limit = 3*re
    r = np.arange(0, radius_limit, annuli_width) #up to 3*re
    
    #an elliptical aperture at the center + series of elliptical ring
    center_annuli = [EllipticalAperture(center, r[1], r[1]*q , theta= -(90-pa) * u.deg.to(u.rad))]
    ring_annuli = [EllipticalAnnulus(center, a_in=r[i], a_out=r[i+1], b_in=r[i]*q, b_out=r[i+1]*q, theta= -(90 -pa) * u.deg.to(u.rad) ) for i in range(1, len(r)-1)]
    final_aperture = center_annuli + ring_annuli
    
    return (r[:-1]+r[1:])/2, final_aperture
