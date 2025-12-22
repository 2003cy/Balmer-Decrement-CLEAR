from    scipy.ndimage           import binary_dilation
import  sbcontrast              as     sbc
import  numpy                   as     np
from    astropy.io              import fits


#compute the surface brightness limit for an image, mask, andspatial scale
def surface_brightness_limit(image, seg, area):
    seg = np.where((seg.data > 0) | (image.data==0), 1, 0) #convert to binary mask, mask out sources or outside of detetor crop
    # make the segmentation more aggressive by dilating the mask
    seg = binary_dilation(seg, iterations=5)

    # Pixel area in (arcsec/pixel)²
    pix_scale = image.header['PIXASEC']  # arcsec/pixel
    pix_area_scale = pix_scale ** 2
    
    zeropoint = 1 #here it can be set to arbitrary, limit will be converted back to flux.
    # Calculate the surface brightness limit using sbcontrast
    limit = sbc.sblimit(image = image.data,
                zeropoint= zeropoint,
                pix_scale = pix_scale,
                mask = seg,
                sigma = 2,
                scale_arcsec = pix_scale * np.max([1,int(area**0.5)]), #set min area to 1 pix^2
                verbose=False
                )[1][0]
    # Convert to flux limit (unit: 1e-17 erg/s/cm²/pix2)
    return limit