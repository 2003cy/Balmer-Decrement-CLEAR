import  utils.rprofiles.sbcontrast as sbc
from scipy.ndimage import binary_dilation
from astropy.stats import sigma_clip
from    astropy.table       import Table
from    scripts.tools       import *
import  numpy               as     np
from    astropy.io          import fits
from    astropy.cosmology   import Planck18
import  astropy.units       as     u
from    scipy.constants     import arcsec

from    functools           import reduce
from    tqdm             import tqdm
import  gc                                         
from    scripts.tools       import *
from    utils.rprofiles.packman     import double_packman
import  sys, os
from    astropy.wcs         import WCS
from   photutils.aperture  import EllipticalAnnulus, EllipticalAperture, aperture_photometry, CircularAperture, CircularAnnulus




#generate a seiries of ellliptical apertures based on object properties
def ellip_aperture_series(obj, image):
    # The axis ratio
    q = obj['q']
    # The position angle
    pa = obj['pa']
    # The center of the object
    center = ((image.data.shape[0]-2) / 2, (image.data.shape[1]-2) / 2)
    # We now generate a series of elliptical annuli with the same width
    r = np.linspace(0, image.data.shape[0]/2, int((image.data.shape[0]+1)/2))
    r = r[r/q < ((image.data.shape[0]-1)/2)] 

    #an elliptical aperture at the center + series of elliptical annuli
    center_annuli = [EllipticalAperture(center, r[1], r[1]*q , theta= -(90-pa) * u.deg.to(u.rad))]
    ellip_annuli = [EllipticalAnnulus(center, a_in=r[i], a_out=r[i+1], b_in=r[i]*q, b_out=r[i+1]*q, theta= -(90 -pa) * u.deg.to(u.rad) ) for i in range(1, len(r)-2)]
    final_aperture = center_annuli + ellip_annuli
    
    return (r[:-1]+r[1:])/2, final_aperture[:-1]

def surface_brightness_limit(path_to_mosaic, line_name, pixel_length, mask):
    with fits.open(path_to_mosaic) as hdu:
        for image in hdu:
            moasic = None
            if image.header.get('EXTTYPE') == line_name and (image.name == 'LINE'):
            #if image.name == 'DSCI':
                moasic = image
                break
            
        seg = find_data('SEG',hdu)[1].data
        seg = np.where(seg > 0, 1, 0) #convert to binary mask
        seg[moasic.data == 0] = 1  
    # make the segmentation more aggressive by dilating the mask
        seg = binary_dilation(seg, iterations=5)

        # Pixel area in (arcsec/pixel)²
        pix_scale = moasic.header['PIXASEC']  # arcsec/pixel
        pix_area_scale = pix_scale ** 2
        flux_lim_per_pixel = []
        for i,area in enumerate(surface_area):
            #try:                
                #create a mask to get only the central 80x80 pixels
                start_indexx = (moasic.data.shape[1] - 80) // 2
                start_indexy = (moasic.data.shape[0] - 80) // 2
                end_indexx = start_indexx + 80
                end_indexy = start_indexy + 80

                # Calculate the surface brightness limit using sbcontrast
                limit = sbc.sblimit(image = moasic.data[start_indexy:end_indexy, start_indexx:end_indexx],
                            zeropoint= 1,
                            pix_scale = pix_scale,
                            mask = seg[start_indexy:end_indexy, start_indexx:end_indexx],
                            sigma = 2,#np.min([2,np.nanmedian((ha_r/ ha_r_err))]),
                            scale_arcsec = pix_scale * np.max([1,int(area**0.5)])) #if the area scale is less than 1 pixel, set it to 1 pixel
                
                # Convert to flux limit (unit: 1e-17 erg/s/cm²/pix2)
                flux_lim_per_pixel.append(10**( (limit[0][0] - 1) / -2.5)* pix_area_scale)



def radial_profile(obj, linemap, weight, seg, pixel_length,ha=False):
    
    #print(f"Processing {obj['subfield']}-{obj['ID']} for radial profile extraction")
    
    # Here we try to use elliptical annuli to extract the radial profile
    # The semi major axis in arcsec, transformed to pixel:
    r, apertures = ellip_aperture_series(obj, linemap)
    # Initialize arrays to store the results
    ha_r = np.zeros(len(apertures))
    ha_r_err = np.ones(len(apertures))
    seg_out = np.logical_or(seg != obj['ID'],linemap.data <=0)
    #-------------------------------------------
    # Mask to select the object
    mask_pack = double_packman(linemap.data.shape[0],45,45)
    mask = np.logical_or(seg_out,mask_pack)
    surface_area = [aperture.area_overlap(linemap.data,mask=mask,method='subpixel',subpixels=4) for aperture in apertures]

    # Loop over each annulus and calculate the surface brightness
    for i, annulus in enumerate(apertures):
        phot_table  = aperture_photometry(
            linemap.data, annulus, 
            error=weight.data**-0.5, mask=mask, method='subpixel', subpixels=4)
        ha_r[i]     = phot_table['aperture_sum'][0] / float(surface_area[i])
        ha_r_err[i] = phot_table['aperture_sum_err'][0] / float(surface_area[i])


    #--------------------------------------------
    #Here I use the Pieter van Dokkum's method to calculate the detection limit

    

            #except Exception as e:
            #    print(f"Error calculating flux limit for {obj['subfield']}-{obj['ID']} at radius {r[i]}: {e}")
            #    flux_lim_per_pixel.append(np.inf)
        
        # Convert to flux limit per pixel (unit: 1e-17 erg/s/cm²)
        #clear_output(wait=True)
    print(flux_lim_per_pixel)
    return r*pixel_length, ha_r, ha_r_err, flux_lim_per_pixel


#this function will generate the             print(r, ha_r, ha_r_err, hb_r, hb_r_err, balmer_r, balmer_r_err)
#radial table for a given object
def gen_radial_table_ellip(obj,
                     LINE_HA='LINE_HA',      LINE_HB='LINE_HB_CONV',
                     LINEWHT_HA='LINEWHT_HA',LINEWHT_HB='LINEWHT_HB_CONV'):
    #try:
        #print(f"Processing {obj['subfield']}-{obj['ID']} ")
        path = f"data_extracted/{file_name(obj,prefix='extracted')}"
        with fits.open(path,mode='update') as hdu:
            table_name = 'RAD_PROFILE_ELLIP'
            #check if the file already has the radial table
            #if find_data(table_name,hdu) != None:
            #    return f"{obj['subfield']}-{obj['ID']} already exists"
            
            if find_data('SEG_MOD',hdu) != None:
                seg_map = find_data('SEG_MOD',hdu)[1].data
            else:
                seg_map = find_data('SEG',hdu)[1].data

            #extract the radial profile surface brightness
            r, ha_r, ha_r_err, ha_limit = radial_profile(obj,
                                            linemap      = find_data(LINE_HA,hdu)[1],
                                            weight       = find_data(LINEWHT_HA,hdu)[1],
                                            seg          = seg_map,
                                            pixel_length = obj['pixel_length'],
                                            ha           = True)      
            
            r, hb_r, hb_r_err, hb_limit = radial_profile(obj,
                                            linemap      = find_data(LINE_HB,hdu)[1],
                                            weight       = find_data(LINEWHT_HB,hdu)[1],
                                            seg          = seg_map,
                                            pixel_length = obj['pixel_length'])


            #calculate the balmer decrement
            balmer_r     = ha_r/hb_r
            balmer_r_err =  balmer_r * np.sqrt((ha_r_err / ha_r) ** 2 + (hb_r_err / hb_r) ** 2)
            
            #now calculate the extinction
            #color excess
            E_ba = 2.5*np.log10(balmer_r/2.86)
            #attenutation
            A_ba = E_ba / (K_lambda('Hb')-K_lambda('Ha')) * K_lambda('Ha')

            #columns for the radial table
            cols = [
                fits.Column(name='DISTANCE [kpc]',                       format='E', array=r),
                fits.Column(name='Ha_SURF_BRIGHT [1e-17 erg/s/cm2]',     format='E', array=ha_r),
                fits.Column(name='Ha_SURF_BRIGHT_err [1e-17 erg/s/cm2]', format='E', array=ha_r_err),
                fits.Column(name='Ha_SURF_BRIGHT_LIMIT [1e-17 erg/s/cm2]', format='E', array=ha_limit),
                fits.Column(name='Hb_SURF_BRIGHT [1e-17 erg/s/cm2]',     format='E', array=hb_r),
                fits.Column(name='Hb_SURF_BRIGHT_err [1e-17 erg/s/cm2]', format='E', array=hb_r_err),
                fits.Column(name='Hb_SURF_BRIGHT_LIMIT [1e-17 erg/s/cm2]', format='E', array=hb_limit),
                fits.Column(name='BALMER_DECREM',                        format='E', array=balmer_r),
                fits.Column(name='BALMER_DECREM_ERR',                    format='E',array=balmer_r_err),
                fits.Column(name='E_BV',                                  format='E', array=E_ba),
                fits.Column(name='A_Ha',                                  format='E', array=A_ba)
            ]
            new_table = fits.BinTableHDU.from_columns(cols, name=table_name)
            #save or update table
            save_update(new_table,hdu)
            #save_update(new_table_pix,hdu)
            hdu.flush()


            #clear_output(wait=True)
            return f"{obj['subfield']}-{obj['ID']} processed"
    #except Exception as e:
    #        return f"! {obj['subfield']}-{obj['ID']} failed, error{e}"

