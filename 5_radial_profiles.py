
import  numpy                    as     np
from    astropy.io               import fits
import  astropy.units            as     u
from    astropy.table            import Table
from    astropy.wcs              import WCS
from    photutils.aperture       import  aperture_photometry 
from    utils.rprofiles.aperture import ellip_aperture_series
from    utils.rprofiles.sblimit  import surface_brightness_limit
from    utils.rprofiles.packman  import double_packman
from    utils.atten_curve        import Calzetti
import  os
import  argparse

def gen_radial_profile(obj, image, weight, seg, annuli_width=1):
    '''
    source function to perform elliptical annuli measurement -> radial profile
    '''
    # The semi major axis in arcsec, transformed to pixel:
    r, apertures = ellip_aperture_series(obj, image, annuli_width=annuli_width)

    print(f"\n Generated {len(apertures)} elliptical annuli ")

    # Masking: segmentation mask
    #          image padding mask
    #          double packman mask (Nelson 2014.)
    #          bad pixel mask (bad weight)
    seg_mask = seg.data != obj['ID']
    padding_mask = image.data ==0
    pack_mask = double_packman(image.data.shape[0],45,45)
    bad_pixel = weight.data <=0

    mask = np.logical_or.reduce((seg_mask, padding_mask, bad_pixel, pack_mask))
    error_array = np.where(mask,0, weight.data**-0.5) #will have infinite warning, inf in masked region
    # Initialize arrays to store the results
    sb_lis = np.array([])
    sb_err_lis = np.array([])
    surface_area = np.array([])
    sb_limit_lis = np.array([])
    
    # Loop over each annulus and calculate the surface brightness and detection limit
    for i, annulus in enumerate(apertures):
        phot_table  = aperture_photometry(
            image.data, annulus, error=error_array, 
            mask=mask, method='subpixel', subpixels=4)
        
        area = annulus.area_overlap(
            image.data,
            mask=mask,method='subpixel',subpixels=4)

        sb_limit = surface_brightness_limit(image, seg, area)
        
        sb_lis = np.append(sb_lis, phot_table['aperture_sum'][0] / float(area))
        sb_err_lis = np.append(sb_err_lis, phot_table['aperture_sum_err'][0] / float(area))
        surface_area = np.append(surface_area, area)
        sb_limit_lis = np.append(sb_limit_lis, float(sb_limit))
        
        print(f"Annulus {i}: area = {area}, SB = {sb_lis[i]}, SB_err = {sb_err_lis[i]}, SB_limit = {sb_limit}")
        
    print(f"radius: {r}\nsb:     {sb_lis}\nsb_err: {sb_err_lis}\nlimit:  {sb_limit_lis}\n\n")
    return r*image.header['PIXASEC'], sb_lis, sb_err_lis, sb_limit_lis

#radial table for a given object
def compute_bd_profile(row_fits_path,extracted_fits_path, profile_fits_path,seg_path=None,  annuli_widths=1):
    
    table = Table.read(row_fits_path)
    obj = table[0]
    
    with fits.open(extracted_fits_path,mode='update') as hdul:
        
        image_ha = hdul[4]
        weight_ha= hdul[5]
        image_hb = hdul[6]
        weight_hb= hdul[7]
        seg      = hdul[2]

        print(f"gen radial profile for Ha")
        #extract the radial profile surface brightness
        r, ha, ha_err, ha_limit = gen_radial_profile(
            obj=obj,image=image_ha,weight=weight_ha,seg=seg,annuli_width = annuli_widths)      
        
        print(f"gen radial profile for Hb")
        r, hb, hb_err, hb_limit = gen_radial_profile(
            obj=obj,image=image_hb,weight=weight_hb,seg=seg,annuli_width = annuli_widths)
        
        #calculate the balmer decrement
        balmer_r     = ha/hb
        balmer_r_err = np.sqrt(ha_err**2 + (ha * hb_err / hb)**2) / np.abs(hb)
        
        #now calculate the extinction
        #color excess
        E_bv = 2.5*np.log10(balmer_r/2.86)
        #attenutation
        A_ha = E_bv / (Calzetti('Hb')-Calzetti('Ha')) * Calzetti('Ha')

        #first manual selection based on sn_limit
        mask_r = (r < 1*obj['re']/0.1*obj['pixel_length'])
        mask_sb_limit = (ha > ha_limit) & (hb > hb_limit)
        if len(r[mask_r & mask_sb_limit]) > 0:
                has_profile = 1
        else:
                has_profile = 0
        
        table['distance'] = [r]
        table['ha'] = [ha]
        table['ha_err'] = [ha_err]
        table['ha_limit'] = [ha_limit]
        table['hb'] = [hb]
        table['hb_err'] = [hb_err]
        table['hb_limit'] = [hb_limit]
        table['balmer'] = [balmer_r]
        table['balmer_err'] = [balmer_r_err]
        table['E_bv'] = [E_bv]
        table['A_ha'] = [A_ha]
        table['has_profile'] = [has_profile]
        
        table.write(profile_fits_path, overwrite=True)
        
        print(f"Saved radial profile to {profile_fits_path}")
        
def main():
    parser = argparse.ArgumentParser(description="Extract Ha and Hb radial profiles from extracted.fits")
    parser.add_argument('--extracted_fits_path', type=str, required=True, help='Path to the extracted data product FITS file')
    parser.add_argument('--row_fits_path', type=str, required=True, help='Path to table row info')
    parser.add_argument('--profile_fits_path', type=str, required=True, help='Path to save the radial profile table')
    parser.add_argument('--annuli_width', type=int, default=1, help='width of each elliptical ring apereture, in unit of pixels')
    args = parser.parse_args()
    
    if os.path.exists(args.profile_fits_path):
        print(f"File {args.profile_fits_path} already exists. Skipping extraction.")
        return 
    else:
        try:
            compute_bd_profile(
                row_fits_path=args.row_fits_path,
                extracted_fits_path=args.extracted_fits_path,
                profile_fits_path=args.profile_fits_path,
                annuli_widths=args.annuli_width
                )
        except Exception as e:
            print(f"Error during gen radial profiles for file {args.extracted_fits_path}: {e}")
            with open(args.profile_fits_path, 'w') as f:
                f.write(f"Error in generating radial profile: {e}\n")
                
if __name__ == '__main__':
    main()