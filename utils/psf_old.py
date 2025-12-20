import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from utils.make_wcs_header import make_wcsheader
from utils.psf import make_subsampled_model_psf
from   drizzlepac        import adrizzle

#get tinytim path either from env variable or set default
if 'TINYTIM' in os.environ:
    tinytim_path = os.environ['TINYTIM']
else:
    tinytim_path = '/u/yacheng/projects/tinytim-7.5'
    print(f"Warning: TINYTIM environment variable not set, using default tinytim_path set in script: {tinytim_path}")


def find_filter_for_ha_hb(redshift):
    '''
    Find the appropriate filter for the given redshift for both HA and HB lines. for g102 or g141
    Since G102 and G141 have overlapping wavelength coverage of ~1075-1150 nm, some redshifts will have both filters available for one given line.
    output: dictionary with 'ha' and 'hb' keys, each containing a list of filters
    '''
    ha_obs, hb_obs = np.array([656.28,486.13]) * (1+redshift) 
    ha_filter = []; hb_filter = []
    if ha_obs  < 1150:
        ha_filter.append('G102')
    if ha_obs > 1075:
        ha_filter.append('G141')
    if hb_obs  < 1150:
        hb_filter.append('G102')
    if hb_obs > 1075:
        hb_filter.append('G141')
    return {'ha':ha_filter, 'hb':hb_filter}

def gen_drizzled_psfs_from_beam(beam_fits_path, row_fits_path,save_psf=True):
    '''
    this section generates tinytim psf for individual observations:
    with unique pixel coordinates, pas, and wavelength for Ha or Hb
    '''
    output_dir = os.path.dirname(beam_fits_path)

    table_row = Table.read(row_fits_path, format='ascii')[0]
    filter_dict_for_line = find_filter_for_ha_hb(table_row['z_MAP'])

    ha_psf_lis = fits.HDUList(); ha_psf_lis.append(fits.PrimaryHDU())
    hb_psf_lis = fits.HDUList(); hb_psf_lis.append(fits.PrimaryHDU())
    wcs_ha_lis = []; wcs_hb_lis = [] #updated wcs info for each psf frame for drizzling
    int_time_ha = []; int_time_hb = [] #integration time for weighting during drizzling

    print(f"Generating PSFs for {table_row['subfield']}_{table_row['ID']} with filters {filter_dict_for_line}")
    #use beam fits to get individual observation info
    with fits.open(beam_fits_path) as hdul:
        for image in hdul:
            if image.name == 'SCI':

                #image cutout center coord
                coord = (image.header['ORIGINX'] - image.header['PAD'] + image.header['NAXIS1']//2,
                         image.header['ORIGINY'] - image.header['PAD'] + image.header['NAXIS2']//2)
                #identifier/rootname for each exposure
                identifier = f"{table_row['subfield']}_{table_row['ID']}_{image.header['ROOTNAME']}"
                rootname = image.header['ROOTNAME']
                
#------------------------------generate individual psfs for each observation------------------------------
                print(f"Processing {identifier} with filter {image.header['FILTER']}")
                #now generate mono psf for both ha hb
                if image.header['FILTER'] in filter_dict_for_line['ha']:
                    tmp_psf_ha = f"{output_dir}/{table_row['subfield']}_{table_row['ID']}_{identifier.split('_')[2]}_ha_tmp.fits"
                    make_subsampled_model_psf(
                                #unique temporary psf path, needed for parallel processing,
                                psf_size = 6,
                                filter_name = image.header['FILTER'],
                                focus = -0.14,  
                                psf_position = coord,
                                mono = 656.28*(1+table_row['z_MAP']),
                                subsampling_factor = 1,
                                tinytim_path = tinytim_path,
                                exist_skip=False)
                    #open tmp psf, update wcs info, and append to hdulist
                    with fits.open(tmp_psf_ha) as hdu_ha:
                        ha_psf = fits.ImageHDU(data=hdu_ha[0].data, header=hdu_ha[0].header)
                        modified_header, wcs_ha = make_wcsheader(
                            ra = table_row['ra'], 
                            dec = table_row['dec'], 
                            size = ha_psf.header['PIXSCALE']*ha_psf.header['NAXIS1'], 
                            pixscale = ha_psf.header['PIXSCALE'], 
                            theta = image.header['PA_APER'])  # 使用实际的PA角度
                        wcs_ha_lis.append(WCS(modified_header))
                        #update header with WCS info, preserving required FITS keywords
                        ha_psf.header.update(modified_header)
                        ha_psf.name = f'{rootname}_ha'
                        ha_psf_lis.append(ha_psf)
                    os.remove(tmp_psf_ha)
                    int_time_ha.append(image.header['DELTATIM'])
                    print(f"  HA PSF generated for {identifier}")
                    
                if image.header['FILTER'] in filter_dict_for_line['hb']:
                    tmp_psf_hb = f"{output_dir}/{table_row['subfield']}_{table_row['ID']}_{identifier.split('_')[2]}_hb_tmp.fits"
                    make_subsampled_model_psf(filename=tmp_psf_hb,
                                psf_size = 6,
                                filter_name = image.header['FILTER'],
                                focus = -0.14,  
                                psf_position = coord,
                                mono = 486.13*(1+table_row['z_MAP']),
                                subsampling_factor = 1,
                                tinytim_path = tinytim_path,
                                exist_skip=False)
                    with fits.open(tmp_psf_hb) as hdu_hb:
                        hb_psf = fits.ImageHDU(data=hdu_hb[0].data, header=hdu_hb[0].header)
                        modified_header, wcs_hb = make_wcsheader(
                            ra = table_row['ra'], 
                            dec = table_row['dec'], 
                            size = hb_psf.header['PIXSCALE']*hb_psf.header['NAXIS1'], 
                            pixscale = hb_psf.header['PIXSCALE'], 
                            theta = image.header['PA_APER']) #positional angle of the beam tile
                        wcs_hb_lis.append(WCS(modified_header))
                        #update header with WCS info, preserving required FITS keywords
                        hb_psf.header.update(modified_header)
                        hb_psf.name = f'{rootname}_hb'
                        hb_psf_lis.append(hb_psf)
                    os.remove(tmp_psf_hb)
                    int_time_hb.append(image.header['DELTATIM'])
                    print(f"  HB PSF generated for {identifier}")
                    
    if save_psf:
        os.makedirs(f'{output_dir}/individual_psf', exist_ok=True)
        ha_psf_lis.writeto(f'{output_dir}/individual_psf/{table_row["subfield"]}_{table_row["ID"]}_ha.fits', overwrite=True)
        hb_psf_lis.writeto(f'{output_dir}/individual_psf/{table_row["subfield"]}_{table_row["ID"]}_hb.fits', overwrite=True)


#----------------------------------------this part drizzles PSFs together------------------------------
    print('Starting to drizzle PSFs together...')
    output_shape = (30, 30)
    target_header, target_wcs = make_wcsheader(
        ra=table_row['ra'],
        dec=table_row['dec'],
        size=output_shape[0] * 0.1,
        pixscale=0.1,
        theta=0,
    )
    drizzled_ha = np.zeros(output_shape, dtype=np.float32)
    drizzled_hb = np.zeros(output_shape, dtype=np.float32)

    # ------------ Ha drizzle ------------
    if len(int_time_ha) > 0:
        weight_ha = np.array(int_time_ha) / np.sum(int_time_ha)
        weight_maps_ha = np.ones_like(ha_psf_lis[1].data)

        for i, (ha, ha_wcs) in enumerate(zip(ha_psf_lis[1:], wcs_ha_lis)):
            print(f"Drizzling HA PSF {i+1}/{len(ha_psf_lis)-1}")
            temp_ha_single =     np.zeros(output_shape, dtype=np.float32)
            temp_ha_wht_single = np.zeros(output_shape, dtype=np.float32)
            temp_ha_con_single = np.zeros(output_shape, dtype=np.int32)

            adrizzle.do_driz(
                insci=ha.data,
                input_wcs=ha_wcs,
                inwht=weight_maps_ha,
                output_wcs=target_wcs,
                outsci=temp_ha_single,
                outwht=temp_ha_wht_single,
                outcon=temp_ha_con_single,
                in_units='cps',
                expin=1.0,
                wt_scl=1.0,
                wcslin_pscale=1,
                pixfrac=0.8,
                kernel='square',
                fillval=0.0,
            )

            drizzled_ha += temp_ha_single * weight_ha[i]

    # ------------ Hb drizzle ------------
    if len(int_time_hb) > 0:
        weight_hb = np.array(int_time_hb) / np.sum(int_time_hb)
        weight_maps_hb = np.ones_like(hb_psf_lis[1].data)

        for i, (hb, hb_wcs) in enumerate(zip(hb_psf_lis[1:], wcs_hb_lis)):
            print(f"Drizzling HB PSF {i+1}/{len(hb_psf_lis)-1}")
            temp_hb_single =     np.zeros(output_shape, dtype=np.float32)
            temp_hb_wht_single = np.zeros(output_shape, dtype=np.float32)
            temp_hb_con_single = np.zeros(output_shape, dtype=np.int32)

            adrizzle.do_driz(
                insci=hb.data,
                input_wcs=hb_wcs,
                inwht=weight_maps_hb,
                output_wcs=target_wcs,
                outsci=temp_hb_single,
                outwht=temp_hb_wht_single,
                outcon=temp_hb_con_single,
                in_units='cps',
                expin=1.0,
                wt_scl=1.0,
                wcslin_pscale=1,
                pixfrac=0.8,
                kernel='square',
                fillval=0.0,
            )

            drizzled_hb += temp_hb_single * weight_hb[i]


    ha_drizzled = fits.ImageHDU(data=drizzled_ha)
    ha_drizzled.header.update(target_header)
    hb_drizzled = fits.ImageHDU(data=drizzled_hb)
    hb_drizzled.header.update(target_header)
    ha_drizzled.writeto(f'{output_dir}/{table_row["subfield"]}_{table_row["ID"]}_psf_ha.fits',overwrite=True)
    hb_drizzled.writeto(f'{output_dir}/{table_row["subfield"]}_{table_row["ID"]}_psf_hb.fits',overwrite=True)
