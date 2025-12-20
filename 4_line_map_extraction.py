from    astropy.table         import Table
from   astropy.convolution    import convolve_fft 
import  numpy                 as     np
from    astropy.io            import fits
from    utils.NII_corrections import nii_ha_ratio_zahid
import  argparse
import  os

def extract_HaHb(full_fits_path, kernel_fits_path, line_fits_path, table_row_path, center_crop=50):

    """
    pass objs from obj_lis to extract ha hb lines
    if objects are downloaded correctely + ha hb lines are present in the data product file,
    return: HDUlist with the following entry:
    0 primary extension, same as original file
    1 line-fit results
    2 segmentation map
    3 clear filter maps
    4,5 Ha line map & line weight + NII correction based on Zahid+2014
    6,7 Hb line map & line weight + PSF matched (also error propagated for weightmap), with convolution kernel generated in upstream workflow
    """
    try:
        with fits.open(full_fits_path) as hdu:
            #set up a crop of 50x50 pix in the center
            shape = hdu[5].shape[0]
            #start index: si and end index: ei
            si = (shape - center_crop) // 2; 
            ei = si + center_crop

            extracted_file = fits.HDUList()
            #save primary extension
            extracted_file.append(hdu[0])
            #save line-fit info
            extracted_file.append(hdu[1])
            
            #select segmentation map [4]
            #also save 1 DSCI image for comparison [5]
            for i in [4,5]: 
                hdu[i].data = hdu[i].data[si:ei,si:ei]
                extracted_file.append(hdu[i])

            #loop to select ha hb line maps
            for image in hdu:
                if image.header.get('EXTTYPE') in ['Ha','Hb'] and (image.name == 'LINE' or image.name == 'LINEWHT'):

                    image.data = image.data[si:ei,si:ei]
                    image.name = f"{image.name}_{image.header['EXTTYPE']}"
                    #NII correction for Ha line map only
                    if image.header.get('EXTTYPE') == 'Ha':
                        if image.name == 'LINE':
                            table_row = Table.read(table_row_path, format='ascii')[0]
                            stellar_mass = table_row['mass']
                            redshift     = table_row['z_MAP']
                            nii_ha_ratio = nii_ha_ratio_zahid(stellar_mass, redshift)
                            image.data = image.data / (1 + nii_ha_ratio)
                            image.header['NII_CORR'] = (str(round(nii_ha_ratio,8)), '[NII]/Ha ratio applied for correction')
                            image.name = 'LINE_HA'
                            
                        elif image.name == 'LINEWHT':
                            image.header['NII_CORR'] = ('applied to LINE_HA', '[NII]/Ha ratio correction applied to LINE_HA')
                            image.data = image.data / (1 + nii_ha_ratio)**2  #propagate weight correction
                            image.name = 'LINEWHT_HA'
                            
                    if image.header.get('EXTTYPE') == 'Hb':
                        with fits.open(kernel_fits_path) as kernel_hdu:
                            kernel_data = kernel_hdu[1].data
                            if image.name == 'LINE':
                                image.data = convolve_fft(image.data, kernel_data, 
                                                        normalize_kernel=True, 
                                                        boundary='wrap',
                                                        nan_treatment='interpolate')
                                image.header['PSF_MATCH'] = ('applied', 'PSF matching convolution applied to Hb line map')
                                image.name = 'LINE_HB'
                            elif image.name == 'LINEWHT':
                                #erorr propagation for convolution on weigth map
                                image.data = 1.0 / convolve_fft(1.0 / image.data, kernel_data**2, 
                                                                normalize_kernel=True, 
                                                                boundary='wrap',
                                                                nan_treatment='interpolate')
                                image.header['PSF_MATCH'] = ('applied', 'PSF matching convolution applied to Hb line map')
                                image.name = 'LINEWHT_HB'
                    extracted_file.append(image)
            extracted_file.writeto(line_fits_path,overwrite=True)
            print(f"Extracted lines saved to {line_fits_path}")
            return None
        
    except Exception as e:
        print(f"! Error processing {full_fits_path}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Extract Ha and Hb line maps from full data products.")
    parser.add_argument('--full_fits_path', type=str, help='Path to the full data product FITS file')
    parser.add_argument('--table_row_path', type=str, default='table_row', help='Path to table row info')
    parser.add_argument('--kernel_fits_path', type=str, help='Path to save the PSF matching kernel FITS file')
    parser.add_argument('--center_crop', type=int, default=50, help='Size of the center crop')
    parser.add_argument('--exist_skip', default=True, help='Skip extraction if file already exists')
    parser.add_argument('--line_fits_path', type=str, help='Path to save the extracted line FITS file')
    args = parser.parse_args()


    if os.path.exists(args.line_fits_path) and args.exist_skip:
        print(f"File {args.line_fits_path} already exists. Skipping extraction.")
        return None
    else:
        extract_HaHb(
            full_fits_path=args.full_fits_path,
            line_fits_path=args.line_fits_path,
            table_row_path=args.table_row_path,
            kernel_fits_path=args.kernel_fits_path,
            center_crop=args.center_crop
        )
        
if __name__ == '__main__':
    main()
