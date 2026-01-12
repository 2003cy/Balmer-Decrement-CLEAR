import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from utils.make_wcs_header import make_wcsheader
from utils.psf import make_subsampled_model_psf
from drizzlepac import adrizzle
import argparse

# get tinytim path either from env variable or set default
if "TINYTIM" in os.environ:
    tinytim_path = os.environ["TINYTIM"]
else:
    tinytim_path = "/u/yacheng/projects/tinytim-7.5"
    print(
        f"Warning: TINYTIM environment variable not set, using default tinytim_path set in script: {tinytim_path}"
    )


def find_filters_for_wavelength(redshift, rest_wavelength):
    """
    given rest-frame wavelength (nm) and redshift, decide which grism the emission line lies in.
    """
    obs_wavelength = rest_wavelength * (1.0 + redshift)  # nm

    filters = []
    if obs_wavelength < 1150.0:
        filters.append("G102")
    if obs_wavelength > 1075.0:
        filters.append("G141")

    return filters


def gen_drizzled_psf_for_wavelength(
    beam_fits_path,
    row_fits_path,
    wavelength,
    save_fits_path,
    save_individual_psf=True,
    psf_shape=31,
):
    """
    PSF generation & drizzle tool for an arbitrary emission line, for HST grism, Grilzli standart dataset.

    Parameters
    ----------
    beam_fits_path : str
        Path to the beam FITS file (contains SCI extensions per exposure).

    row_fits_path : str
        ASCII table row file; first row must contain at least
        ['subfield', 'ID', 'ra', 'dec', 'z_MAP'].

    wavelength : float
        Rest-frame wavelength in nm (e.g., 656.28 for Hα).

    save_fits_path : str
        Full path for the final drizzled PSF FITS file.

    save_individual_psf : bool, optional
        If True, also save a multi-extension FITS with all individual PSFs
        in a subdirectory 'individual' alongside save_fits_path.
    """
    #try:
    # where outputs (combined + individual) will live
    base_dir = os.path.dirname(save_fits_path)
    if base_dir == "":
        base_dir = "."
    os.makedirs(base_dir, exist_ok=True)

    # read targe table's row
    table_row = Table.read(row_fits_path)[0]
    redshift = table_row["z_MAP"]
    filters_for_line = find_filters_for_wavelength(redshift, wavelength)

    print(
        f"Generating PSFs for {table_row['subfield']}_{table_row['ID']} "
        f"targeting output {save_fits_path} "
        f"at rest wavelength {wavelength:.2f} nm * (1 + {redshift:.4f}) (z={redshift:.4f}), filters={filters_for_line}"
    )
    

    # containers for individual PSFs and their WCS / weights
    psf_list = fits.HDUList()
    psf_list.append(fits.PrimaryHDU())
    wcs_list = []
    int_times = []

    # use beam fits to get individual observation info
    with fits.open(beam_fits_path) as hdul:
        for i,image in enumerate(hdul):
            if image.name != "SCI":
                continue
            
            rootname = f"{image.header['ROOTNAME']}_{int(wavelength * (1.0 + redshift))}_{i}" #unique root extension, pointed to _rates raw observations
            identifier = f"{table_row['subfield']}_{table_row['ID']}_{rootname}"
            
            # check if this exposure/filter covers the line
            if image.header["FILTER"] not in filters_for_line:
                print(f"   Skipping {identifier}: wavelenth not in filter {image.header['FILTER']}")
                continue
            
            print(f"Processing {identifier}, with filter {image.header['FILTER']}")
            # image cutout center coord in detector pixels
            coord = (
                image.header["ORIGINX"] - image.header["PAD"] + image.header["NAXIS1"] // 2,
                image.header["ORIGINY"] - image.header["PAD"] + image.header["NAXIS2"] // 2,
            )
            
            # temporary PSF file (deleted after use)
            tmp_psf = os.path.join(
                base_dir,
                f"{table_row['subfield']}_{table_row['ID']}_{rootname}_tmp.fits",
            )

            make_subsampled_model_psf(
                filename=tmp_psf,
                psf_size=5,  # arcsec; 5arcssec in TINYTIM give 45x45 pixelsc-> later drizzle+crop to 31x31
                filter_name=image.header["FILTER"],
                focus=-0.14,
                psf_position=coord,
                mono=wavelength * (1.0 + redshift),  # observed wavelength in nm
                subsampling_factor=1,
                tinytim_path=tinytim_path,
                exist_skip=False,
            )

            # open tmp psf, update wcs info, and append to hdulist
            with fits.open(tmp_psf) as hdu_psf:
                img_hdu = fits.ImageHDU(data=hdu_psf[0].data, header=hdu_psf[0].header)
                print(f"  PSF shape: {img_hdu.data.shape}")

                # PSF image FOV in arcsec = PIXSCALE * NAXIS1
                size_arcsec = img_hdu.header["PIXSCALE"] * img_hdu.header["NAXIS1"]
                modified_header, wcs_obj = make_wcsheader(
                    ra=table_row["ra"],
                    dec=table_row["dec"],
                    size=size_arcsec,
                    pixscale=img_hdu.header["PIXSCALE"],
                    theta=image.header["PA_APER"],
                )

                wcs_list.append(WCS(modified_header))
                img_hdu.header.update(modified_header)
                img_hdu.name = rootname

                psf_list.append(img_hdu)

            os.remove(tmp_psf)

            # integration time for weighting during drizzle
            int_times.append(image.header["DELTATIM"])
            print(f"  PSF generated for {identifier}")

    if len(psf_list) <= 1:
        print("No PSFs generated (no matching filters / exposures). Nothing to drizzle.")
        return

    # optionally save individual PSFs as a multi-extension FITS
    if save_individual_psf:
        indiv_dir = os.path.join(base_dir, "individual")
        os.makedirs(indiv_dir, exist_ok=True)
        indiv_path = os.path.join(indiv_dir, os.path.basename(save_fits_path))
        psf_list.writeto(indiv_path, overwrite=True)
        print(f"Saved individual PSFs to {indiv_path}")
        
    #except Exception as e:
    #    print(f"Error during PSF generation for file {beam_fits_path}: {e}")
    #    return

    #try:
    # ------------ drizzle all PSFs together ------------
    print("Starting to drizzle PSFs together...")
    output_shape = (psf_shape, psf_shape)  # same as your original choice
    print(f"Output drizzled PSF shape cropped to: {output_shape}")
    target_header, target_wcs = make_wcsheader(
        ra=table_row["ra"],
        dec=table_row["dec"],
        size=output_shape[0] * 0.1,  # 0.1 arcsec/pixel × 30 pixels this is the pixel scale of the CLEAR grism data set
        pixscale=0.1,
        theta=0,  # target alligned to north
    )

    drizzled = np.zeros(output_shape, dtype=np.float32)

    int_times = np.array(int_times, dtype=float)
    if int_times.sum() <= 0:
        weights = np.ones_like(int_times) / len(int_times)
    else:
        weights = int_times / int_times.sum()

    weight_map = np.ones_like(psf_list[1].data, dtype=np.float32)

    for i, (psf_hdu, psf_wcs) in enumerate(zip(psf_list[1:], wcs_list)):
        print(f"Drizzling PSF {i+1}/{len(psf_list) - 1}")

        temp_sci = np.zeros(output_shape, dtype=np.float32)
        temp_wht = np.zeros(output_shape, dtype=np.float32)
        temp_con = np.zeros(output_shape, dtype=np.int32)

        adrizzle.do_driz(
            insci=psf_hdu.data,
            input_wcs=psf_wcs,
            inwht=weight_map,
            output_wcs=target_wcs,
            outsci=temp_sci,
            outwht=temp_wht,
            outcon=temp_con,
            in_units="cps",
            expin=1.0,
            wt_scl=1.0,
            wcslin_pscale=1,
            pixfrac=0.8,
            kernel="square",
            fillval=0.0,
        )

        drizzled += temp_sci * weights[i]

    # save final combined PSF to save_fits_path
    final_hdu = fits.ImageHDU(data=drizzled)
    final_hdu.header.update(target_header)
    final_hdu.writeto(save_fits_path, overwrite=True)
    print(f"Drizzled PSF saved to {save_fits_path}")
    #except Exception as e:
    #    print(f"Error during drizzling process for file {save_fits_path}: {e}")
    #    return

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="Generate drizzled PSF for a given emission line.")
    parser.add_argument("--beam_fits_path", type=str, required=True, help="Path to the beam FITS file.")
    parser.add_argument("--row_fits_path", type=str, required=True, help="Path to the ASCII table row file.")
    parser.add_argument("--wavelength", type=float, required=True, help="Rest-frame wavelength in nm.")
    parser.add_argument("--psf_shape", type=int, default=31, help="Output drizzled PSF shape in pixels.")
    parser.add_argument("--save_fits_path", type=str, required=True, help="Path to save the drizzled PSF FITS file.")
    parser.add_argument("--save_individual_psf", action="store_true", help="Save individual PSFs as a multi-extension FITS.")
    parser.add_argument("--exist_skip", action="store_true", help="Skip generation if output file exists.")
    args = parser.parse_args()
    
    if os.path.exists(args.save_fits_path) and args.exist_skip:
        print(f"File {args.save_fits_path} already exists. Skipping PSF generation.")
    else:
        gen_drizzled_psf_for_wavelength(
            beam_fits_path=args.beam_fits_path,
            row_fits_path=args.row_fits_path,
            wavelength=args.wavelength,
            save_fits_path=args.save_fits_path,
            save_individual_psf=args.save_individual_psf,
        )

if __name__ == "__main__":
    main()