from   astropy.table        import Table, vstack, join
from   astropy.coordinates  import SkyCoord, match_coordinates_sky
import astropy.units        as     u
from   astropy.cosmology    import Planck18
from   utils.chandra        import load_chandra_gn, load_chandra_gs
from   utils.catalogs       import download_catalogs, load_catalogs
import os
import numpy                as     np
import matplotlib.pyplot    as     plt
import argparse


# This is a script to:
#  download the necessary catalogs for object selection from the CLEAR survey.
#  perform catalog cleaning and initial object selection.
#  the downloading of inidividual catalogs is handled in utils/catalogs.py

def generate_master_catalog(reload=False, master_catalog_path='catalogs/master_cat_raw.fits'):
    catalog_dir = os.path.dirname(master_catalog_path)
        
    if not reload and os.path.exists(master_catalog_path):
        print('Master catalog already exists. Skipping generation and loading.')
        return Table.read(master_catalog_path)
    else:
        
        print('Generating master catalog... \n Downloading necessary catalogs if not present.')
        download_catalogs(master_catalog_path=master_catalog_path)

        # Load catalogs: gn gs spectroscopic, gn gs eazy, gn gs morphology
        gn, gs, gn_m, gs_m, gn_gal, gs_gal = load_catalogs(catalog_dir)
        
        ### CLEAR GoodsN+S sample prearation, stacking
        #connect morphology, eazy, spectroscopic catalogs
        gn_all = join(gn_m,gn_gal,keys_left='id',keys_right='NUMBER',join_type='left')
        gn_cat = join(gn,gn_all,keys_left='ID',keys_right='id',join_type='left',metadata_conflicts='silent')
        gn_cat['field'] = 'gn'

        gs_all = join(gs_m,gs_gal,keys_left='id',keys_right='NUMBER',join_type='left')
        gs_cat = join(gs,gs_all,keys_left='ID',keys_right='id',join_type='left',metadata_conflicts='silent')
        gs_cat['field'] = 'gs'
        
        #initial catalog without cleaning
        master_cat = vstack([gn_cat,gs_cat],metadata_conflicts='silent')
        master_cat.write(master_catalog_path,overwrite=True)

'''
The following function cleans the master catalog by applying several quality cuts:
- non-negative redshift
- positive stellar mass
- realistic effective radius (f <= 1)
- positive H-alpha and H-beta fluxes
It also calculates the signal-to-noise ratios.
'''

def clean_catalogs(master_catalog_path='catalogs/master_cat_raw.fits',
                    cleaned_catalog_path='catalogs/master_cat.fits',
                    plt_output_dir=None,
                    chandra_match_radius = 2,
                    test_length=-1):
    
    catalog_dir = os.path.dirname(master_catalog_path)
    cat_lis = Table.read(master_catalog_path)

    #constrain positive redshift
    z_real = cat_lis['z_MAP'] > 0
    #constrain positive mass
    mass_real = cat_lis['mass'] > 0
    #constrain positive re
    re_real = cat_lis['f'] <=  1
    #constrain positie ha and hb
    ha_real = cat_lis['Ha_FLUX'] > 0
    hb_real = cat_lis['Hb_FLUX'] > 0
    
    radius_real = cat_lis['re'] > 0
    good_morph = cat_lis['f'] <= 1

    print('total number of objects:', len(cat_lis))
    print('fraction of bad mass:', len(cat_lis[cat_lis['mass'] <= 0]) / len(cat_lis))
    print('fraction of bad re:', len(cat_lis[cat_lis['f'] > 1]) / len(cat_lis))
    print('fraction of bad z:', len(cat_lis[cat_lis['z_MAP'] <= 0]) / len(cat_lis))
    print('fraction of bad values before bad ha hb:', len(cat_lis[~(z_real & mass_real & re_real)]) / len(cat_lis))
    print('fraction of bad ha:', len(cat_lis[~ha_real]) / len(cat_lis))
    print('fraction of bad hb:', len(cat_lis[~hb_real]) / len(cat_lis))
    print('fraction of bad radius:', len(cat_lis[~radius_real]) / len(cat_lis))
    print('fraction of bad morphology:', len(cat_lis[~good_morph]) / len(cat_lis))
    cat_lis = cat_lis[z_real & mass_real & re_real & ha_real & hb_real & radius_real & good_morph]
    
    #calculate sn ratio
    sn_ha = cat_lis['Ha_FLUX']/cat_lis['Ha_FLUX_ERR']
    cat_lis['sn_ha'] = sn_ha
    sn_hb = cat_lis['Hb_FLUX']/cat_lis['Hb_FLUX_ERR']
    cat_lis['sn_hb'] = sn_hb
    
    #calculate pixel_length; given pixel scale is 0.1 arcsec/pixel
    cat_lis['pixel_length']  =  np.deg2rad(0.1/3600) * Planck18.angular_diameter_distance(cat_lis['z_MAP']).to(u.kpc).value
    
    print('number of objects after cleaning:', len(cat_lis))

    #--- AGN selection-----
    # find agns based on xray catalogs
    chandra_gn = load_chandra_gn(catalog_dir)
    chandra_gs = load_chandra_gs(catalog_dir)
    
    #sample objects
    coord_obj  = SkyCoord(cat_lis['ra'],cat_lis['dec'],unit=(u.deg,u.deg))

    #
    #goodss agn crossmatch
    agn_gs      = chandra_gs[chandra_gs['OType']=='AGN']
    coord_chan_gs = SkyCoord(agn_gs['RAdeg'],agn_gs['DEdeg'],unit=(u.deg, u.deg))
    idx_gs,d2d_gs,d3d_gs = match_coordinates_sky(coord_chan_gs,coord_obj)

    #goodsn agn catalog
    agn_gn      = chandra_gn[chandra_gn['Type']=='AGN']
    #agn_mask_gn = chandra_gn['Type'] != 'ddf'
    coord_chan_gn = SkyCoord(agn_gn['RAdeg'],agn_gn['DEdeg'],unit=(u.deg, u.deg))
    idx_gn,d2d_gn,d3d_gn = match_coordinates_sky(coord_chan_gn,coord_obj)

    # --- Cross-matching and tagging ---
    thresh = chandra_match_radius * u.arcsec
    mask_gs = d2d_gs < thresh; mask_gn = d2d_gn < thresh
    # Combine indices of all matched sources (within threshold) from both fields
    id_matched = np.append(idx_gs[mask_gs], idx_gn[mask_gn])

    # Tag sources as 'agn' if matched, otherwise 'gxy'
    cat_lis['tag'] = ['agn' if index in id_matched else 'gxy' for index in range(len(cat_lis))]
    # Save the full catalog with AGN tags
    cat_lis[:int(test_length)].write(cleaned_catalog_path, overwrite=True)
    print(f"Cleaned catalog saved to {cleaned_catalog_path}, with length {len(cat_lis[:int(test_length)])}.")


    # ------------------- visualization of agn cross-matching ------------------------
    if plt_output_dir is not None:
        os.makedirs(plt_output_dir, exist_ok=True)
        plt.figure(figsize=(20, 10))
        # Make first subplot: GOODS-South
        plt.subplot(121)
        # Plot all GOODS-South sources (blue, small)
        gs_cat = cat_lis[cat_lis['field'] == 'gs']
        plt.scatter(gs_cat['ra'], gs_cat['dec'], s=1, alpha=0.5, color='blue', label='GOODS-South, Ha/Hb flux > 0')
        # Plot cross-matched Chandra AGN (red 'x')
        plt.scatter(
            agn_gs['RAdeg'][mask_gs],
            agn_gs['DEdeg'][mask_gs],
            s=100, alpha=1, color='red', label='GOODS-South Chandra catalog',
            marker='x', edgecolor='orange'
        )
        plt.xlim(np.min(gs_cat['ra']), np.max(gs_cat['ra']))
        plt.ylim(np.min(gs_cat['dec']), np.max(gs_cat['dec']))
        plt.legend(loc='upper right', fontsize=12)
        
        # Make second subplot: GOODS-North
        plt.subplot(122)
        # Plot all GOODS-North sources (blue, small)
        gn_cat = cat_lis[cat_lis['field'] == 'gn']
        plt.scatter(gn_cat['ra'], gn_cat['dec'], s=1, alpha=0.5, color='blue', label='GOODS-North, Ha/Hb flux > 0')
        # Plot cross-matched Chandra AGN (red 'x', orange edge)
        plt.scatter(
            agn_gn['RAdeg'][mask_gn],
            agn_gn['DEdeg'][mask_gn],
            s=100, alpha=1, color='red', label='GOODS-North Chandra catalog',
            marker='x', edgecolor='orange'
        )
        plt.xlim(np.min(gn_cat['ra']), np.max(gn_cat['ra']))
        plt.ylim(np.min(gn_cat['dec']), np.max(gn_cat['dec']))
        plt.legend(loc='upper right', fontsize=12)
        plt.suptitle('Cross-matching CLEAR catalog with Chandra X-ray AGN catalogs', fontsize=16)
        plt.savefig(f'{plt_output_dir}/chandra_agn_crossmatch.png', dpi=300)
        plt.close()
    
    return cat_lis



def main(master_catalog_raw='catalogs/master_cat_raw.fits',
                 master_catalog_clean='catalogs/master_cat.fits',
                 reload_master=False,
                 plt_output_dir='plots/',
                 test_length=-1):
    """Run download, generate raw catalog, and clean catalog in one call."""

    os.makedirs(os.path.dirname(master_catalog_raw) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(master_catalog_clean) or '.', exist_ok=True)

    generate_master_catalog(
        reload=reload_master,
        master_catalog_path=master_catalog_raw,
    )

    clean_catalogs(
        master_catalog_path=master_catalog_raw,
        cleaned_catalog_path=master_catalog_clean,
        plt_output_dir=plt_output_dir,
        test_length=test_length,
    )
    
    


def cli():
    parser = argparse.ArgumentParser(description="Generate and clean CLEAR master catalog")
    parser.add_argument("--raw", dest="master_catalog_raw", default='catalogs/master_cat_raw.fits',
                        help="Path to write the raw master catalog FITS file")
    parser.add_argument("--clean", dest="master_catalog_clean", default='catalogs/master_cat.fits',
                        help="Path to write the cleaned master catalog FITS file")
    parser.add_argument("--plt_dir", dest="plt_output_dir", default='plots/',
                        help="Directory to save output plots")
    parser.add_argument("--reload", dest="reload_master", action="store_true",
                        help="Force re-generation of the raw master catalog even if it exists")
    parser.add_argument("--test_length", dest="test_length", type=int, default=-1,
                        help="Number of objects to use for a test run; set to -1 for full catalog")
    args = parser.parse_args()

    main(
        master_catalog_raw=args.master_catalog_raw,
        master_catalog_clean=args.master_catalog_clean,
        reload_master=args.reload_master,
        plt_output_dir=args.plt_output_dir,
        test_length=args.test_length,
    )


if __name__ == "__main__":
    cli()