from astropy.io import fits
from astropy.table import Table
from tqdm.auto import tqdm
import numpy as np
import os

def collect_results():
    profiles_dir = 'output/data_radial_profiles/'
    fits_files = [os.path.join(profiles_dir, f) for f in os.listdir(profiles_dir) if f.endswith('.fits')]
    rows = []

    for file_path in fits_files:
        with fits.open(file_path) as hdul:
            obj_row = Table(hdul[1].data)[0]
            rows.append(obj_row)
    master_table = Table(rows=rows)

    has_profile_lis = []
    for obj in tqdm(master_table):
            r = obj['distance'] #arcsec
            ha = obj['ha']; ha_limit=obj['ha_limit']
            hb = obj['hb']; hb_limit=obj['hb_limit']
            #first manual selection based on sn_limit
            mask_r = (r < 1*obj['re']) #within 1 re
            mask_sb_limit = (ha > ha_limit) & (hb > hb_limit)
            if len(r[mask_r & mask_sb_limit]) >= 2:
                    has_profile = 1
            else:
                    has_profile = 0
            has_profile_lis.append(has_profile)

    master_table['has_profile'] = has_profile_lis
    print('initial selection based on radial 2 sigma sb limit',len(master_table[master_table['has_profile'] == 1]))
    master_table.write('master_table.fits', overwrite=True)
    return master_table


def copy_diagnostic_plots(master_table):
    src_dir   = 'output/plots/diagnostic_plot_all'
    trget_dir = 'output/plots/diagnostic_plot_selected'
    if os.path.exists(trget_dir):
        #remove existing directory
        os.system(f"rm -r {trget_dir}")
    os.makedirs(trget_dir, exist_ok=True)

    for obj in tqdm(master_table[master_table['has_profile'] == 1]):
            #copy diagnostic plots to a separate folder for easier manual inspection
            src_plot = f"{src_dir}/{obj['subfield'].lower()}_{str(obj['ID']).zfill(5)}.png"
            trget_plot = f"{trget_dir}/{obj['subfield'].lower()}_{str(obj['ID']).zfill(5)}.png"
            os.system(f"cp {src_plot} {trget_plot}")

def main():
    print("Collecting results from radial profile fitting...")
    master_table = collect_results()
    print("Copying diagnostic plots for selected objects...")
    copy_diagnostic_plots(master_table)
    
if __name__ == "__main__":
    main()