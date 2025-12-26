from   astropy.table import Table
from   astropy.io import fits
import numpy as np
import argparse
import os

def profile_auto_selection(obj):
    
    r = obj['distance']
    ha = obj['ha']
    ha_err = obj['ha_err']
    ha_lim = obj['ha_limit']
    hb = obj['hb']
    hb_err = obj['hb_err']
    hb_lim = obj['hb_limit']
    re = obj['re']
    
    mask_r = (r<1*re) #within 1 re
    mask_limit= (ha>ha_lim) & (hb>hb_lim)  #both lines,spatially detected

    if len(r[mask_r & mask_limit]) > 0:
        return True
    else:
        return False
        
        
def main():
    parser = argparse.ArgumentParser(description="Auto select radial profiles based on detection criteria")
    parser.add_argument('--all_radial_profiles', nargs='+', required=True, help='List of all radial profile FITS files')
    parser.add_argument('--selected_catalog', type=str, required=True, help='Path to save the selected catalog')
    args = parser.parse_args()
    
    #first open and contatenate all radial profiles
    selected_objs = []
    for profile_fits_path in args.all_radial_profiles:
        try:
            table = Table.read(profile_fits_path)
            obj = table[0]
            if profile_auto_selection(obj):
                selected_objs.append(obj)
                print(f"Object {obj['field']}_{obj['id']} passed the selection criteria.")
            else:
                print(f"Object {obj['field']}_{obj['id']} did not pass the selection criteria.")
        except Exception as e:
            print(f"Error reading radial profile {profile_fits_path}: {e}")
            continue
    if len(selected_objs) > 0:
        selected_table = Table(rows=selected_objs)
        selected_table.write(args.selected_catalog, overwrite=True)
        print(f"Saved selected catalog to {args.selected_catalog}")
    else:
        print("No objects passed the selection criteria. No catalog saved.")