#! /usr/bin/env python

# Import packages
import os
import argparse
import subprocess
from functools import partial
import numpy as np
from astropy.table import Table
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# Download product with rsync
def download_spectrum(extract, output_dir, data_products, extract_rows):
    """Download product from remote server."""

    # Product field + id
    field = extract['subfield'].lower()
    id    = str(extract['ID']).zfill(5)
    
    if extract_rows:
        row_path = f"{output_dir}/{field}_{id}_row.fits"
        if os.path.exists(row_path):
            print(f'{field}-{id} rows exists')
        else:
            row = Table(extract)
            row.write(row_path, format='fits', overwrite=True)
            print(f'{field}-{id} rows saved')
    
    # Files
    files = [];remote_urls = []; file_name_for_save = []
    for product in data_products:
        # Remote URL
        remote_urls.append(f"https://archive.stsci.edu/hlsps/clear/data/{product}/{field}/")
        files.append(f"hlsp_clear_hst_wfc3_{field}-{id}_g102-g141_v4_{product}.fits")
        file_name_for_save.append(f"{field}_{id}_{product}.fits")

    # Execute command
    for file,remote_url, save_name in zip(files,remote_urls, file_name_for_save):
        # Download command
        out_path = f'{output_dir}/{save_name}'
        command = [
            'curl',
            '-f',
            '--output',
            out_path,
            f'{remote_url}/{file}',
        ]

        if os.path.exists(out_path):
            print(f'{file} exists')
            continue

        try:
            subprocess.run(command, check=True)
            print(f'\n download {file} downloaded \n')
        except subprocess.CalledProcessError as e:
            print(f'Failed to download {file}. Error: {e}')


# Main Function
def main():

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('extracted', type=str, help='Path to extracted table')
    parser.add_argument('--ncpu', type=int, default=1, help='Number of download threads')
    parser.add_argument('--output-dir', default='data_products', help='Directory for downloaded spectra')
    parser.add_argument('--data-products', nargs='+', default=['full', 'beams'], help='Data product types to download (e.g. 1d full stack beams)')
    parser.add_argument('--extract-rows', type=bool, default=True, help='Whether to extract individual rows')
    args = parser.parse_args()

    # Load extracted
    extracted = Table.read(args.extracted)

    # Number of CPUs
    ncpu = args.ncpu
    
    # Create directories
    os.makedirs(args.output_dir, exist_ok=True)


# Multi-threaded download
    if ncpu > 1:
        with ThreadPoolExecutor(ncpu) as executor:
            list(tqdm(
                executor.map(partial(download_spectrum, output_dir=args.output_dir, data_products=args.data_products, extract_rows=args.extract_rows), extracted),
                total=len(extracted),
                desc="Downloading spectra"
            ))

    # Single-threaded download
    else:
        for extract in tqdm(extracted, desc="Downloading spectra"):
            download_spectrum(extract, args.output_dir, args.data_products, args.extract_rows)


if __name__ == '__main__':
    main()
