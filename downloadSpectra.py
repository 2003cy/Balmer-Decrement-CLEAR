#! /usr/bin/env python

# Import packages
import os
import argparse
import subprocess
import numpy as np
from astropy.table import Table
from concurrent.futures import ThreadPoolExecutor


# Download product with rsync
def download_spectrum(extract):
    """Download product from remote server."""

    # Product field + id
    field = extract['subfield'].lower()
    id    = str(extract['ID']).zfill(5)
    # Files
    files =[]
    remote_urls = []
    for type in ['full','stack','beams']:
        # Remote URL
        remote_urls.append(f"https://archive.stsci.edu/hlsps/clear/data/{type}/{field}/")
        files.append(f"hlsp_clear_hst_wfc3_{field}-{id}_g102-g141_v4_{type}.fits")

    # Execute command
    for file,remote_url in zip(files,remote_urls):
        # Download command
        command = [
            'curl',
            '-f',
            '--output',
            f'data_products/{file}',
            f'{remote_url}/{file}',
        ]

        if os.path.exists(f'data_products/{file}'):
            print(f'{file}exist')
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

    parser.add_argument('--ncpu', type=int, default=1)
    args = parser.parse_args()

    # Load extracted
    extracted = Table.read(args.extracted)

    # Number of CPUs
    ncpu = args.ncpu

    # Create directories
    home = os.getcwd()
    os.makedirs(os.path.join(home,'data_products'), exist_ok=True)


    # Multi-threaded download
    if ncpu > 1:
        with ThreadPoolExecutor(ncpu) as executor:
            executor.map(
                lambda e: download_spectrum(e),
                extracted,
            )

    # Single-threaded download
    else:
        for extract in extracted:
            download_spectrum(extract)


if __name__ == '__main__':
    main()
