from astropy.table import Table
import argparse
import os


def extract_rows_from_table(path_to_table, output_dir):
    """
    Extract individual rows from a table and save them as separate FITS files.

    Parameters
    ----------
    path_to_table : str
        Path to the input table file.
    output_dir : str
        Directory where the extracted row files will be saved.
    """
    # Read the full table
    table = Table.read(path_to_table)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop over each row in the table
    for extract in table:
        field = extract['subfield'].lower()
        id    = str(extract['ID']).zfill(5)
        row_path = f"{output_dir}/{field}_{id}_row.fits"
        row = Table(extract)
        row.write(row_path, format='fits', overwrite=True)
        print(f'Extracted row saved to {row_path}')\
            
def main():
    parser = argparse.ArgumentParser(description="Extract individual rows from a table.")
    parser.add_argument('table_path', type=str, help='Path to the input table file')
    parser.add_argument('--output-dir', default='extracted_rows', help='Directory to save extracted row files')
    args = parser.parse_args()

    extract_rows_from_table(args.table_path, args.output_dir)