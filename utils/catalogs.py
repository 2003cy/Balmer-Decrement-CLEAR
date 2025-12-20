from astropy.table import Table
from pathlib import Path
import tarfile
import os

def download_catalogs(master_catalog_path='catalogs/master_cat_raw.fits'):
    
    def download_file_from_url(url, save_name,overwrite = False):
        import os
        from urllib.request import urlretrieve

        os.makedirs(os.path.dirname(master_catalog_path), exist_ok=True)
        filename = f"{os.path.dirname(master_catalog_path)}/{save_name}"
        if not os.path.exists(filename) or overwrite:
            print(f"Downloading {filename}...")
            urlretrieve(url, filename)
        else:
            print(f"{filename} already exists. Skipping download.")
    
    catalogs = {
        # spectroscopic catalogs
        "gn_s": "https://archive.stsci.edu/hlsps/clear/catalogs/gdn/spectroscopic_catalog/hlsp_clear_hst_wfc3_gdn_multi_v4.1_clear.fits",
        "gs_s": "https://archive.stsci.edu/hlsps/clear/catalogs/gds/spectroscopic_catalog/hlsp_clear_hst_wfc3_gds_multi_v4.1_clear.fits",
        
        # EAZY catalogs
        "gn_e": "https://archive.stsci.edu/hlsps/clear/catalogs/gdn/eazy/hlsp_clear_hst_wfc3-acs_gdn-3dhst_multi_v4.6_zout.fits",
        
        "gs_e": "https://archive.stsci.edu/hlsps/clear/catalogs/gds/eazy/hlsp_clear_hst_wfc3-acs_gds-3dhst_multi_v4.6_zout.fits",
        
        # morphology catalog
        "morph": "https://users.ugent.be/~avdrwel/data/allfields.tar",
        
        #chandra catalogs
        "chandra_gn": "https://content.cld.iop.org/journals/0067-0049/224/2/15/revision1/apjs523032t3_mrt.txt",
        "chandra_gs": "https://content.cld.iop.org/journals/0067-0049/228/1/2/revision1/apjsaa4dbdt4_mrt.txt"
    }

    for catalog_name, url in catalogs.items():
        #for candra catalog, renew the file name for better readability
        if "chandra" in catalog_name:
            save_name = f"{catalog_name}.txt"
        else:
            save_name = url.split('/')[-1]
        

        if os.path.exists(f"{os.path.dirname(master_catalog_path)}/{save_name}"):
            print(f"{save_name} already exists. Skipping download.")
        else:
            download_file_from_url(url, save_name, overwrite=False)
            
            
            
            
def load_catalogs(catalog_dir):
        gn = Table.read(f'{catalog_dir}/hlsp_clear_hst_wfc3_gdn_multi_v4.1_clear.fits')
        gs = Table.read(f'{catalog_dir}/hlsp_clear_hst_wfc3_gds_multi_v4.1_clear.fits')

        gn_m = Table.read(f'{catalog_dir}/hlsp_clear_hst_wfc3-acs_gdn-3dhst_multi_v4.6_zout.fits')
        gs_m = Table.read(f'{catalog_dir}/hlsp_clear_hst_wfc3-acs_gds-3dhst_multi_v4.6_zout.fits')
        
        
        # unzip morphology catalog into catalog_dir and load galfit tables
        morph_dir = Path(f'{catalog_dir}/goodsn')
        if not morph_dir.exists():
                with tarfile.open(f'{catalog_dir}/allfields.tar', 'r') as tar:
                    tar.extractall(path=catalog_dir)

        gn_gal = Table.read(f"{catalog_dir}/goodsn/goodsn_3dhst.v4.1_f125w.galfit", format='ascii')
        gs_gal = Table.read(f"{catalog_dir}/goodss/goodss_3dhst.v4.1_f125w.galfit", format='ascii')
        
        return gn, gs, gn_m, gs_m, gn_gal, gs_gal