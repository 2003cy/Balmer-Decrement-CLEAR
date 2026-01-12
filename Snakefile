num_threads = config["num_threads"]
conda_env_flag = config["env_flag"]
home_dir = config["home_dir"]
output_dir = config["output_dir"]
import os
import glob
run_dir = os.path.join(home_dir, output_dir)

path_lis = config["paths"]
data_products_dir = os.path.join(run_dir, path_lis["data_products_dir"])
data_extracted_dir = os.path.join(run_dir, path_lis["line_extraction_dir"])
radial_profiles_dir = os.path.join(run_dir, path_lis["radial_profiles_dir"])
psf_dir = os.path.join(run_dir, path_lis["psf_dir"])
plt_dir  = os.path.join(run_dir, path_lis["plt_dir"])

#helper functrion to initialize every individual galaxy data 
#gen data field and id lists from master catalog
from snakemake.io import glob_wildcards, expand
def _ensure_download_checkpoint(wildcards=None):
    # Force Snakemake to wait for the download checkpoint before globbing results.
    return checkpoints.download_data.get(**wildcards) if wildcards else checkpoints.download_data.get()
def list_downloaded_objs(wildcards=None):
    _ensure_download_checkpoint(wildcards)
    pattern = f"{data_products_dir}/{{field}}_{{id}}_full.fits"
    return glob_wildcards(pattern)

#---------------the following are output file lists based on downloaded objects----------------------------
def psf_ha(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{psf_dir}/{{field}}_{{id}}_psf_ha.fits",
        zip,field=wc.field, id=wc.id,)

def psf_hb(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{psf_dir}/{{field}}_{{id}}_psf_hb.fits",
        zip,field=wc.field, id=wc.id,)

def psf_kernel(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{psf_dir}/{{field}}_{{id}}_psf_kernel.fits",
        zip,field=wc.field, id=wc.id,)

def extracted_objs(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{data_extracted_dir}/{{field}}_{{id}}_extracted.fits",
        zip,field=wc.field, id=wc.id,)

def radial_profiles(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{radial_profiles_dir}/{{field}}_{{id}}_profile.fits",
        zip,field=wc.field, id=wc.id,)

def diagnostic_plots(wildcards=None):
    wc = list_downloaded_objs(wildcards)
    return expand(f"{plt_dir}/diagnostic_plot_all/{{field}}_{{id}}.png",
        zip,field=wc.field, id=wc.id,)

# top-level rule
rule all:
    input:
        f"{run_dir}/{path_lis['master_catalog_raw']}",
        f"{run_dir}/{path_lis['master_catalog_clean']}",
        f"{data_products_dir}/.download_completed",
        psf_ha,
        psf_hb,
        psf_kernel,
        extracted_objs,  
        radial_profiles,
        diagnostic_plots,
        
#generate master catalog
cfg = config["gen_master_catalog"]
rule master_catalog:
    threads: 1
    input:
        script="1_catalogs_and_initial_select.py",
    output:
        raw=f"{run_dir}/{path_lis['master_catalog_raw']}",
        clean=f"{run_dir}/{path_lis['master_catalog_clean']}",
    log:
        f"{run_dir}/logs/master_catalog.log",
    params:
        conda_env_flag=conda_env_flag,
        reload_flag=cfg["reload"],
        test_length=cfg["test_length"],
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --raw {output.raw} \
            --clean {output.clean} \
            --plt_dir {plt_dir} \
            --test_length {params.test_length} \
            2>&1 | tee {log}
        """

#download data based on master catalog
cfg_dl = config["download_data"]
checkpoint download_data:
    threads: num_threads
    input:
        script="2_download_spectra.py",
        catalog_to_download=f"{run_dir}/{path_lis['master_catalog_clean']}",
    output:
        flag=f"{data_products_dir}/.download_completed",
    log:
        f"{run_dir}/logs/download_data.log",
    params:
        conda_env_flag=conda_env_flag,
        output_dir=lambda wildcards: f"{data_products_dir}",
        data_products=lambda wildcards: " ".join(cfg_dl.get("data_products", [])),
        extract_rows=lambda wildcards: cfg_dl.get("extract_individual_rows", False),
    shell:
        """
        mkdir -p $(dirname {log}) {params.output_dir}
        conda run {params.conda_env_flag} python {input.script} \
            {input.catalog_to_download} \
            --ncpu {threads} \
            --output-dir {params.output_dir} \
            --data-products {params.data_products} \
            --extract-rows {params.extract_rows} \
            2>&1 | tee {log}
        touch {output.flag}
        """

cfg_psfha = config["psf_ha"]
rule gen_psf_ha:
    input:
        script="3_psf.py",
        beam_fits=f"{data_products_dir}/{{field}}_{{id}}_beams.fits",
        row_fits=f"{data_products_dir}/{{field}}_{{id}}_row.fits",
    output:
        psf_ha_fits=f"{psf_dir}/{{field}}_{{id}}_psf_ha.fits",
    log:
        f"{run_dir}/logs/gen_psf_ha/{{field}}_{{id}}.log",
    params:
        conda_env_flag=conda_env_flag,
        wavelength=cfg_psfha["wavelength"],
        psf_shape=cfg_psfha.get("psf_shape", 31),
        save_individual=lambda wildcards: "--save_individual_psf" if cfg_psfha.get("save_individual_psf", False) else "",
        exist_skip= lambda wildcards: "--exist_skip" if cfg_psfha["exist_skip"] else "",
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --beam_fits_path {input.beam_fits} \
            --row_fits_path {input.row_fits} \
            --save_fits_path {output.psf_ha_fits} \
            --wavelength {params.wavelength} \
            --psf_shape {params.psf_shape} \
            {params.exist_skip} \
            {params.save_individual} \
            2>&1 | tee {log}
        """

cfg_psfhb = config["psf_hb"]
rule gen_psf_hb:
    input:
        script="3_psf.py",
        beam_fits=f"{data_products_dir}/{{field}}_{{id}}_beams.fits",
        row_fits=f"{data_products_dir}/{{field}}_{{id}}_row.fits",
    output:
        psf_hb_fits=f"{psf_dir}/{{field}}_{{id}}_psf_hb.fits",
    log:
        f"{run_dir}/logs/gen_psf_hb/{{field}}_{{id}}.log",
    params:
        conda_env_flag=conda_env_flag,
        wavelength=cfg_psfhb["wavelength"],
        psf_shape=cfg_psfhb.get("psf_shape", 31),
        save_individual=lambda wildcards: "--save_individual_psf" if cfg_psfhb.get("save_individual_psf", False) else "",
        exist_skip="--exist_skip" if cfg_psfhb["exist_skip"] else "",
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --beam_fits_path {input.beam_fits} \
            --row_fits_path {input.row_fits} \
            --save_fits_path {output.psf_hb_fits} \
            --wavelength {params.wavelength} \
            --psf_shape {params.psf_shape} \
            {params.exist_skip} \
            {params.save_individual} \
            2>&1 | tee {log}
        """


cfg_psfm= config["psf_match"]
rule psf_match:
    input:
        script="3_psf_match.py",
        psf_ha=f"{psf_dir}/{{field}}_{{id}}_psf_ha.fits",
        psf_hb=f"{psf_dir}/{{field}}_{{id}}_psf_hb.fits",
    output:
        kernel_fits=f"{psf_dir}/{{field}}_{{id}}_psf_kernel.fits",
    log:
        f"{run_dir}/logs/psf_match/{{field}}_{{id}}.log",
    params:
        conda_env_flag=conda_env_flag,
        window=cfg_psfm.get("window", "None"),
        alpha=cfg_psfm.get("alpha", 3.0),
        beta=cfg_psfm.get("beta", 0.9),
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --psf_ha_path {input.psf_ha} \
            --psf_hb_path {input.psf_hb} \
            --save_kernel_path {output.kernel_fits} \
            --window {params.window} \
            --alpha {params.alpha} \
            --beta {params.beta} \
            2>&1 | tee {log}
        """



#extract lines for inidividual objects
cfg_el = config["line_extraction"]
rule extract_lines:
    input:
        script="4_line_map_extraction.py",
        full_fits=f"{data_products_dir}/{{field}}_{{id}}_full.fits",
        row_fits=f"{data_products_dir}/{{field}}_{{id}}_row.fits",
        kerenel_fits=f"{psf_dir}/{{field}}_{{id}}_psf_kernel.fits",
    output:
        extracted_fits=f"{data_extracted_dir}/{{field}}_{{id}}_extracted.fits",
    log:
        f"{run_dir}/logs/extract_lines/{{field}}_{{id}}.log",
    params:
        conda_env_flag=conda_env_flag,
        exist_skip=cfg_el.get("exist_skip", False),
        center_crop=cfg_el.get("window_size", 50),
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --full_fits_path {input.full_fits} \
            --table_row_path {input.row_fits} \
            --kernel_fits_path {input.kerenel_fits} \
            --center_crop {params.center_crop} \
            --exist_skip {params.exist_skip} \
            --line_fits_path {output.extracted_fits} \
            2>&1 | tee {log} 
        """

cfg_rp = config["radial_profiles"]
rule gen_radial_profiles:
    input:
        script="5_radial_profiles.py",
        extracted_fits=f"{data_extracted_dir}/{{field}}_{{id}}_extracted.fits",
        row_fits=f"{data_products_dir}/{{field}}_{{id}}_row.fits",
    output:
        profile_fits = f"{radial_profiles_dir}/{{field}}_{{id}}_profile.fits",
    params:
        conda_env_flag=conda_env_flag,
        annuli_width = cfg_rp.get("annuli_width", 1),
    log:
        f"{run_dir}/logs/radial_profiles/{{field}}_{{id}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        conda run {params.conda_env_flag} python {input.script} \
            --extracted_fits_path {input.extracted_fits} \
            --row_fits_path {input.row_fits} \
            --profile_fits_path {output.profile_fits} \
            --annuli_width {params.annuli_width} \
            2>&1 | tee {log}
        """

#diagnostic plots for individual objects
rule diagnostic_plots:
    input:
        script="6_diag_plot.py",
        extracted_fits=f"{data_extracted_dir}/{{field}}_{{id}}_extracted.fits",
        profile_fits=f"{radial_profiles_dir}/{{field}}_{{id}}_profile.fits",
    output:
        plot_png=f"{plt_dir}/diagnostic_plot_all/{{field}}_{{id}}.png",
    log:
        f"{run_dir}/logs/diagnostic_plots/{{field}}_{{id}}.log",
    params:
        conda_env_flag=conda_env_flag,
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.plot_png})
        conda run {params.conda_env_flag} python {input.script} \
            --extracted_fits_path {input.extracted_fits} \
            --profile_fits_path {input.profile_fits} \
            --save_plot_path {output.plot_png} \
            2>&1 | tee {log}
        """