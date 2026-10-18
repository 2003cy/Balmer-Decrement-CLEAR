# Measuring the Spatially Resolved Balmer Decrement of Galaxies in the CLEAR Survey

Code and analysis pipeline for measuring the spatially resolved Balmer decrement
(Hα/Hβ) — and hence dust attenuation — as a function of galactocentric radius for
galaxies in the [CLEAR](https://archive.stsci.edu/hlsp/clear) (CANDELS Lyman-alpha
Emission at Reionization) HST/WFC3 grism survey. This supports an ongoing draft paper.

The end-to-end workflow (catalog → download → PSF → line maps → radial profiles →
diagnostics) is orchestrated with [Snakemake](https://snakemake.readthedocs.io/);
final science plots are produced in [`analysis.ipynb`](analysis.ipynb).

<p align="center">
<img src="diagnosticplot.png" alt="Per-galaxy diagnostic plot" width="1000">
</p>

*Example per-galaxy diagnostic: (top) PSF-matched Hα, Hβ and F105W maps; (bottom)
Hα and Hβ surface-brightness profiles with detection limits, and the resulting
Balmer decrement vs. radius (dashed line = Case B value 2.86).*

---

## 1. Setup — install with pixi

The environment is fully specified in [`pixi.toml`](pixi.toml) / `pixi.lock`
(Python, numpy, astropy, photutils, matplotlib, drizzlepac, snakemake, …).

```bash
# install pixi once (https://pixi.sh)
curl -fsSL https://pixi.sh/install.sh | bash

# from the repo root — installs the locked environment
pixi install
```

**One external dependency:** PSF generation uses [TinyTim](https://www.stsci.edu/hst/instrumentation/focus-and-pointing/tiny-tim).
Point the `TINYTIM` environment variable at your TinyTim install (a copy lives in
`../tinytim-7.5`), otherwise `3_psf.py` falls back to a hard-coded default path:

```bash
export TINYTIM=/path/to/tinytim-7.5
```

Before running, edit [`snakemake/config.yaml`](snakemake/config.yaml) and set
`home_dir` to your machine (all outputs are written under `home_dir/output_dir`).

---

## 2. End-to-end run — data → analysis

The whole pipeline is a single Snakemake DAG. Run it with the pixi tasks:

```bash
pixi run run       # full pipeline, 1 core (good for a first/test run)
pixi run run-all   # same pipeline, 72 cores (per-galaxy steps run in parallel)
```

On a cluster, submit [`slurm`](slurm) (SLURM script; runs Snakemake twice to settle
the download checkpoint, then collects results).

### Pipeline stages

| Step | Script | What it does | Output |
|------|--------|--------------|--------|
| 1 | `1_catalogs_and_initial_select.py` | Download CLEAR/CANDELS catalogs, clean & make initial selection | `catalogs/master_cat.fits` |
| 2 | `2_download_spectra.py` | rsync the `full` / `beams` grism data products per object *(checkpoint)* | `data_products/` |
| 3 | `3_psf.py` + `3_psf_match.py` | Build TinyTim PSFs at Hα & Hβ wavelengths, then a matching kernel | `data_psf/` |
| 4 | `4_line_map_extraction.py` | Extract Hα & Hβ line maps, NII-correct (Zahid+14), PSF-match Hβ | `data_extracted/` |
| 5 | `5_radial_profiles.py` | Elliptical-annulus photometry → radial SB & Balmer-decrement profiles | `data_radial_profiles/` |
| 6 | `6_diag_plot.py` | Per-galaxy diagnostic plot (the figure above) | `plots/diagnostic_plot_all/` |
| 7 | `7_collect_results.py` | Stack all per-object profiles into one master table | `master_table.fits` |

### Final analysis

Open [`analysis.ipynb`](analysis.ipynb) (register the kernel first with
`pixi run register-kernel`). It reads the collected `master_table` and produces the
sample-level science figures:

- Sample properties (mass / redshift / selection)
- Balmer decrement vs. physical radius (kpc)
- Dust attenuation profiles in normalized radius (r/rₑ)
- Hα surface brightness → SFR surface density profiles

Reusable helpers live in [`utils/`](utils/) (attenuation curves, NII corrections,
radial-profile apertures, surface-brightness limits, catalog I/O).

---

## Repository layout

```
1_… 7_*.py          numbered pipeline stages (run via Snakemake)
Snakefile            workflow DAG
snakemake/config.yaml  paths, threads, per-step parameters
pixi.toml / pixi.lock  reproducible environment + task shortcuts
analysis.ipynb       sample-level science plots
diagnosticplot.png   example per-galaxy diagnostic
utils/               shared analysis code
output/              generated catalogs, data products, profiles, plots
```

---

## Diagnostic & suggested plots

`diagnosticplot.png` in the repo root is the reference per-galaxy diagnostic
(re-generated for every object in `output/plots/diagnostic_plot_all/`). It is the
first thing to check when validating a run — it shows whether the Hα/Hβ maps, PSF
matching, and profile extraction all behaved for a given galaxy.

Additional plots worth adding for QA and for the paper:

- **Sample selection plot** — mass–redshift (or UVJ) diagram with the parent
  catalog vs. the final selection highlighted, to make the selection reproducible.
- **PSF-matching residual** — Hβ PSF convolved with the kernel minus the Hα PSF
  (and encircled-energy curves), to quantify match quality across the sample.
- **Stacked / median Balmer-decrement profile** — all galaxies overlaid plus a
  mass-binned median with scatter, the key population-level result.
- **Detection-limit / completeness plot** — fraction of annuli above the Hβ limit
  vs. radius, so the outer-profile reliability is explicit.
- **AGN / contamination check** — Chandra X-ray cross-match
  (`output/plots/chandra_agn_crossmatch.png` already exists) shown alongside the
  sample to flag AGN hosts.
```
