# Balmer-decremen-CLEAR
Study Balmer decrement from CLEAR data

public data available at https://archive.stsci.edu/hlsps/clear/


## package requirement

numpy; astropy; matplotlib; webbpsf; os; tqdm; gc; warnings; pathlib; photutils; magpie(recommend in linux system)

## Current progress

Most primitive function is realized, incl. 

1. loading spectra
2. extract Ha Hb maps, error maps generated from line weight maps
3. pf calculation (ongoing)
4. radial profile of surface brightness & balmer decrement

### progress on psf check

what to expect from psfs: the hbeta psf should be slightly smaller than halpha psf

A sanity check: fit gaussian to each psfs and compare the FWHM; also can compute the growth curve



What can be done in the future:

(+++) SFMS plot on going: theoretical model[Whitaker 2014](https://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdfhttps://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdf) needed

(+++) matching psf kernel sensitive to psfs' high freq noise, window function needs to be tweaked

(+++) the result from balmer decrement radial needs optimization

(+) 1D(not spatially resolved) diagram plot (e.g. balmer decrement, EW against mass SFR)

(+) Ha line maps rescaling to correct for NII line emission [paper](https://iopscience.iop.org/article/10.1088/0004-637X/792/1/75/pdf)

chandry goodsn xray catalog 
>selection.



sfr binning + sfr/m

2d ha/hb

