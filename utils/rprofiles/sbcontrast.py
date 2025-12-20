# -----------------------------------------------------------------------------
# Overview
# -----------------------------------------------------------------------------

# Authors: Keim, M. A., van Dokkum, P., Li, J.

# Email: michael [dot] keim [at] yale [dot] edu

# Description: The following code is a method to obtain surface brightness
# limits that accurately reflect an images depth at a given spatial scale.
# Rather than relying on the naive Poisson expectation, this code will 
# estimate limitations presented by large scale variations. Here we present an
# example application. We hope to release this as a pip installable code that 
# may be run from command line in the near future. A full description of this
# method is given in Keim et al. 2022. For more information, please refer to: 
# https://www.pietervandokkum.com/software

# Usage: Note that genuine sources should be masked - otherwise the reported
# limit will not be as deep as the true value. Note also that the accuracy of
# this limit is limited by data reduction - for scales exceeding that of 
# background subtraction, genuine features above the calculated limit may 
# have been removed. Moreover, at the single pixel scale sbcontrast may 
# diverge from the true rms due to correlations introduced by re-sampling.

# Notes on this application: By replacing 'image' and 'masks', and providing
# you may run this code to find the Nsigma depth of your image at X scales. 

# -----------------------------------------------------------------------------
# Imports & Dependencies
# -----------------------------------------------------------------------------

import numpy as np
from astropy.stats import biweight_midvariance
from astropy.stats import biweight_location

# -----------------------------------------------------------------------------
# Specifications
# -----------------------------------------------------------------------------

# Image to Calculate Limit On
image = 'data.fits'

# Masks
masks = 'masks.fits'

# Pixel Size in arcseconds
pix   = 0.371

# Zeropoint
zp    = 30.000

# Desired Scale in arcseconds
s     = 60.

# Nsigma Limit
N     = 1. 

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

def sbcontrast(image, mask, pixel_scale, zeropoint, sigma=1.0, scale_arcsec=60, minfrac=0.8, minback=6, verbose=True):
	"""
	A method to calculate the surface brightness detection limit on a given angular scale.

	Parameters:
		image (numpy 2-D array): input image.
		mask (numpy 2-D array): if you want to mask out a pixel, set its value to 1; otherwise set to zero.
		pixel_scale (float): pixel scale of the input image, in the unit of ``arcsec/pixel``.
		zeropoint (float): photometric zeropoint of the input image.
		sigma (float): indicates the detection threshold on which the SB limit will be calculated.
		scale_arcsec (float): on which scale we calculate SB limit, in the unit of ``arcsec``.
			If ``scale_arcsec=60``, this function prints out SB limit on the scale of 60 arcsec * 60 arcsec square.
		minfrac (float): Must be less than 1.0. We discard super-pixels in which less than ``minfrac`` fraction of pixels are available.
			Hence super-pixels with too many pixels masked out are discarded.
		minback (int): Given a super-pixel, we discard it (set to zero) if there are less than ``minback`` non-zero super-pixels surrounding it.
		verbose (bool): whether print out results.
		logger (``logging.logger`` object): logger for this function. Default is ``None``.
	
	"""

	##  Image Dimensions
	ny, nx = image.shape
	scale_pix = scale_arcsec / pixel_scale

	# Image Binning
	scale_x = np.array([scale_pix, int(scale_pix), int(scale_pix), int(scale_pix) + 1])
	scale_y = np.array([scale_pix, int(scale_pix), int(scale_pix) + 1, int(scale_pix) + 1])
	area = scale_x * scale_y
	area -= area[0]
	area = abs(area)[1:] # d1, d2, d3
	bin_x = int(scale_x[np.argmin(area) + 1])
	bin_y = int(scale_y[np.argmin(area) + 1])
	area_ratio = bin_x * bin_y / scale_pix**2
	if verbose:
		print('Binning factors: dx = {0}, dy = {1}'.format(bin_x, bin_y))
		print('Used bin area / True bin area = {:.5f}'.format(area_ratio))
	nbins_x = int(nx / bin_x)
	nbins_y = int(ny / bin_y)

	# Initialize Maps
	im_loc = np.zeros((nbins_y, nbins_x)) # Binned Map
	im_frac = np.zeros((nbins_y, nbins_x)) # Masked Fraction
	im_fluct = np.zeros((nbins_y, nbins_x)) # Contrast/Fluctuation Map

	# Calculate Bin Values
	for i in range(nbins_x - 1):
		for j in range(nbins_y - 1):
			x1, x2, y1, y2 = i * bin_x, (i + 1) * bin_x, j * bin_y, (j + 1) * bin_y
			im_sec = image[y1:y2, x1:x2]
			im_mask_sec = mask[y1:y2, x1:x2]
			im_sec_in = im_sec[(im_mask_sec == 0)]
			if im_sec_in.size > 0:\
				im_loc[j, i] = biweight_location(im_sec_in)
			im_frac[j, i] = 1 - float(im_sec_in.size) / float(im_sec.size)

	# Calculate Local Background and Fluctuation
	for i in range(1, nbins_x - 1):
		for j in range(1, nbins_y - 1):
			backvec = im_loc[j-1:j+2, i-1:i+2]
			backvec = np.delete(backvec.flatten(), 4)
			maskvec = im_frac[j-1:j+2, i-1:i+2]
			maskvec = np.delete(maskvec.flatten(), 4) # fraction of being masked out
			backvec_in = backvec[(maskvec < 1 - minfrac)]
			if len(backvec_in) > minback:
				im_fluct[j,i] = im_loc[j,i] - biweight_location(backvec_in)

	# Ignore pixels with too many masks
	im_fluct_in = im_fluct[im_fluct != 0]

	# Find limit
	sig_adu = np.sqrt(biweight_midvariance(im_fluct_in)) * 1./np.sqrt(1.125)
	sig_adu *= sigma

	# For the standard deviation of standard deviation, see:
	# https://stats.stackexchange.com/questions/631/standard-deviation-of-standard-deviation
	dsig_adu = sig_adu / np.sqrt(2 * (im_fluct_in.size - 1))
	dsig_adu *= sigma

	# Convert limit in ADU to magnitudes
	sb_lim = zeropoint - 2.5 * np.log10(sig_adu / pixel_scale**2)
	dsb_lim = 2.5 * np.log10(1 + 1/np.sqrt(im_fluct_in.size))

	if verbose:
			print('{0:.1f}-sigma variation in counts = {1:.4f} +- {2:.4f}'.format(sigma, sig_adu, dsig_adu))
			print('{0:.1f}-sigma surface brightness limit on {1} arcsec scale is {2:.4f} +- {3:.4f}'.format(sigma, scale_arcsec, sb_lim, dsb_lim))

	return (sb_lim, dsb_lim), (sig_adu, dsig_adu), [im_fluct, im_loc]

# Open Files
from astropy.io import fits
hduli  = fits.open(image)
hdulm  = fits.open(masks)
datai  = hduli[0].data
datam  = hdulm[0].data
hduli.close()
hdulm.close()

# Run
temp = sbcontrast(datai, datam, pix, zp, sigma=N, scale_arcsec=s)




