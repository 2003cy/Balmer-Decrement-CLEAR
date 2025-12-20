from    astropy.io          import fits
from    astropy.wcs         import WCS
import    numpy               as     np


def make_wcsheader(
    ra=40.07293, dec=-1.6137748, size=2, pixscale=0.1, get_hdu=False, theta=0
):
    """
    Make a celestial WCS header

    Parameters
    ----------
    ra, dec : float
        Celestial coordinates in decimal degrees

    size, pixscale : float or 2-list
        Size of the thumbnail, in arcsec, and pixel scale, in arcsec/pixel.
        Output image will have dimensions `(npix,npix)`, where

            >>> npix = size/pixscale

    get_hdu : bool
        Return a `~astropy.io.fits.ImageHDU` rather than header/wcs.

    theta : float
        Position angle of the output thumbnail (degrees)

    Returns
    -------
    hdu : `~astropy.io.fits.ImageHDU`
        HDU with data filled with zeros if `get_hdu=True`.

    header, wcs : `~astropy.io.fits.Header`, `~astropy.wcs.WCS`
        Header and WCS object if `get_hdu=False`.

    Examples
    --------

        >>> from grizli.utils import make_wcsheader
        >>> h, wcs = make_wcsheader()
        >>> print(wcs)
        WCS Keywords
        Number of WCS axes: 2
        CTYPE : 'RA---TAN'  'DEC--TAN'
        CRVAL : 40.072929999999999  -1.6137748000000001
        CRPIX : 10.0  10.0
        CD1_1 CD1_2  : -2.7777777777777e-05  0.0
        CD2_1 CD2_2  : 0.0  2.7777777777777701e-05
        NAXIS    : 20 20

        >>> from grizli.utils import make_wcsheader
        >>> hdu = make_wcsheader(get_hdu=True)
        >>> print(hdu.data.shape)
        (20, 20)
        >>> print(hdu.header.tostring)
        XTENSION= 'IMAGE   '           / Image extension
        BITPIX  =                  -32 / array data type
        NAXIS   =                    2 / number of array dimensions
        PCOUNT  =                    0 / number of parameters
        GCOUNT  =                    1 / number of groups
        CRPIX1  =                   10
        CRPIX2  =                   10
        CRVAL1  =             40.07293
        CRVAL2  =           -1.6137748
        CD1_1   = -2.7777777777777E-05
        CD1_2   =                  0.0
        CD2_1   =                  0.0
        CD2_2   = 2.77777777777777E-05
        NAXIS1  =                   20
        NAXIS2  =                   20
        CTYPE1  = 'RA---TAN'
        CTYPE2  = 'DEC--TAN'
    """

    if np.isscalar(pixscale):
        cdelt = [pixscale / 3600.0] * 2
    else:
        cdelt = [pixscale[0] / 3600.0, pixscale[1] / 3600.0]

    if np.isscalar(size):
        npix = np.asarray(np.round([size/pixscale, size/pixscale]),dtype=int)
    else:
        npix = np.asarray(np.round([size[0]/pixscale, size[1]/pixscale]),dtype=int)

    hout = fits.Header()
    hout["CRPIX1"] = (npix[0] - 1) / 2 + 1
    hout["CRPIX2"] = (npix[1] - 1) / 2 + 1
    hout["CRVAL1"] = ra
    hout["CRVAL2"] = dec
    hout["CD1_1"] = -cdelt[0]
    hout["CD1_2"] = hout["CD2_1"] = 0.0
    hout["CD2_2"] = cdelt[1]
    hout["NAXIS1"] = npix[0]
    hout["NAXIS2"] = npix[1]
    hout["CTYPE1"] = "RA---TAN"
    hout["CTYPE2"] = "DEC--TAN"

    hout["RADESYS"] = "ICRS"
    hout["EQUINOX"] = 2000
    hout["LATPOLE"] = hout["CRVAL2"]
    hout["LONPOLE"] = 180

    hout["PIXASEC"] = pixscale, "Pixel scale in arcsec"

    wcs_out = WCS(hout)

    theta_rad = np.deg2rad(theta)
    mat = np.array(
        [
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad), np.cos(theta_rad)],
        ]
    )

    rot_cd = np.dot(mat, wcs_out.wcs.cd)

    for i in [0, 1]:
        for j in [0, 1]:
            hout["CD{0:d}_{1:d}".format(i + 1, j + 1)] = rot_cd[i, j]
            wcs_out.wcs.cd[i, j] = rot_cd[i, j]

    cd = wcs_out.wcs.cd
    wcs_out.pscale = np.sqrt((cd[0,:]**2).sum())*3600.

    if get_hdu:
        hdu = fits.ImageHDU(
            header=hout, data=np.zeros((npix[1], npix[0]), dtype=np.float32)
        )
        return hdu
    else:
        return hout, wcs_out