from   astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
import astropy.units as u
from   astropy.io import fits
from   astropy.table import Table
from   matplotlib import cm
import numpy as np
from   photutils.aperture import EllipticalAperture
from   astropy.visualization import LogStretch, ImageNormalize
from   utils.colors import color1, color2, color3, color_notation
import matplotlib
matplotlib.use('Agg')
from   matplotlib import pyplot as plt

# Unified paper plot style (aligned with analysis.ipynb)
plt.rcParams.update({
    'font.size'         : 14,
    'axes.labelsize'    : 16,
    'axes.titlesize'    : 16,
    'xtick.labelsize'   : 13,
    'ytick.labelsize'   : 13,
    'legend.fontsize'   : 14,
    'xtick.direction'   : 'in',
    'ytick.direction'   : 'in',
    'xtick.top'         : True,
    'ytick.right'       : True,
    'xtick.major.size'  : 8,
    'ytick.major.size'  : 8,
    'xtick.minor.size'  : 6,
    'ytick.minor.size'  : 6,
    'xtick.major.width' : 1.8,
    'ytick.major.width' : 1.8,
    'xtick.minor.width' : 1,
    'ytick.minor.width' : 1,
    'axes.linewidth'    : 1.2,
    'savefig.dpi'       : 300,
    'savefig.bbox'      : 'tight',
})

PROFILE_TICK_KWARGS = dict(
    direction='in',
    which='both',
    top=True,
    right=True,
    left=True,
    bottom=True,
    labelsize=13,
    width=1.8,
    length=8,
)

IMAGE_TICK_KWARGS = dict(PROFILE_TICK_KWARGS, labelbottom=False, labelleft=False)
ERRORBAR_KWARGS = dict(
    fmt='-',
    linewidth=3,
    marker='o',
    markersize=4,
    capsize=8,
    capthick=2,
    elinewidth=2,
)


def style_image_axis(ax, name):
    ax.tick_params(axis='both', **IMAGE_TICK_KWARGS)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.text(
        0.03, 0.95, f"{name}", transform=ax.transAxes,
        fontsize=22, color=color_notation, va='top', ha='left',
        bbox=dict(facecolor='black', alpha=0.35, edgecolor='none', pad=3)
    )


def style_profile_axis(ax, title, xlabel='Distance (kpc)'):
    ax.tick_params(axis='both', **PROFILE_TICK_KWARGS)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=13, loc='best')


def plot_for_obj(extracted_fits_path, profile_fits_path, save_plot_path):
    with fits.open(profile_fits_path) as hdul:
        obj = Table(hdul[1].data)[0]
        r=obj['distance']/0.1 * obj['pixel_length'] #in kpc
        ha=obj['ha']; ha_err=obj['ha_err']; ha_lim=obj['ha_limit']
        hb=obj['hb']; hb_err=obj['hb_err']; hb_lim=obj['hb_limit']
        balmer = obj['balmer']; balmer_err = obj['balmer_err']
        xlim = max(10,2*obj['re']/0.1 * obj['pixel_length']) #in kpc
        mask_detected = (ha>np.sqrt(ha_err**2 + ha_lim**2)) & (hb>np.sqrt(hb_err**2 + hb_lim**2))
        
    with fits.open(extracted_fits_path) as hdul:
        
        re = obj['re']/0.1; q = obj['q']; pa = obj['pa'] #re in pixels
        
        #prepare segmentation map
        seg = hdul[2].data == obj['ID']
        #prepare 2D Balmer map
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=300)
        for plot_index, hdul_index, name in [
                (0, 4, r'H$\alpha$'),
                (1, 6, r'H$\beta$'),
                (2, 3, f"{hdul[3].header['FILTER']}")]:
            ax = axes.flat[plot_index]
            style_image_axis(ax, name)
            data = hdul[hdul_index].data
            # avoid zeros and negatives for log stretch
            norm = ImageNormalize(np.where(data > 0, data, np.nan), interval=PercentileInterval(99.5), stretch=LogStretch())

            # plot the segmentation map
            cmap = plt.get_cmap('plasma')
            grey_cmap = cm.colors.ListedColormap(
                np.mean(cmap(np.linspace(0, 1, 256))[:, :3], axis=1, keepdims=True).repeat(3, axis=1))
            grey_cmap._init()

            ax.imshow(np.where(seg, data, np.nan), origin='lower', norm=norm, cmap=cmap)
            ax.imshow(np.where(~seg, data, np.nan), origin='lower', norm=norm, cmap=grey_cmap, alpha=0.7)
            # set the title

            # plot the ellipse at the effective radius for all images, use photutils aperture
            center = ( (data.shape[1]-2)/2, (data.shape[0]-2)/2)
            aperture = EllipticalAperture(center, a=re, b=re * q, theta=np.deg2rad(pa*u.deg - 90*u.deg))
            aperture.plot(ax=ax, color='gray', lw=2.5, alpha=0.85, label=r'$r_{eff}$')
            #ax.scatter(center[1], center[0], color='white', s=90, marker='x', linewidths=2, label='Center')
            
            # Set the axis limits to center a 30x30 region around the image center
            y, x = center
            ax.set_xlim(x - 20, x + 20)
            ax.set_ylim(y - 20, y + 20)

            # physical scale bar anchored to the bottom-left of the visible region
            #set a physical scale bar for 5 kpc
            scale_bar_length_kpc = 5
            scale_bar_length_pix = scale_bar_length_kpc / obj['pixel_length'] #/ 0.1

            bar_x0 = (x - 20) + 3
            bar_y0 = (y - 20) + 3
            ax.plot([bar_x0, bar_x0 + scale_bar_length_pix], [bar_y0, bar_y0],
                    color=color_notation, linewidth=4, solid_capstyle='butt')
            ax.text(bar_x0 + scale_bar_length_pix / 2, bar_y0 + 1.2,
                    f'{scale_bar_length_kpc} kpc', color=color_notation, fontsize=23,
                    ha='center', va='bottom')
            
#--------------plot radial profiles-----------------
        #plot the radial profile of ha and ha s/n
        ax = axes[1, 0]
        #plot a grey area representing the effective radius
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3, label=r'$r_{eff}$')
        ax.errorbar(r, ha, yerr=ha_err, label=r'H$\alpha$', color=color2, **ERRORBAR_KWARGS)
        #plot ha limit as horizontal line
        ax.plot(r, ha_lim, color=color2, linestyle='-.', alpha=0.6, linewidth=2.5, label=r'H$\alpha$ limit')
        #plot detected ha
        ax.scatter(r[mask_detected], ha[mask_detected], facecolors='none', edgecolors='black', s=55, linewidths=1.5, label='Detected',zorder=10)
        
        # shared y-limits for Ha and Hb so both panels use the same ytick labels
        sb_vals = np.concatenate([ha + ha_err, hb + hb_err, ha_lim, hb_lim])
        sb_vals = sb_vals[np.isfinite(sb_vals) & (sb_vals > 0)]
        sb_ylim = (np.nanmin(sb_vals) * 0.7, np.nanmax(sb_vals) * 1.5)

        ax.set_xlim(-0.2, xlim)
        ax.set_yscale('log')
        ax.set_ylim(*sb_ylim)
        ax.set_ylabel('[erg/s/cm$^2$/arcsec$^2$]')
        style_profile_axis(ax, r'H$\alpha$ Surface Brightness')

        #plot the radial profile of hb and hb s/n
        ax = axes[1, 1]
        #plot a grey area representing the effective radius
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3)#, label=r'$r_{eff}$')
        ax.errorbar(r, hb, yerr=hb_err, label=r'H$\beta$', color=color1, **ERRORBAR_KWARGS)
        #plot hblimit as horizontal line
        ax.plot(r, hb_lim, color=color1, linestyle='-.', alpha=0.6, linewidth=2.5, label=r'H$\beta$ limit')
        #plot detected hb
        ax.scatter(r[mask_detected], hb[mask_detected], facecolors='none', edgecolors='black', s=55, linewidths=1.5, zorder=10)#, label='Detected')
        ax.set_xlim(-0.2, xlim)
        ax.set_yscale('log')
        ax.set_ylim(*sb_ylim)
        #ax.set_ylabel('[erg/s/cm$^2$/arcsec$^2$]')
        style_profile_axis(ax, r'H$\beta$ Surface Brightness')

        #plot the radial profile of balmer decrement
        ax = axes[1, 2]
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3)#, label=r'$r_{eff}$')
        ax.errorbar(
            r, balmer, yerr=balmer_err, label=r'H$\alpha$/H$\beta$',
            color=color3, **ERRORBAR_KWARGS
        )
        #plot used points
        ax.errorbar(
            r[mask_detected], balmer[mask_detected], yerr=balmer_err[mask_detected],
            fmt='o', color='black', markersize=6, elinewidth=2, capsize=8, capthick=2,
            label=r'> 1-$\sigma$', zorder=10
        )
        ax.axhline(y=2.86, color='black', linestyle='--', linewidth=2.5, alpha=0.8, label='2.86')
        ax.fill_between(r, 2.86*0.8, 2.86*1.2, color=color1, alpha=0.18)
        ax.set_xlim(-0.2, xlim)
        ax.set_ylim(0, np.nanmin([np.nanmax(balmer)+3,25]))
        style_profile_axis(ax, r'Balmer Decrement H$\alpha$/H$\beta$')
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.22, hspace=0.25)
        plt.savefig(save_plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

def cli():
    import argparse
    parser = argparse.ArgumentParser(description="Plot diagnostic plots for a given object.")
    parser.add_argument('--extracted_fits_path', type=str, help='Path to the extracted fits file.')
    parser.add_argument('--profile_fits_path', type=str, help='Path to the profile fits file.')
    parser.add_argument('--save_plot_path', type=str, help='Path to save the diagnostic plot.')
    args = parser.parse_args()
    
    plot_for_obj(args.extracted_fits_path, args.profile_fits_path, args.save_plot_path)
    print(f"Diagnostic plot saved to {args.save_plot_path}")
    
if __name__ == "__main__":
    cli()
