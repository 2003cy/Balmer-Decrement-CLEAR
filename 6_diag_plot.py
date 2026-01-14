from   astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
import astropy.units as u
from   astropy.io import fits
from   astropy.table import Table
from   matplotlib import cm
import numpy as np
from   photutils.aperture import EllipticalAperture
from   astropy.visualization import LogStretch, ImageNormalize
from   utils.colors import color1, color2, color3
import matplotlib
from   matplotlib import pyplot as plt
matplotlib.use('Agg')

def plot_for_obj(extracted_fits_path, profile_fits_path, save_plot_path):
    with fits.open(profile_fits_path) as hdul:
        obj = Table(hdul[1].data)[0]
        r=obj['distance']/0.1 * obj['pixel_length'] #in kpc
        ha=obj['ha']; ha_err=obj['ha_err']; ha_lim=obj['ha_limit']
        hb=obj['hb']; hb_err=obj['hb_err']; hb_lim=obj['hb_limit']
        balmer = obj['balmer']; balmer_err = obj['balmer_err']
        
    with fits.open(extracted_fits_path) as hdul:
        
        re = obj['re']/0.1; q = obj['q']; pa = obj['pa'] #re in pixels
        
        #prepare segmentation map
        seg = hdul[2].data == obj['ID']
        #prepare 2D Balmer map
        
        plt.figure(figsize=(15,10))        
        for plot_index, hdul_index, name in [
                (1, 3, f"{hdul[3].header['FILTER']}"),
                (2, 4, 'Hα'),
                (3, 6, 'Hβ')]:
            ax = plt.subplot(2, 3, plot_index)
            ax.tick_params(axis='both', labelbottom=False, direction='in', which='both', 
            top=True, right=True, left=True, bottom=True, labelsize=13)
            
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            
            ax.text(
                0.02, 0.98, f"{name}", transform=ax.transAxes,
                fontsize=18, color='white', va='top', ha='left',
                bbox=dict(facecolor='black', alpha=0.3, edgecolor='none', pad=2)
            )
            data = hdul[hdul_index].data
            # avoid zeros and negatives for log stretch
            norm = ImageNormalize(np.where(data > 0, data, np.nan), interval=PercentileInterval(99.5), stretch=LogStretch())

            # plot the segmentation map
            viridis_cmap = plt.get_cmap('viridis')
            grey_viridis_cmap = cm.colors.ListedColormap(
                np.mean(viridis_cmap(np.linspace(0, 1, 256))[:, :3], axis=1, keepdims=True).repeat(3, axis=1))
            grey_viridis_cmap._init()

            ax.imshow(np.where(seg, data, np.nan), origin='lower', norm=norm, cmap='viridis')
            ax.imshow(np.where(~seg, data, np.nan), origin='lower', norm=norm, cmap=grey_viridis_cmap, alpha=0.7)
            # set the title


            # plot the ellipse at the effective radius for all images, use photutils aperture
            center = ( (data.shape[1]-2)/2, (data.shape[0]-2)/2)
            aperture = EllipticalAperture(center, a=re, b=re * q, theta=np.deg2rad(pa*u.deg - 90*u.deg))
            aperture.plot(ax=ax, color='gray', lw=2, alpha=0.8, label='Effective Radius')
            ax.scatter(center[1], center[0], color='white', s=100, marker='x', label='Center')

            # Set the axis limits to center a 30x30 region around the image center
            y, x = center
            ax.set_xlim(x - 20, x + 20)
            ax.set_ylim(y - 20, y + 20)
            
#--------------plot radial profiles-----------------
        ax = plt.subplot(2, 3, 4)
        ax.tick_params(axis='both', labelbottom=True, direction='in', which='both', 
                top=True, right=True, left=True, bottom=True, labelsize=13)
        #plot a grey area representing the effective radius
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3, label='Effective Radius')
        ax.errorbar(r, ha, yerr=ha_err, label='Ha', fmt='o:', color=color2, markersize=2,
            capsize=4, capthick=1.5, elinewidth=1.5)
        #plot halimit as horizontal line
        ax.plot(r, ha_lim, color='darkred', linestyle='-.', label='Ha Limit (alt)')

        ax.set_xlim(-0.2, 10)
        ax.set_xlabel('Distance (kpc)')
        ax.set_title('$H\\alpha$ Surface Brightness')
        ax.set_yscale('log')
        ax.legend()

        #plot the radial profile of ha and hb s/n
        ax = plt.subplot(2, 3, 5)
        ax.tick_params(axis='both', labelbottom=True, direction='in', which='both', 
                top=True, right=True, left=True, bottom=True, labelsize=13)
        #plot a grey area representing the effective radius
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3, label='Effective Radius')
        ax.errorbar(r, hb, yerr=hb_err, label='Hb', fmt='o:', color=color1, markersize=2,
            capsize=4, capthick=1.5, elinewidth=1.5)
        #plot hblimit as horizontal line
        ax.plot(r, hb_lim, color='darkgreen', linestyle='-.', label='Hb Limit (alt)')
                                        
        ax.set_xlim(-0.2, 10)
        ax.set_xlabel('Distance (kpc)')
        ax.set_title('$H\\beta$ Surface Brightness')
        ax.set_yscale('log')
        ax.legend()

        #plot the radial profile of balmer decrement
        ax = plt.subplot(2, 3, 6)
        ax.tick_params(axis='both', labelbottom=True, direction='in', which='both', 
                top=True, right=True, left=True, bottom=True, labelsize=13)
        ax.axvspan(0, re*obj['pixel_length'], color='grey', alpha=0.3, label='Effective Radius')
        ax.errorbar(
            r, balmer, yerr=balmer_err, label='Balmer Decrement',
            fmt='o:', color=color3, markersize=4,
            capsize=4, capthick=1.5, elinewidth=1.5
        )
        ax.axhline(y=2.86, color='green', linestyle='--', label='theory')
        ax.fill_between(r, 2.86*0.8, 2.86*1.2, color='green', alpha=0.2)
        ax.set_xlim(-0.2, 8)
        ax.set_ylim(0, np.nanmin([np.nanmax(balmer)+3,25]))
        ax.set_xlabel('Distance [kpc]')
        ax.set_title('Balmer Decrement Ha/Hb')
        plt.tight_layout()
        plt.savefig(save_plot_path, dpi=300)
        plt.show()

def cli():
    import argparse
    parser = argparse.ArgumentParser(description="Plot diagnostic plots for a given object.")
    parser.add_argument('--extracted_fits_path', type=str, help='Path to the extracted fits file.')
    parser.add_argument('--profile_fits_path', type=str, help='Path to the profile fits file.')
    parser.add_argument('--save_plot_path', type=str, help='Path to save the diagnostic plot.')
    args = parser.parse_args()
    
    plot_for_obj(args.extracted_fits_path, args.profile_fits_path, args.save_plot_path)
    
if __name__ == "__main__":
    cli()