from   astropy.io          import fits
import numpy               as     np
from   photutils.psf       import matching as match
from   astropy.convolution  import convolve_fft 
import argparse
import inspect




#generate matching kernel for psf_ha & psf_hb
def gen_kernel(psf_hb,psf_ha,window=match.CosineBellWindow,alpha=3,beta=0.9):
        
        if window==None:
                kernel = fits.ImageHDU(
                data = match.create_matching_kernel(psf_hb.data,psf_ha.data),
                name = 'PSF_MATCH')
                kernel.data = kernel.data/np.sum(kernel.data)
                return kernel 
        
        if len(inspect.signature(window).parameters) == 0:
            window = window()
        elif len(inspect.signature(window).parameters) == 1:
            window = window(alpha)
        elif len(inspect.signature(window).parameters) == 2:
            window = window(alpha,beta)

        kernel = fits.ImageHDU(
                data = match.create_matching_kernel(psf_hb.data,psf_ha.data,window=window),
                name = 'PSF_MATCH')
        kernel.data = kernel.data/np.sum(kernel.data)
        
        return kernel 
    
    
def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="Generate PSF matching kernel between H-alpha and H-beta.")
    parser.add_argument("--psf_ha_path", type=str, required=True, help="Path to the H-alpha PSF FITS file.")
    parser.add_argument("--psf_hb_path", type=str, required=True, help="Path to the H-beta PSF FITS file.")
    parser.add_argument("--save_kernel_path", type=str, required=True, help="Path to save the PSF matching kernel FITS file.")
    parser.add_argument("--window", type=str, default="CosineBellWindow", help="Type of window function to use.")
    parser.add_argument("--alpha", type=float, default=3.0, help="Alpha parameter for the window function.")
    parser.add_argument("--beta", type=float, default=0.9, help="Beta parameter for the window function.")
    args = parser.parse_args()
    
    try:
        # load PSFs
        psf_ha = fits.open(args.psf_ha_path)[1]
        psf_hb = fits.open(args.psf_hb_path)[1]
        
        #window function from string
        window_dict = {
            "CosineBellWindow": match.CosineBellWindow,
            "HanningWindow": match.HanningWindow,
            "SplitCosineBellWindow": match.SplitCosineBellWindow,
            "TopHatWindow": match.TopHatWindow,
            "TukeyWindow": match.TukeyWindow,
            "CosineBellWindow": match.CosineBellWindow,
            "None": None
        }
        window_func = window_dict.get(args.window, match.CosineBellWindow)
        
        # generate matching kernel
        kernel_hdu = gen_kernel(psf_hb, psf_ha, window=window_func, alpha=args.alpha, beta=args.beta)
        
        # save kernel to file
        kernel_hdu.writeto(args.save_kernel_path, overwrite=True)
        print(f"PSF matching kernel saved to {args.save_kernel_path}")
        
    except Exception as e:
        print(f"Error in generating PSF matching kernel for {args.psf_ha_path} and {args.psf_hb_path}: {e}")
        #generate a placeholder emptyfile, with error message
        with open(args.save_kernel_path, 'w') as f:
            f.write(f"Error in generating PSF matching kernel: {e}\n")

if __name__ == "__main__":
    main()