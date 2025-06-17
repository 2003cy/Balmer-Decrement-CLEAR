""" @file make_psf.py

    Created 5 Jul 2019

    This module contains the needed functions to generate a Tiny Tim PSF, using either
    the best-fit optical parameters found by Gillis et al. (2019) or the linear
    relationship with the focus-secondary-mirror despace.

    ---------------------------------------------------------------------

    Copyright (C) 2019 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program.  If not, see
    <http://www.gnu.org/licenses/>.
"""

__all__ = ['make_subsampled_model_psf']

import os
import subprocess as sbp


# Default values
default_psf_position = (2048, 1024)  # Center of the detector by default
default_focus = -1.0  # Approximately the middle of expected values
default_chip = 1
default_spec_file_path = '' # Use spectrum for a K-type star by default must be specified
default_filter_name = 'F105W'
default_detector = 23  # WFC3/IR
default_psf_size = 3.0
default_tinytim_path = "/Users/yangcheng/tinytim"  # Adjust as needed for your own purposes
default_subsampling_factor = 8
default_mono = False

# Default optical parameters
optical_params_means = {"z2": -0.0042,
                        "z3":                 0.0046,
                        "astigmatism_0":      0.0241,
                        "astigmatism_45":     0.0300,
                        "coma_x":             0.0159,
                        "coma_y":             0.0000,
                        "clover_x":           0.0074,
                        "clover_y":           0.0163,
                        "spherical_3rd": -0.0217,
                        "z12":                0.0037,
                        "z13":                0.0001,
                        "z14":                0.0043,
                        "z15":                0.0059,
                        "z16": -0.0061,
                        "z17":                0.0059,
                        "z18":                0.0039,
                        "z19":                0.0020,
                        "z20": -0.0008,
                        "z21":                0.0072,
                        "spherical_5th":      0.0101,
                        "kernel_adjustment":  0.9978}

# Linear fit for optical parameters in (intercept, slope)
optical_params_int_and_slopes = {"z2":                (-0.0043,  0.0000),
                                 "z3":                (0.0048, -0.0002),
                                 "astigmatism_0":     (0.0246, -0.0005),
                                 "astigmatism_45":    (0.0302, -0.0002),
                                 "coma_x":            (0.0172, -0.0012),
                                 "coma_y":            (0.0000,  0.0000),
                                 "clover_x":          (0.0079, -0.0004),
                                 "clover_y":          (0.0169, -0.0005),
                                 "spherical_3rd":     (-0.0207, -0.0008),
                                 "z12":               (0.0002,  0.0005),
                                 "z13":               (0.0003, -0.0002),
                                 "z14":               (0.0039,  0.0003),
                                 "z15":               (0.0061, -0.0002),
                                 "z16":               (-0.0069,  0.0007),
                                 "z17":               (0.0059,  0.0000),
                                 "z18":               (0.0034, -0.0004),
                                 "z19":               (0.0026, -0.0005),
                                 "z20":               (-0.0014,  0.0005),
                                 "z21":               (0.0062,  0.0008),
                                 "spherical_5th":     (0.0088,  0.0011),
                                 "kernel_adjustment": (0.9978,  0.0000), }


def replace_multiple_in_file(input_filename, output_filename, input_strings, output_strings):
    """ Replaces every occurence of an input_string in input_filename with the corresponding
        output string and prints the results to $output_filename.

        @param[in]  input_filename <str>
        @param[out] output_filename <str>
        @param[in]  input_strings <iterable of strs>
        @param[in]  output_strings <iterable of strs>

        @return None
    """

    with open(output_filename, "w") as fout:
        with open(input_filename, "r") as fin:
            for line in fin:
                new_line = line
                for input_string, output_string in zip(input_strings, output_strings):
                    if((input_string is None) or (output_string is None)):
                        continue
                    new_line = new_line.replace(input_string, output_string)
                fout.write(new_line)

    return


def make_subsampled_model_psf(filename,
                              psf_position=default_psf_position,
                              focus=default_focus,
                              spec_file_path=default_spec_file_path,
                              detector=default_detector,
                              filter_name=default_filter_name,
                              psf_size=default_psf_size,
                              tinytim_path=default_tinytim_path,
                              subsampling_factor=default_subsampling_factor,
                              linear_fit=False,
                              exist_skip=True,
                              mono = default_mono,
                              **optical_params):
    if psf_position[0]>1023 or psf_position[1]>1023:
        return f'{filename} has a problem with psf position'
    """ Generates a subsampled model PSF, using the desired (or default) optical parameters.
        For input parameters spec_type and detector, the allowed options can be seen through
        running tiny1

        @param[out] filename <str> Desired filename of the generated PSF. If it already exists, the
                                   'clobber' parameter will determine whether or not it will be
                                   overwritten.
        @param[in]  psf_position <(float, float)> Position on the detector of the PSF in x, y
        @param[in]  focus <float> Focus-secondary-mirror despace of the PSF
        @param[in]  chip <int> Which chip the model PSF is for. Allowed values are 1 and 2
        @param[in]  spec_type <(int, *)> Spectral type of the PSF to generate. First value chooses type
                                         of spectrum, second chooses from options for this type
        @param[in]  detector <int> Index of detector to be used.
        @param[in]  filter_name <str> Name of the filter to use (eg. F606W)
        @param[in]  psf_size <float> Size of the PSF image in arcseconds
        @param[in]  tinytim_path <str> Location of the Tiny Tim executables
        @param[in]  subsampling_factor <int> Factor by which to subsample the PSF
        @param[in]  linear_fit <bool> If False, unspecified optical parameters will be given values
                                      based on the mean from Gillis et al. (2018)'s testing. If True,
                                      will use the linear fit from the analysis instead
        @param[in]  clobber <bool> Whether or not to overwrite the target file if it already exists.
        @param[in]  optical_params <dict> Optical parameters aside from focus for this PSF. If not
                                      specified here, defaults will be used based on the
                                      linear_fit parameter.

        @return None


    """
    os.environ["TINYTIM"] = tinytim_path


    if exist_skip:
            if os.path.exists(filename):
                return f"{filename} skipped"
    # Create a directory to contain this project
    try:
        os.makedirs(os.path.split(filename)[0])
    except OSError as e:
        if not ("[Errno 17] File exists:" in str(e) or "[Errno 2] No such file or directory: ''" in str(e)):
            raise
        else:
            pass  # No need to raise if the directory already exists

    filename_base = filename.replace(".fits", "")
    if not filename_base + ".fits" == filename:
        raise ValueError("Filename (" + filename + ") must end in '.fits'.")

    par_file = filename_base + ".par"
    tmp_par_file = filename_base + ".par.tmp"

    # Set up the command to call tiny1 and execute it
    if mono == False:
        cmd = f"{tinytim_path}/tiny1" + f' {tmp_par_file}' + " << EOF \n" + \
                str(detector) + "\n" + \
                str(psf_position[0]) + " " + str(psf_position[1]) + "\n" + \
                str(filter_name) + "\n" + \
                "5 \n" + str(spec_file_path) + "\n" + \
                str(psf_size) + "\n" + \
                str(focus) + "\n" + \
                filename_base + "\nEOF"
    else:
        cmd = f"{tinytim_path}/tiny1" + f' {tmp_par_file}' + " << EOF \n" + \
                str(detector) + "\n" + \
                str(psf_position[0]) + " " + str(psf_position[1]) + "\n" + \
                'MONO ' + "\n" + \
                str(1000) + "\n" + \
                str(psf_size) + "\n" + \
                str(focus) + "\n" + \
                filename_base + "\nEOF"

    #sbp.run(cmd, shell=True,executable="/bin/bash")
    sbp.run(cmd, shell=True,executable="/bin/bash",
            stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
            #timeout=10,check=True)

#============================
    # Determine which optical parameters we'll be using
    optical_params_to_use = {}
    for param in optical_params_means:
        if param in optical_params:
            optical_params_to_use[param] = optical_params[param]
        elif linear_fit:
            intercept, slope = optical_params_int_and_slopes[param]
            optical_params_to_use[param] = intercept + focus * slope
        else:
            optical_params_to_use[param] = optical_params_means[param]


    # Edit the parameter file to adjust optical parameters
    strs_to_replace = []
    replacements = []

    strs_to_replace.append("0.       # Z2 = X (V2) tilt")
    replacements.append(str(optical_params_to_use["z2"]) + "       # Z2 = X (V2) tilt")

    strs_to_replace.append("0.       # Z3 = Y (V3) tilt")
    replacements.append(str(optical_params_to_use["z3"]) + "       # Z3 = Y (V3) tilt")

    strs_to_replace.append("0.031    # Z5 = 0 degree astigmatism")
    replacements.append(str(optical_params_to_use["astigmatism_0"]) + "    # Z5 = 0 degree astigmatism")

    strs_to_replace.append("0.028    # Z6 = 45 degree astigmatism")
    replacements.append(str(optical_params_to_use["astigmatism_45"]) + "    # Z6 = 45 degree astigmatism")

    strs_to_replace.append("0.003    # Z7 = X (V2) coma")
    replacements.append(str(optical_params_to_use["coma_x"]) + "    # Z7 = X (V2) coma")

    strs_to_replace.append("0.001    # Z8 = Y (V3) coma")
    replacements.append(str(optical_params_to_use["coma_y"]) + "    # Z8 = Y (V3) coma")

    strs_to_replace.append("0.007    # Z9 = X clover")
    replacements.append(str(optical_params_to_use["clover_x"]) + "    # Z9 = X clover")

    strs_to_replace.append("0.018    # Z10 = Y clover")
    replacements.append(str(optical_params_to_use["clover_y"]) + "    # Z10 = Y clover")

    strs_to_replace.append("-0.025    # Z11 = 3rd order spherical")
    replacements.append(str(optical_params_to_use["spherical_3rd"]) + "    # Z11 = 3rd order spherical")

    strs_to_replace.append("0.       # Z12 = 0 degree Spherical astigmatism")
    replacements.append(str(optical_params_to_use["z12"]) + "       # Z12 = 0 degree Spherical astigmatism")

    strs_to_replace.append("0.       # Z13 = 45 degree Spherical astigmatism")
    replacements.append(str(optical_params_to_use["z13"]) + "       # Z13 = 45 degree Spherical astigmatism")

    strs_to_replace.append("0.       # Z14 = X Ashtray")
    replacements.append(str(optical_params_to_use["z14"]) + "       # Z14 = X Ashtray")

    strs_to_replace.append("0.       # Z15 = Y Ashtray")
    replacements.append(str(optical_params_to_use["z15"]) + "       # Z15 = Y Ashtray")

    strs_to_replace.append("0.       # Z16")
    replacements.append(str(optical_params_to_use["z16"]) + "       # Z16")

    strs_to_replace.append("0.       # Z17")
    replacements.append(str(optical_params_to_use["z17"]) + "       # Z17")

    strs_to_replace.append("0.       # Z18")
    replacements.append(str(optical_params_to_use["z18"]) + "       # Z18")

    strs_to_replace.append("0.       # Z19")
    replacements.append(str(optical_params_to_use["z19"]) + "       # Z19")

    strs_to_replace.append("0.       # Z20")
    replacements.append(str(optical_params_to_use["z20"]) + "       # Z20")

    strs_to_replace.append("0.       # Z21")
    replacements.append(str(optical_params_to_use["z21"]) + "       # Z21")

    strs_to_replace.append("0.009    # Z22 = 5th order spherical")
    replacements.append(str(optical_params_to_use["spherical_5th"]) + "    # Z22 = 5th order spherical")

    replace_multiple_in_file(tmp_par_file, par_file, strs_to_replace, replacements)
    #------------------------------------------------

    # this part is very important!!!
    #very interesting that tinytim1 param file generation has a weird limitation in wavelength, thats narrower than WFC3 response curve
    #so this part we manually change the wavelength in the par file
    with open(par_file,'r') as file:
        lines = file.readlines()
    parts = lines[18].split()
    parts[0] = str(round(mono/1000,6)) #change the wavlength (in microns)
    lines[18] = ' '.join(parts) + '\n'
    with open(par_file,'w') as output_file:
        output_file.writelines(lines)
    #------------------------------------------------

    # Set up the command to call tiny2
    cmd = f"{tinytim_path}/tiny2 " + par_file
    # Run the command to call tiny2
    #sbp.run(cmd, shell=True,executable="/bin/bash")
    sbp.run(cmd, shell=True,executable="/bin/bash",
            stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
            #timeout=10,check=True)
    


    # Set up the command to call tiny3
    cmd = f"{tinytim_path}/tiny3 " + par_file  + " SUB=" + str(int(subsampling_factor))
    # Run the command to call tiny3
    #sbp.run(cmd, shell=True,executable="/bin/bash")
    sbp.run(cmd, shell=True,executable="/bin/bash",
            stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)
            #timeout=10,check=True)




    # PSF should be generated, now move it to the desired filename
    init_filename = filename_base + "00.fits"
    os.rename(init_filename, filename)

    # Clean up unneeded files. Silently suppress any errors here
    try:
        os.remove(filename_base + "00_psf.fits")
    except OSError as _e:
        pass
    try:
        os.remove(filename_base + ".tt3")
    except OSError as _e:
        pass
    try:
        os.remove(init_filename)
    except OSError as _e:
        pass
    try:
        os.remove(par_file)
    except OSError as _e:
        pass
    try:
        os.remove(tmp_par_file)
    except OSError as _e:
        pass

    return f"{filename_base} processed"


    
