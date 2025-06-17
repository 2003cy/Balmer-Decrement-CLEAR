#this is just a handy little function to return the desired file path
#give one entry in the object list, return the desired file path
def file_name(obj,prefix,filetype='fits'):
    field = obj['subfield'].lower()
    id    = str(obj['ID']).zfill(5)
    return f"hlsp_clear_hst_wfc3_{field}-{id}_g102-g141_v4_{prefix}.{filetype}"


#save updataed or add hduimage to the extracted file
def save_update(image_to_save,extracted):
        for i,image in enumerate(extracted):
            if image.name == image_to_save.name:
                extracted[i] = image_to_save
                return extracted
        extracted.append(image_to_save)
        extracted.flush()
        return extracted

#this is just a small function to count error from catagory process
def errorcounting(results):
    number = 0
    for result in results:
        if 'error' in result:
            number +=1
    print('total number of obj processed:',len(results))
    print('number of failed obj',number)


def find_data(name,hdu):
    for i,image in enumerate(hdu):
        if name == image.name:
            return i,image

def K_lambda(line='Ha'):
    """
    Calculate the dust attenuation value k(lambda) from the Calzetti et al. (2000) attenuation curve
    for either H-alpha (656.3 nm) or H-beta (486.1 nm).
    
    Parameters:
        line (str): 'Ha' for H-alpha (656.3 nm) or 'Hb' for H-beta (486.1 nm).
        
    Returns:
        float: Attenuation value k(λ)
    """
    if line == 'Ha':
        wavelength_um = 0.6563  # Hα (656.3 nm)
    elif line == 'Hb':
        wavelength_um = 0.4861  # Hβ (486.1 nm)
    else:
        raise ValueError("Invalid input! Use 'Ha' for H-alpha or 'Hb' for H-beta.")

    if 0.12 <= wavelength_um <= 0.63:
        # UV to optical range (Calzetti 2000)
        k_lambda = 2.659 * (-2.156 + 1.509 / wavelength_um - 0.198 / (wavelength_um ** 2) + 0.011 / (wavelength_um ** 3)) + 4.05
    elif 0.63 < wavelength_um <= 2.2:
        # Near-Infrared range (Calzetti 2000)
        k_lambda = 2.659 * (-1.857 + 1.040 / wavelength_um) + 4.05
    else:
        raise ValueError("Wavelength out of range (0.12 - 2.2 μm).")

    return k_lambda
