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

