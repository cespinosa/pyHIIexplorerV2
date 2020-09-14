import os
import csv
from astropy.io import fits

def get_slice_from_flux_elines(fe_file, n_param=45, header=True,
                               verbose=False):
    if header:
        if verbose:
            print('Getting data and header from fits file')
        head, data = read_flux_elines_cubes(fe_file,
                                            header=header, verbose=verbose)
    else:
        if verbose:
            print('Getting only data from fits file')
        data = read_flux_elines_cubes(fe_file,
                                      header=header, verbose=verbose)
    
    if data is not None:
        if verbose:
            print('Getting cube slice {}'.format(n_param))
        map_slice = data[n_param, :, :]
        if header:
            return head, map_slice
        else:
            return map_slice
    else:
        if header:
            return None, None
        else:
            return None

def read_fits_file(fits_file, header=True, verbose=False):
    print("Reading {}".format(fits_file))
    if os.path.isfile(fits_file):
        if header:
            data, header = fits.getdata(fits_file, header=header)
            if verbose:
                print('Read fits file Done')
            return header, data
        else:
            data = fits.getdata(fits_file, header=header)
            if verbose:
                print('Read fits file Done')
            return data
    else:
        print('fits file not found')
        if header:
            return None, None
        else:
            return None

def read_flux_elines_cubes(fe_file, header=True, verbose=False):
    print("Reading flux eline cube: {}".format(fe_file))
    if os.path.isfile(fe_file):
        if header:
            data, header = fits.getdata(fe_file, header=header)
            if verbose:
                print('Read flux elines done')
            return header, data
        else:
            data = fits.getdata(fe_file, header=header)
            if verbose:
                print('Read flux elines done')
            return data
    else:
        print('File not found: {}'.format(fe_file))
        if header:
            return None, None
        else:
            return None

def get_center(obj_name, verbose=False):
    XC = 0
    YC = 0
    dir_path = '/home/espinosa/CALIFA_DATA/eCALIFA/'
    file_name = 'get_proc_elines_CALIFA.clean.csv'
    if verbose:
        print("File path {}".format(dir_path + file_name))
    try:
        with open(dir_path + file_name, 'r') as fileCAL:
            reader = csv.reader(fileCAL)
            for row in reader:
                if row[0] == obj_name:
                    XC = float(row[167])
                    YC = float(row[168])
                    break
        if XC == 0:
            print('{} is not in'.format(obj_name) +
                  'get_proc_elines_CALIFA.csv file')
            print('Setting XC=YC=0')
        else:
            print('{} found in get_proc_elines'.format(obj_name))
    except FileNotFoundError:
        print('get_proc_elines_CALIFA.csv not found')
    except:
        print('Something wrong  with get center function')
    return XC, YC

def read_seg_map(seg_map, header=False, verbose=True):
    file_path = seg_map
    if os.path.isfile(file_path):
        if header:
            data, header = fits.getdata(file_path, header=header)
        else:
            data = fits.getdata(file_path, header=header)
    else:
        print('No exist segmentation map file: {}'.format(seg_map))
        if header:
            return None, [0]
        else:
            [0]
    print('Read Segmentation map file done')
    if header:
        return header, data
    else:
        return data

def read_mask_map(mask_map, header=False, verbose=False):
    file_path = mask_map
    if os.path.isfile(file_path):
        if header:
            data, header = fits.getdata(file_path, header=header)
        else:
            data = fits.getdata(file_path, header=header)
    else:
        print('No exist masked map file')
        if header:
            return None, [0]
        else:
            [0]
    print('Read masked map file done')
    if header:
        return header, data
    else:
        return data

    
def read_SSP_fits(ssp_file, header=True, verbose=False):
    print("Reading SSP cube: {}".format(ssp_file))
    file_path = ssp_file
    if os.path.isfile(file_path):
        try:
            if header:
                data, head = fits.getdata(file_path, header=header)
                if verbose:
                    print('Read SSP done')
                return head, data
            else:
                data = fits.getdata(file_path, header=header)
                if verbose:
                    print('Read SSP done')
                return data
        except Exception:
            print('Something wrong with SSP FITS file for' +
                  '{}'.format(obj_name))
            if header:
                return None, None
            else:
                return None
    else:
        print('SSP Fits file not found for {}'.format(obj_name))
        if header:
            return None, None
        else:
            return None

def read_SFH_fits(sfh_file, header=True, verbose=False):
    print("Reading SFH cube: {}".format(sfh_file))
    file_path = sfh_file
    if os.path.isfile(file_path):
        try:
            if header:
                data, head = fits.getdata(file_path, header=header)
            else:
                data = fits.getdata(file_path, header=header)
            if header:
                return head, data
            else:
                return data
        except Exception:
            print('Something wrong with SFH FITS file')
            if header:
                return None, None
            else:
                return None
    else:
        print('SFH Fits file not found')
        if header:
            return None, None
        else:
            return None
