import os
import logging
import csv
from astropy.io import fits

def get_slice_from_flux_elines(fe_file, n_param=45, header=True,
                               log_level='info'):
    logger = logging.getLogger('get slice from fe file')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(ch)
    if header:
        logger.debug('Getting data and header from fits file')
        head, data = read_flux_elines_cubes(fe_file,
                                            header=header, log_level=log_level)
    else:
        logger.debug('Getting only data from fits file')
        data = read_flux_elines_cubes(fe_file,
                                      header=header, log_level=log_level)
    
    if data is not None:
        logger.debug('Getting cube slice {}'.format(n_param))
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

def read_fits_file(fits_file, header=True, log_leve='info'):
    logger = logging.getLogger('read fits file')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(ch)
    logger.info("Reading {}".format(fits_file))
    if os.path.isfile(fits_file):
        if header:
            data, header = fits.getdata(fits_file, header=header)
            logger.debug('Read fits file Done')
            return header, data
        else:
            data = fits.getdata(fits_file, header=header)
            logger.debug('Read fits file Done')
            return data
    else:
        logger.error('fits file not found')
        if header:
            return None, None
        else:
            return None


def read_flux_elines_cubes(fe_file, header=True, log_level='info'):
    logger = logging.getLogger('read fe file')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(ch)
    logger.info("Reading {}".format(fe_file))
    if os.path.isfile(fe_file):
        if header:
            data, header = fits.getdata(fe_file, header=header)
            logger.debug('Read flux elines done')
            return header, data
        else:
            data = fits.getdata(fe_file, header=header)
            logger.debug('Read flux elines done')
            return data
    else:
        logger.error('file not found')
        if header:
            return None, None
        else:
            return None

def get_center(obj_name, log_level='info'):
    XC = 0
    YC = 0
    dir_path = '/home/espinosa/CALIFA_DATA/eCALIFA/'
    file_name = 'get_proc_elines_CALIFA.clean.csv'    
    logger = logging.getLogger('get center')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(ch)
    logger.debug("File path {}".format(dir_path + file_name))
    try:
        with open(dir_path + file_name, 'r') as fileCAL:
            reader = csv.reader(fileCAL)
            for row in reader:
                if row[0] == obj_name:
                    XC = float(row[167])
                    YC = float(row[168])
                    break
        if XC == 0:
            logger.warn('{} is not in'.format(obj_name) + 
                           'get_proc_elines_CALIFA.csv file')
            logger.warn('Setting XC=YC=0')
        else:
            logger.debug('{} found in get_proc_elines'.format(obj_name))
    except FileNotFoundError:
        logger.error('get_proc_elines_CALIFA.csv not found')
    except:
        logger.warn('Something wrong  with get center function')
    return XC, YC

def read_seg_map(seg_map, header=False, log_level='info'):
    logger = logging.getLogger('read seg map')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(ch)
    file_path = seg_map
    if os.path.isfile(file_path):
        if header:
            data, header = fits.getdata(file_path, header=header)
        else:
            data = fits.getdata(file_path, header=header)
    else:
        logger.warn('No exist segmentation map file')
        return [0]
    logger.info('Read Segmentation map file done')
    if header:
        return header, data
    else:
        return data

def read_SSP_fits(obj_name, ssp_file, header=True, log_level='info'):
    logger = logging.getLogger('read ssp file')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("{} SSP file ".format(obj_name)
                  + "path: {}".format(ssp_file))
    file_path = ssp_file
    if os.path.isfile(file_path):
        try:
            if header:
                data, head = fits.getdata(file_path, header=header)
                logger.debug('Read SSP done')
                return head, data
            else:
                data = fits.getdata(file_path, header=header)
                logger.debug('Read SSP done')
                return data
        except Exception:
            logger.warn('Something wrong with SSP FITS file for' +
                           '{}'.format(obj_name))
            if header:
                return None, None
            else:
                return None
    else:
        logger.error('SSP Fits file not found for {}'.format(obj_name))
        if header:
            return None, None
        else:
            return None

def read_SFH_fits(obj_name, sfh_file, header=True, log_level='info'):
    logger = logging.getLogger('read sfh file')
    logger.propagate = False
    ch = logging.StreamHandler()
    if log_level == 'info':
        logger.setLevel(level=logging.INFO)
        ch.setLevel(logging.INFO)
    if log_level == 'debug':
        logger.setLevel(level=logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(levelname)s %(name)s %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("{} SFH file ".format(obj_name)
                  + "path: {}".format(sfh_file))
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
            logger.warn('Something wrong with SFH FITS file for' +
                           '{}'.format(obj_name))
            if header:
                return None, None
            else:
                return None
    else:
        logger.error('SFH Fits file not found for {}'.format(obj_name))
        if header:
            return None, None
        else:
            return None
