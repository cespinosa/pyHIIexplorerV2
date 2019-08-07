import os
import logging
import csv
from astropy.io import fits

def get_slice_from_flux_elines(obj_name, fe_file, n_param=45, header=True,
                               log_level='info'):
    logger = logging.getLogger('get slice from fe file')
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
    if header:
        logger.debug('Getting data and header from fits file')
        head, data = read_flux_elines_cubes(obj_name, fe_file,
                                            header=header, log_level=log_level)
    else:
        logger.debug('Getting only data from fits file')
        data = read_flux_elines_cubes(obj_name, fe_file,
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
        

def read_flux_elines_cubes(obj_name, fe_file, header=True,
                           log_level='info'):
    logger = logging.getLogger('read fe file')
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
    logger.info("{} Flux elines file ".format(obj_name)
                  + "path: {}".format(fe_file))
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
        logger.error('fe file not found for {}'.format(obj_name))
        if header:
            return None, None
        else:
            return None

def get_center(obj_name, log_level='info'):
    XC = 0
    YC = 0
    dir_path = '/home/espinosa/data/'
    file_name = 'get_proc_elines_CALIFA.clean.csv'    
    logger = logging.getLogger('get center')
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
    logger.debug("File path".format(dir_path + file_name))
    try:
        with open(dir_path + file_name, 'r') as fileCAL:
            reader = csv.reader(fileCAL)
            for row in reader:
                if row[0] == obj_name:
                    XC = float(row[167])
                    YC = float(row[168])
                    break
        if XC == 0:
            logger.warning('{} is not in'.format(obj_name) + 
                           'get_proc_elines_CALIFA.csv file')
            logger.warning('Setting XC=YC=0')
        else:
            logger.debug('{} found in get_proc_elines')
    except FileNotFoundError:
        logger.error('get_proc_elines_CALIFA.csv not found')
    except:
        logger.warning('Something wrong  with get center function')
    return XC, YC
