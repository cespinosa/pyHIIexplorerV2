import glob
import numpy as np
import pandas as pd
from multiprocessing import Pool
from pyHIIexplorer import pyHIIexplorer
from extract_SFH_table import extract_SFH_table
from extract_SSP_table import extract_SSP_table
from CALIFA_utils import get_slice_from_flux_elines, get_center
from pyHIIexplorer import pyHIIexplorer
from extract_flux_elines_tables import extract_flux_elines_table

eCALIFA_path = '/home/espinosa/CALIFA_DATA/eCALIFA/'

class pyHIIexplorer_mp():
    def __init__(self, fe_file,  max_dist, frac_peak, F_max, dist,
                       min_flux, output_path, index=None):
        self.fe_file = fe_file
        self.index = index
        self.obj_name = self.get_obj_name()
        self.Ha_map = self.get_Ha_map()
        _nx, _ny = self.get_size()
        self.nx = _nx
        self.ny = _ny
        self.F_max = F_max
        self.frac_peak = frac_peak
        self.dist = dist
        self.min_flux = min_flux
        self.max_dist = max_dist
        self.output_path = output_path
        _xc, _yc = self.get_centers()
        self.XC = _xc
        self.YC = _yc

    def get_obj_name(self):
        obj_name = self.fe_file.split('.', 1)[1][:-13]
        return obj_name 
        
    def get_Ha_map(self):
        if index is None:
            print('Reading Ha map file')
            Ha_map = read_flux_elines_cubes(self.obj_name, self.fe_file,
                                             header=False)
        else:
            print('Reading index {} from fits file'.format(self.index))
            Ha_map = get_slice_from_flux_elines(self.obj_name,
                                                self.fe_fil,
                                                n_param=self.index,
                                                header=False)
        return Ha_map

    def get_size(self):
        if index is None:
            hdr, Ha_map = read_flux_elines_cubes(self.obj_name, self.fe_file)
        else:
            hdr, Ha_map = get_slice_from_flux_elines(self.obj_name,
                                                     self.fe_file,
                                                     n_param=index)
                                                     
        nx = hdr['NAXIS1']
        ny = hdr['NAXIS2']
        return nx, ny
    def get_centers(self):
        xc, yc = get_center(self.obj_name)
        return xc, yc

    

def get_names_eCALIFA():
    path_files = glob.glob(eCALIFA_path + 'fe_files/'  + 'flux_elines.*')
    obj_names = [path_file.split('.', 1)[1][:-13]
             for path_file in path_files ]
    return obj_names

def get_names_get_proc_files(path=None, clean_catalog=True):
    if path is None:
        path = eCALIFA_path
    if clean_catalog is True:
        name_file = path + 'get_proc_elines_CALIFA.clean.csv'
        df_clean = pd.read_csv(name_file, comment='#', header=None,
                               usecols=[0])
        galaxy_names = df_clean.values.ravel()
    elif clean_catalog is False:
        name_file = path + 'get_proc_elines_CALIFA.csv'
        df = pd.read_csv(name_file, comment='#', header=None,
                         usecols=[0])
        galaxy_names = df.values.ravel()
    elif clean_catalog is None:
        print('No option for clean_catalog. Clean catalog is set by default')
        name_file = path + 'get_proc_elines_CALIFA.clean.csv'
        df_clean = pd.read_csv(name_file, comment='#', header=None,
                               usecols=[0])
        galaxy_names = df_clean.values.ravel()
    return galaxy_names

def clean_sample():
    fe_names = get_names_eCALIFA()
    galaxy_clean_names = get_names_get_proc_files()
    obj_names = np.intersect1d(galaxy_clean_names, fe_names)
    obj_names_diff = np.setdiff1d(galaxy_clean_names, fe_names)
    # print("# missing fe files: {}".format(len(obj_names_diff)))
    return obj_names

def get_files_paths():
    obj_names = clean_sample()
    fe_files = [eCALIFA_path + 'fe_files/flux_elines.' + obj_name +
                '.cube.fits.gz' for obj_name in obj_names]
    return fe_files

def get_centers():
    XC = np.array([])
    YC = np.array([])
    obj_names = clean_sample()
    for name in obj_names:
        xc, yc = get_center(name)
        XC = np.append(XC, xc)
        YC = np.append(YC, yc)
    return XC, YC

def get_size():
    NX = np.array([])
    NY = np.array([])
    obj_names = clean_sample()
    fe_files = get_files_paths()
    for obj_name, fe_file in zip(obj_names, fe_files):
        hdr, Ha_map = get_slice_from_flux_elines(obj_name, fe_file)
        nx = hdr['NAXIS1']
        ny = hdr['NAXIS2']
        NX = np.append(NX, nx)
        NY = np.append(NY, ny)
    return NX, NY

def get_Ha_maps():
    obj_names = clean_sample()
    fe_files = get_files_paths()
    Ha_maps = []
    for obj_name, fe_file in zip(obj_names, fe_files):
        hdr, Ha_map = get_slice_from_flux_elines(obj_name, fe_file)
        Ha_maps.append(Ha_map[45])
    return Ha_maps
