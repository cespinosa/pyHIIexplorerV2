import glob
import numpy as np
import pandas as pd
from multiprocessing import Pool
from run_pyHIIexp import extract_HIIregions
from run_pyHIIexp import extract_HIIregions_p

eCALIFA_path = '/home/espinosa/CALIFA_DATA/eCALIFA/'
pCALIFA_path = '/home/espinosa/CALIFA_DATA/pCALIFA/'

def make_catalog_eCALIFA(obj_names, N_CORE=4):
    pool = Pool(N_CORE)
    pool.map(extract_HIIregions, obj_names, chunksize=1)
    pool.close()
    pool.join()

def make_catalog_pCALIFA(obj_names, N_CORE=4):
    pool = Pool(N_CORE)
    pool.map(extract_HIIregions_p, obj_names, chunksize=1)
    pool.close()
    pool.join()

    
def get_names_eCALIFA():
    path_files = glob.glob(eCALIFA_path + 'fe_files/'  + 'flux_elines.*')
    obj_names = [path_file.split('.', 1)[1][:-13]
             for path_file in path_files ]
    return obj_names

def get_names_get_proc_files_eCALIFA(path=None, clean_catalog=True):
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

def get_names_pCALIFA():
    path_files = glob.glob(pCALIFA_path + 'fe_files_corr/'  + 'flux_elines.*')
    obj_names = [path_file.split('.', 1)[1][:-13]
             for path_file in path_files ]
    return obj_names

def get_names_get_proc_files_pCALIFA(path=None, clean_catalog=True):
    # For PILOT sample, only one get_proc_elines table exists
    if path is None:
        path = pCALIFA_path
    if clean_catalog is True:
        name_file = path + 'get_proc_elines_PILOT.csv'
        df_clean = pd.read_csv(name_file, comment='#', header=None,
                               usecols=[0])
        galaxy_names = df_clean.values.ravel()
    elif clean_catalog is False:
        name_file = path + 'get_proc_elines_PILOT.csv'
        df = pd.read_csv(name_file, comment='#', header=None,
                         usecols=[0])
        galaxy_names = df.values.ravel()
    elif clean_catalog is None:
        print('No option for clean_catalog. Clean catalog is set by default')
        name_file = path + 'get_proc_elines_PILOT.clean.csv'
        df_clean = pd.read_csv(name_file, comment='#', header=None,
                               usecols=[0])
        galaxy_names = df_clean.values.ravel()
    return galaxy_names

def clean_sample_eCALIFA():
    fe_names = get_names_eCALIFA()
    galaxy_clean_names = get_names_get_proc_files_eCALIFA()
    obj_names = np.intersect1d(galaxy_clean_names, fe_names)
    obj_names_diff = np.setdiff1d(galaxy_clean_names, fe_names)
    # print("# missing fe files: {}".format(len(obj_names_diff)))
    return obj_names

def clean_sample_pCALIFA():
    fe_names = get_names_pCALIFA()
    galaxy_clean_names = get_names_get_proc_files_pCALIFA()
    obj_names = np.intersect1d(galaxy_clean_names, fe_names)
    obj_names_diff = np.setdiff1d(galaxy_clean_names, fe_names)
    # print("# missing fe files: {}".format(len(obj_names_diff)))
    return obj_names

