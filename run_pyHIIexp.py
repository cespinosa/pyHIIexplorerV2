from pyHIIexplorer import pyHIIexplorer
from extract_flux_elines_tables import extract_flux_elines_table
from extract_SFH_table import extract_SFH_table
from extract_SSP_table import extract_SSP_table

log_level='info'

fe_dir = '/home/espinosa/CALIFA_DATA/dataproducts/fe_files/'
sfh_dir = '/home/espinosa/CALIFA_DATA/dataproducts/sfh_files/'
ssp_dir = '/home/espinosa/CALIFA_DATA/dataproducts/ssp_files/'

def extract_HIIregions(obj_name):
    max_dist=5.5
    frac_peak=0.05
    F_max=0.3
    dist=0
    min_flux=0.05
    output_path='/home/espinosa/tmp/'
    fe_file = fe_dir + 'flux_elines.{}.cube.fits.gz'.format(obj_name)
    sfh_file = sfh_dir + '{}.SFH.cube.fits.gz'.format(obj_name)
    ssp_file = ssp_dir + '{}.SSP.cube.fits.gz'.format(obj_name)

    a = pyHIIexplorer(fe_file, output_path, max_dist, frac_peak, F_max, dist,
                      min_flux, log_level=log_level, n_index=45)
    a.HIIrecover()

    seg_map = output_path + "seg_Ha_EW.{}.fits.gz".format(obj_name)

    extract_flux_elines_table(seg_map, fe_file, output_path, log_level)
    extract_SFH_table(seg_map, sfh_file, output_path, log_level)
    extract_SSP_table(seg_map, ssp_file, fe_file, output_path, log_level)
