from pyHIIexplorer import pyHIIexplorer
from extract_flux_elines_tables import extract_flux_elines_table
from extract_SFH_table import extract_SFH_table
from extract_SSP_table import extract_SSP_table

log_level='info'

fe_dir = '/home/espinosa/CALIFA_DATA/dataproducts/fe_files/'
sfh_dir = '/home/espinosa/CALIFA_DATA/dataproducts/sfh_files/'
ssp_dir = '/home/espinosa/CALIFA_DATA/dataproducts/ssp_files/'
fe_file = fe_dir + 'flux_elines.NGC5947.cube.fits.gz'
sfh_file = sfh_dir + 'NGC5947.SFH.cube.fits.gz'
ssp_file = ssp_dir + 'NGC5947.SSP.cube.fits.gz'


max_dist = 5.5
frac_peak = 0.05
F_max = 0.3
dist = 0
min_flux = 0.05
output_path = '/home/espinosa/tmp/'

a = pyHIIexplorer(fe_file, output_path, max_dist, frac_peak, F_max, dist, min_flux,
                  log_level=log_level, n_index=45)
a.HIIrecover()

# seg_map = output_path + "seg_Ha_EW.{}.fits.gz".format(obj_name)

# extract_flux_elines_table(obj_name, seg_map, fe_file, output_path, log_level)
# extract_SFH_table(obj_name, seg_map, sfh_file, output_path, log_level)
# extract_SSP_table(obj_name, seg_map, ssp_file, fe_file, output_path, log_level)
