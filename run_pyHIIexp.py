from CALIFA_utils import get_slice_from_flux_elines, get_center
from pyHIIexplorer import pyHIIexplorer

log_level='info'

obj_name='NGC5947'
fe_dir = '/home/espinosa/CALIFA_DATA/dataproducts/fe_files/'
fe_file = fe_dir + 'flux_elines.NGC5947.cube.fits.gz'

hdr, Ha_map = get_slice_from_flux_elines(obj_name, fe_file, log_level=log_level)
if Ha_map is None:
    logger.critical('Ha map error')
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']
XC, YC = get_center(obj_name, log_level=log_level)

max_dist = 5.5
frac_peak = 0.05
F_max = 0.3
dist = 0
min_flux = 0.05
output_path = '/home/espinosa/'

a = pyHIIexplorer(Ha_map, nx, ny, max_dist, frac_peak, F_max, dist, min_flux,
                  obj_name, output_path, XC, YC, log_level=log_level)
a.HIIrecover()
