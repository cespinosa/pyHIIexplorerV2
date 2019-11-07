from pyHIIexplorer import pyHIIexplorer
from extract_flux_elines_tables import extract_flux_elines_table
from extract_SFH_table import extract_SFH_table
from extract_SSP_table import extract_SSP_table

log_level='info'
#eCALIFA
fe_dir = '/home/espinosa/CALIFA_DATA/dataproducts/fe_files/'
sfh_dir = '/home/espinosa/CALIFA_DATA/dataproducts/sfh_files/'
ssp_dir = '/home/espinosa/CALIFA_DATA/dataproducts/ssp_files/'
#pCALIFA
p_fe_dir = '/home/espinosa/CALIFA_DATA/pCALIFA/fe_files_corr/'
p_sfh_dir = '/home/espinosa/CALIFA_DATA/pCALIFA/sfh_files/'
p_ssp_dir = '/home/espinosa/CALIFA_DATA/pCALIFA/ssp_files/'


def extract_HIIregions(obj_name):
    max_dist=5.5
    frac_peak=0.05
    F_max=0.3
    dist=0
    min_flux=0.05
    output_path='/home/espinosa/CALIFA_DATA/eCALIFA/catalog/'
    fe_file = fe_dir + 'flux_elines.{}.cube.fits.gz'.format(obj_name)
    sfh_file = sfh_dir + '{}.SFH.cube.fits.gz'.format(obj_name)
    ssp_file = ssp_dir + '{}.SSP.cube.fits.gz'.format(obj_name)

    
    a = pyHIIexplorer(fe_file, output_path, max_dist, frac_peak, F_max, dist,
                      min_flux, log_level=log_level, n_index=45)
    a.HIIrecover()
    
    seg_map = output_path + "seg_Ha.{}.fits.gz".format(obj_name)
    print(seg_map)
    extract_flux_elines_table(seg_map, fe_file, output_path, log_level)
    extract_SFH_table(seg_map, sfh_file, output_path, log_level)
    extract_SSP_table(seg_map, ssp_file, fe_file, output_path, log_level)

def extract_HIIregions_p(obj_name):
    if obj_name in ['ngc1058', ' ngc1637', 'ngc3184', 'ngc4625', 'ngc628',
                       'ngc7771']:
            max_dist=5.5
            frac_peak=0.05
            F_max=0.3*15
            dist=0
            min_flux=0.05*15
    else:
        max_dist=5.5
        frac_peak=0.05
        F_max=0.3
        dist=0
        min_flux=0.05
    p_flag = True
    output_path='/home/espinosa/CALIFA_DATA/pCALIFA/catalog/'
    fe_file = p_fe_dir + 'flux_elines.{}.cube.fits.gz'.format(obj_name)
    sfh_file = p_sfh_dir + '{}.SFH.cube.fits.gz'.format(obj_name)
    ssp_file = p_ssp_dir + '{}.SSP.cube.fits.gz'.format(obj_name)

    a = pyHIIexplorer(fe_file, output_path, max_dist, frac_peak, F_max, dist,
                      min_flux, p_flag, log_level=log_level, n_index=45)
    a.HIIrecover()

    seg_map = output_path + "seg_Ha.{}.fits.gz".format(obj_name)

    extract_flux_elines_table(seg_map, fe_file, output_path, log_level)
    extract_SFH_table(seg_map, sfh_file, output_path, log_level)
    extract_SSP_table(seg_map, ssp_file, fe_file, output_path, log_level)
