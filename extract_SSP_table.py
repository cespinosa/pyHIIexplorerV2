#!/usr/bin/env python

# ./extract_SSP_table.py /home/espinosa/tmp/seg_Ha_EW.NGC5947.fits.gz /home/espinosa/CALIFA_DATA/eCALIFA/ssp_files/NGC5947.SSP.cube.fits.gz /home/espinosa/CALIFA_DATA/eCALIFA/fe_files/flux_elines.NGC5947.cube.fits.gz /home/espinosa/tmp/
 

import numpy as np
import logging
import argparse
import csv
from datetime import date
from CALIFA_utils import read_seg_map, read_SSP_fits, get_slice_from_flux_elines

def extract_SSP_table(seg_map, ssp_file, fe_file,  output, log_level,
                      obj_name_from_header=False):
    logger = logging.getLogger('extract_SSP_table')
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
    
    head, data = read_SSP_fits(ssp_file, header=True, log_level=log_level)
    hd_fe, dt_fe = get_slice_from_flux_elines(fe_file, 45,
                                              header=True, log_level=log_level)
    hdr, seg_map = read_seg_map(seg_map, header=True, log_level='info')
    ns = int(np.max(seg_map))
    name_elines = np.loadtxt('emission_lines.LIST', comments='#',
                             usecols=[1], dtype=np.str)
    nz_Ha = np.where(name_elines == 'Ha')[0].item()
    if obj_name_from_header:
        name_ssp = head['OBJECT']
    else:
        name_ssp = fe_file.split('/')[-1][12:-13]
    if ns > 0:
        name_seg = hdr['OBJECT']
        if name_ssp != name_seg:
            logger.warn('OBJECT in header files does not match')
        obj_name = name_ssp
        logger.info('Starting extract SSP table for galaxy {}'.format(obj_name))
        
        (nz_dt, ny_dt, nx_dt) = data.shape
        crval1 = hd_fe['CRVAL1']
        cdelt1 = hd_fe['CDELT1']
        crpix1 = hd_fe['CRPIX1']
        crval2 = hd_fe['CRVAL2']
        cdelt2 = hd_fe['CDELT2']
        crpix2 = hd_fe['CRPIX2']
        NAME = name_ssp
        CALIFAID = hd_fe['CALIFAID']
        nr = np.unique(seg_map)[1:]
        nz = data.shape[0]
        seg = np.copy(seg_map)
        a_out = np.zeros((nz_dt, ns))
        a_out_med = np.zeros((nz_dt, ns))
        a_out_sq = np.zeros((nz_dt, ns))
        seg = np.copy(seg_map)
        x = np.array([])
        y = np.array([])
        weights = dt_fe[ :, :]
        weights = weights**2

        for i, region in enumerate(nr):
            npt = (seg == region).sum()
            data_region_mask = seg == region
            weights_masked_ = weights[data_region_mask]
            if weights_masked_.sum() == 0:
                weights_masked = np.ones_like(weights_masked_)
                weights_tmask = np.ones_like(weights)
            else:
                weights_masked = weights_masked_
                weights_tmask = weights*data_region_mask
            x_all = np.where(data_region_mask)[1]
            x = np.append(x, np.average(x_all, weights=weights_masked))
            y_all = np.where(data_region_mask)[0]
            y = np.append(y, np.average(y_all, weights=weights_masked))
            mean_array = np.array([])
            mean_array_sq = np.array([])
            for slide in np.arange(data.shape[0]):
                data_sec = np.ma.array(data[slide],
                                       mask=~data_region_mask)
                ave = np.ma.average(data_sec, weights=weights_tmask)
                mean_array = np.append(mean_array, ave)
                ave = np.ma.average(data_sec**2, weights=weights_tmask)
                mean_array_sq = np.append(mean_array_sq, ave)
            a_out[:, i] = mean_array * npt
            a_out_med[:, i] = mean_array
            a_out_sq[:, i] = mean_array_sq/npt
        logger.debug('Output path: ' + output)
        output_name = "HII.{}.SSP.csv".format(obj_name)
        with open(output + output_name, "wt") as fp:
            logger.debug("Writing csv file in:{}".format(output
                                                         + output_name))
            today = date.today()
            day = today.day
            month = today.month
            year = today.year
            fp.write("# AUTHOR: C. Espinosa\n")
            fp.write("# SOURCE: QCVC\n")
            fp.write("# DATE: {}-{}-{}\n".format(day, month,
                                                 year))
            fp.write("# VERSION: 1.0\n")
            fp.write("# COLAPRV: S.F.Sanchez\n")
            fp.write("# PUBAPRV: None\n")
            fp.write("#  COLUMN1:  HIIREGID           , string, "
                     + ",  HII region ID\n")
            fp.write("#  COLUMN2:  GALNAME            , string, "
                     + ",  The name of the galaxy\n")
            fp.write("#  COLUMN3:  CALIFAID           , int   , "
                     + ",  CALIFAID\n")
            fp.write("#  COLUMN4:  X                  , float, arcsec "
                     + ", X Centroid of the HII region in the cube\n")
            fp.write("#  COLUMN5:  Y                  , float, arcsec "
                     + ", Y Centroid of the HII region in the cube\n")
            fp.write("#  COLUMN6:  RA                 , float, degree "
                     + ", RA of the Centroid of the HII region\n")
            fp.write("#  COLUMN7:  DEC                , float, degree "
                     + ", DEC of the Centroid of the HII region\n")
            fp.write("#  COLUMN8:  log_age_LW         , float,"
                     + "log(age/yr), "
                     + "luminosity weighted age of"
                     + "the stellar population\n")
            fp.write("#  COLUMN9:  log_age_MW         , float,"
                     + "log(age/yr), "
                     + "mass weighted age of the stellar population \n")
            fp.write("#  COLUMN10: e_log_age          , float,"
                     + "dex,  error of the age of"
                     + "the stellar population \n")
            fp.write("#  COLUMN11: ZH_LW              , float,"
                     + "[Z/H], "
                     + "luminosity weighted metallicity of the"
                     + "stellar population \n")
            fp.write("#  COLUMN12: ZH_MW              , float,"
                     + "[Z/H],  mass weighted metallicity of the stellar"
                     + "population \n")
            fp.write("#  COLUMN13: e_ZH               , float,"
                     + "dex, "
                     + "error metallicity of the stellar population \n")
            fp.write("#  COLUMN14: AV_ssp             , float,"
                     + "mag,  average dust attnuation of the stellar"
                     + "population \n")
            fp.write("#  COLUMN15: e_AV_ssp           , float,"
                     + "mag,  error of the average dust attnuation of"
                     + "the stellar population \n")
            fp.write("#  COLUMN16: vel_ssp            , float,"
                     + "km/s,  velocity of the stellar population \n")
            fp.write("#  COLUMN17: e_vel_ssp          , float,"
                     + "km/s,  error in the velocity of the stellar"
                     + "population \n")
            fp.write("#  COLUMN18: disp_ssp           , float,"
                     + "km/s,  velocity dispersion of the stellar"
                     + "population \n")
            fp.write("#  COLUMN19: e_disp_ssp         , float,"
                     + "km/s,  error in velocity dispersion of the"
                     + "stellar population \n")
            fp.write("#  COLUMN20: log_ML             , float,"
                     + "log(M_sun/L_sun),  average mass-to-light"
                     + "ratio of the stellar population \n")
            fp.write("#  COLUMN21: log_Sigma_Mass     , float,"
                     + "log(M_sun/arcsec^2),  stellar mass density \n")
            fp.write("#  COLUMN22: log_Sigma_Mass_corr, float,"
                     + "log(M_sun/arcsec^2),  stellar mass density"
                     + "dust corrected\n")
            writer = csv.writer(fp)
            for i, n in enumerate(nr):
                Row_lines = []
                DEC = crval2+cdelt2*(y[i]-(crpix2-1))/3600
                cos_ang = np.cos(np.pi*DEC/180)
                RA = crval1-cdelt1*(x[i]-(crpix1-1))/3600/cos_ang
                HIIREGID = str(NAME)+'-'+str(int(n))
                Row_lines.append(HIIREGID)
                Row_lines.append(NAME)
                Row_lines.append(CALIFAID)
                Row_lines.append(x[i])
                Row_lines.append(y[i])
                Row_lines.append(RA)
                Row_lines.append(DEC)
                for j in np.arange(5, nz):
                    val_sq = a_out_sq[j, i]
                    val_med = a_out_med[j, i]
                    if head['DESC_{}'.format(j)].split()[0] == 'error':
                        val = val_sq
                    else:
                        val = val_med
                    Row_lines.append(val)
                writer.writerow(Row_lines)
    else:
        obj_name = name_ssp
        logger.info("No regions detected for {}".format(obj_name))
    logger.info('Extract SSP table finish for {}'.format(obj_name))

if __name__ == "__main__":
    description = "Extract the mean SSP values for every ionized regions from segmentation map and fe file of CALIFA survey"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('SEG_MAP', help='Segmentation path file')
    parser.add_argument('SSP_FILE', help='SSP path file')
    parser.add_argument('FE_FILE', help='flux elines path file')
    parser.add_argument('OUTPUT', help='outpat path directory')
    parser.add_argument('--log_level', help="Level of verbose: 'info'"+
                        " or 'debug'", default='info', metavar='level')
    args = parser.parse_args()
    seg_map = args.SEG_MAP
    ssp_file = args.SSP_FILE
    fe_file = args.FE_FILE
    output = args.OUTPUT
    log_level = args.log_level
    extract_SSP_table(seg_map, ssp_file, fe_file, output, log_level)
