#!/usr/bin/env python
import logging
import numpy as np
import scipy.constants as scts
from datetime import date
from CALIFA_utils import read_flux_elines_cubes, read_seg_map

def extract_flux_elines_table(obj_name, seg_map, fe_file, output):
    logger = logging.getLogger('extract_flux_elines_table')
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
    spol = scts.speed_of_light/1000
    head, data = read_flux_elines_cubes(obj_name, fe_file,
                                        header=True, log_level=log_level)
    (nz_dt, ny_dt, nx_dt) = data.shape
    crval1 = header['CRVAL1']
    cdelta1 = header['CDELT1']
    crpix1 = header['CRPIX1']
    crval2 = header['CRVAL2']
    cdelta2 = header['CDELT2']
    crpix2 = header['CRPIX2']
    name = header['OBJECT']
    califaID = header['CALIFAID']
    label_list = ('HIERARCH V500 PPAK P1 GRAT_ID',
                  'HIERARCH V500 PPAK P3F1 GRAT_ID',
                  'HIERARCH V500 PPAK P3F2 GRAT_ID',
                  'HIERARCH V500 PPAK P3F3 GRAT_ID',
                  'HIERARCH V500 PPAK P6F2 GRAT_ID',
                  'HIERARCH V500 PPAK P7F11 GRAT_ID',
                  'HIERARCH V500 PPAK P9F2 GRAT_ID')
    for label in label_list:
        try:
            grat = header[label]
            if verbose:
                logger.debug("Using GRAT_ID label:{}".format(label))
            break
        except Exception:
            grat = None
    if grat is None:
        logger.warn("GRAT ID label is not in the label list")
        logger.warn("Table for {} is not created".format(obj_name))
        return None
    fwhm_inst = 2.3
    if grat == 9:
        fwhm_inst = 6.0
    wavelengths = np.loadtxt('emission_lines.LIST',
                             comments='#', usecols=[0])
    name_elines = np.loadtxt('emission_lines.LIST', comments='#',
                             usecols=[1], dtype=np.str)
    nz_Ha = np.where(name_elines == 'Ha')[0].item()
    NZ = wavelengths.size
    logger.debug('Index of Halpha = {}'.format(nz_Ha))
    NZ_CUBE = nz_dt/8
    if NZ_CUBE != NZ:
        looger.warn("Wrong emission line table! ({} != {})".format(NZ, NZ_CUBE))
        return None
    data_seg_map = read_seg_map(seg_map, log_level='info')
    (ny_seg, nx_seg) = data_seg_map.shape
    if nx_dt > nx_seg or ny_dt > ny_seg:
        logger.warn("Dimension doesn't match ({}x{})!=({}x{})".format(nx_dt,
                                                                ny_dt,
                                                                nx_seg,
                                                                ny_seg))
        return None
    ns = int(np.max(data_seg_map))
    regions_id = np.unique(data_seg_map)
    if ns > 0:
        logger.debug("{} sources detected in {}".format(ns, obj_name))
        a_out = np.zeros((nz_dt, ns))
        a_out_med = np.zeros((nz_dt, ns))
        a_out_sq = np.zeros((nz_dt, ns))
        seg = np.copy(data_seg_map)
        x = np.array([])
        y = np.array([])
        weights = data[nz_Ha, :, :]
        weights = weights**2
        for i, region in enumerate(regions_id[1:]):
            npt = (seg == region).sum()
            data_region_mask = seg == region
            ave = np.average(np.where(data_region_mask)[1],
                             weights=weights[data_region_mask])
            x = np.append(x, ave)
            ave = np.average(np.where(data_region_mask)[0],
                             weights=weights[data_region_mask])
            y = np.append(y, ave)
            mean_array = np.array([])
            mean_array_sq = np.array([])
            for slide in np.arange(nz_dt):
                data_sec = np.ma.array(data[slide],
                                       mask=~data_region_mask)
                ave = np.ma.average(data_sec,
                                    weights=weights*data_region_mask)
                mean_array = np.append(mean_array, ave)
                ave = np.ma.average(data_sec**2,
                                    weights=weights*data_region_mask)
                mean_array_sq = np.append(mean_array_sq, ave)
            a_out[:, i] = mean_array * npt
            a_out_med[:, i] = mean_array
            a_out_sq[:, i] = mean_array_sq/npt*(1+1.6*np.log(npt))
        logger.debug('Output path: ' + output)
        output_name = "HII." + obj_name + ".flux_elines.csv"
        with open(output + output_name, "wt") as fp:
            if verbose:
                print("Writing csv file in:{}".format(output
                                                      + output_name))
            today = date.today()
            day = today.day
            month = today.month
            year = today.year
            fp.write("# AUTHOR: C. Espinosa\n")
            fp.write("# SOURCE: QCVC\n")
            fp.write("# DATE: {}-{}-{}".format(day, month, year))
            fp.write("# VERSION: 1.0\n")
            fp.write("# COLAPRV: S.F.Sanchez\n")
            fp.write("# PUBAPRV: None\n")
            fp.write("#  COLUMN1:  HIIREGID           ,"
                     + " string,  ,  HII region ID\n")
            fp.write("#  COLUMN2:  GALNAME            ,"
                     + " string,  ,  The name of the galaxy\n")
            fp.write("#  COLUMN3:  CALIFAID           ,"
                     + " int   ,  ,  CALIFAID\n")
            fp.write("#  COLUMN4:  X                  ,"
                     + " float, arcsec  ,"
                     + "X Centroid of the HII region in the cube\n")
            fp.write("#  COLUMN5:  Y                  ,"
                     + " float, arcsec  ,"
                     + " Y Centroid of the HII region in the cube\n")
            fp.write("#  COLUMN6:  RA                 ,"
                     + " float, degree  ,"
                     + " RA of the Centroid of the HII region\n")
            fp.write("#  COLUMN7:  DEC                ,"
                     + " float, degree  ,"
                     + " DEC of the Centroid of the HII region\n")
            fp.write("#  COLUMN8:  Ha_VEL             ,"
                     + "  float, km/s  , Ha velocity in km/s\n")
            fp.write("#  COLUMN9:  Ha_DISP            ,"
                     + "  float, km/s  ,"
                     + " Ha velocity dispersion-sigma in km/s\n")
            
if __name__ == "__main__":
