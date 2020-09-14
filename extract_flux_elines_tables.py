#!/usr/bin/env python

#./extract_flux_elines_tables.py /home/espinosa/tmp/seg_Ha_EW.NGC5947.fits.gz /home/espinosa/CALIFA_DATA/eCALIFA/fe_files/flux_elines.NGC5947.cube.fits.gz /home/espinosa/tmp/

import argparse
import csv
import numpy as np
import scipy.constants as scts
from datetime import date
from CALIFA_utils import read_flux_elines_cubes, read_seg_map

def extract_flux_elines_table(seg_map, fe_file, output, verbose,
                              obj_name_from_header=False):
    spol = scts.speed_of_light/1000
    header, data = read_flux_elines_cubes(fe_file,
                                          header=True, verbose=verbose)
    (nz_dt, ny_dt, nx_dt) = data.shape
    if obj_name_from_header:
        name_fe = header['OBJECT']
    else:
        name_fe = fe_file.split('/')[-1][12:-13] 
    hdr, data_seg_map = read_seg_map(seg_map, header=True, verbose=verbose)
    ns = int(np.max(data_seg_map))
    if verbose:
        print('Clumpy regions found: {}'.format(ns))
    if ns > 0:
        name_seg = hdr['OBJECT']
        if name_fe != name_seg:
            print('OBJECT in header files does not match')
        obj_name = name_fe
        print('Starting extract fe table for galaxy {}'.format(obj_name))
        crval1 = header['CRVAL1']
        cdelta1 = header['CDELT1']
        crpix1 = header['CRPIX1']
        crval2 = header['CRVAL2']
        cdelta2 = header['CDELT2']
        crpix2 = header['CRPIX2']
        # name = header['OBJECT']
        name = name_fe
        califaID = header['CALIFAID']
        # label_list = ('HIERARCH V500 PPAK P1 GRAT_ID',
        #              'HIERARCH V500 PPAK P3F1 GRAT_ID',
        #              'HIERARCH V500 PPAK P3F2 GRAT_ID',
        #              'HIERARCH V500 PPAK P3F3 GRAT_ID',
        #              'HIERARCH V500 PPAK P6F2 GRAT_ID',
        #              'HIERARCH V500 PPAK P7F11 GRAT_ID',
        #              'HIERARCH V500 PPAK P9F2 GRAT_ID')
        # for label in label_list:
        #     try:
        #         grat = header[label]
        #         logger.debug("Using GRAT_ID label:{}".format(label))
        #         break
        #     except Exception:
        #         grat = None
        # if grat is None:
        #     logger.warn("GRAT ID label is not in the label list")
        #     logger.warn("Table for {} is not created".format(obj_name))
        #     return None
        # fwhm_inst = 2.3
        # if grat == 9:
        #     fwhm_inst = 6.0
        wavelengths = np.loadtxt('emission_lines.LIST',
                                 comments='#', usecols=[0])
        name_elines = np.loadtxt('emission_lines.LIST', comments='#',
                                 usecols=[1], dtype=np.str)
        nz_Ha = np.where(name_elines == 'Ha')[0].item()
        NZ = wavelengths.size
        if verbose:
            print('Index of Halpha = {}'.format(nz_Ha))
        NZ_CUBE = nz_dt/8
        if NZ_CUBE != NZ:
            print("Wrong emission line table! ({} != {})".format(NZ, NZ_CUBE))
            return None
        (ny_seg, nx_seg) = data_seg_map.shape
        if nx_dt > nx_seg or ny_dt > ny_seg:
            print("Dimension doesn't match ({}x{})!=({}x{})".format(nx_dt,
                                                                    ny_dt,
                                                                    nx_seg,
                                                                    ny_seg))
            return None

        regions_id = np.unique(data_seg_map)
        if verbose:
            print("{} sources detected in {}".format(ns, obj_name))
        a_out = np.zeros((nz_dt, ns))
        a_out_med = np.zeros((nz_dt, ns))
        a_out_sq = np.zeros((nz_dt, ns))
        seg = np.copy(data_seg_map)
        x = np.array([])
        y = np.array([])
        weights = data[nz_Ha, :, :]
        weights = weights**2
        n_pts = np.array([])
        for i, region in enumerate(regions_id[1:]):
            npt = (seg == region).sum()
            n_pts = np.append(n_pts, npt)
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
            # a_out_sq[:, i] = mean_array_sq/npt*(1+1.6*np.log(npt))
            a_out_sq[:, i] = mean_array_sq
        if verbose:
            print('Output path: ' + output)
        output_name = "HII.{}.flux_elines.csv".format(obj_name)
        # print(n_pts)
        # print(n_pts.shape)
        # return None
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
            fp.write("# DATE: {}-{}-{}\n".format(year, month, day))
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
            # fp.write("#  COLUMN8:  Ha_VEL             ,"
            #          + "  float, km/s  , Ha velocity in km/s\n")
            # fp.write("#  COLUMN9:  Ha_DISP            ,"
            #          + "  float, km/s  ,"
            #          + " Ha velocity dispersion-sigma in km/s\n")
            NC = 8
            writer = csv.writer(fp)
            for i, wl in enumerate(wavelengths):
                for I in [i, i+4*NZ, i+3*NZ, i+7*NZ]:
                          # i+NZ, i+5*NZ, i+2*NZ, i+6*NZ]:
                    header_now = header['NAME{}'.format(I)]
                    for s_char in ' []':
                        header_now = header_now.replace(s_char, '')
                    junk = header_now.split('I')
                    if len(junk) > 1:
                        temp = header_now.split('I')[-1]
                        header_now = header_now.replace(temp, '')
                    else:
                        header_now = header_now
                    if "vel" in header_now:
                        units = "km/s"
                    if "disp" in header_now:
                        # header_now = "sigma corr"+header_now
                        units = "km/s"
                    if "flux" in header_now:
                        units = "10^-16 erg/s/cm^2"
                    if "EW" in header_now:
                        units = "Angstrom"
                    header_now = header_now+str(int(wl))
                    fp.write("#  COLUMN{}: ".format(NC)
                             + "{}      ".format(header_now)
                             + ",float,  {}  ,".format(units)
                             + "{} for".format(header['NAME{}'.format(I)])
                             + "{}\n".format(wl))
                    NC += 1
            for i in np.arange(0, ns):
                # print(i)
                Row_lines = []
                K = i + 1
                DEC = crval2+cdelta2*(y[i]-(crpix2-1))/3600
                cos_ang = np.cos(np.pi*DEC/180)
                RA = crval1-cdelta1*(x[i]-(crpix1-1))/3600/cos_ang
                HIIREGID = str(name)+'-'+str(K)
                a_out_weights = np.copy(a_out[:NZ, i])
                a_out_vel = np.copy(a_out_med[NZ:2*NZ, i])
                a_out_disp = np.copy(a_out_med[2*NZ:3*NZ, i])
                indexs = [0, 26, 28, 45]
                a_out_weights = a_out_weights[indexs]
                a_out_vel = a_out_vel[indexs]
                a_out_disp = a_out_disp[indexs]
                a_lambda = wavelengths[indexs]                
                # vel_Ha = a_out_med[96, i]
                # disp_Ha = a_out_med[147, i]
                # disp_Ha = np.sqrt(abs(disp_Ha**2 - fwhm_inst**2))
                # disp_Ha = (disp_Ha/6563)*scts.speed_of_light/1000
                # a_out_disp = np.sqrt(abs(a_out_disp**2-fwhm_inst**2))
                # a_out_disp = (a_out_disp/a_lambda)*scts.speed_of_light/1000
                # a_out_disp[a_out_disp <= 10] = np.nan
                # vel = vel_Ha
                # disp = disp_Ha
                # if np.nanmin(a_out_disp) < disp:
                #     disp = np.nanmin(a_out_disp)
                # disp = disp/2.354
                Row_lines.append(HIIREGID)
                Row_lines.append(name)
                Row_lines.append(califaID)
                Row_lines.append(x[i])
                Row_lines.append(y[i])
                Row_lines.append(RA)
                Row_lines.append(DEC)
                # Row_lines.append(vel)
                # Row_lines.append(disp)
                for j in np.arange(NZ):
                    for J in [j, j+4*NZ, j+3*NZ, j+7*NZ]:
                              # i+NZ, i+5*NZ, i+2*NZ, i+6*NZ]:
                        val = a_out[J, i]
                        val_sq = a_out_sq[J, i]
                        val_med = a_out_med[J, i]
                        header_now = header['NAME{}'.format(J)]
                        for char in ' []':
                            header_now = header_now.replace(char, '')
                        # if "disp" in header_now:
                        #     val = val_med
                        #     val_now = np.sqrt(np.abs(val**2-fwhm_inst**2))
                        #     val = (val_now/wavelengths[j]) * spol
                        # if "vel" in header_now:
                        #     val = val_med
                        if "EW" in header_now:
                            val = val_med
                        if "e_" in header_now:
                            # val = np.sqrt(np.abs(val_sq))
                            val = np.sqrt(np.abs(val_sq)) * (1 + 1.6*np.log(n_pts[i]))
                        Row_lines.append(val)
                writer.writerow(Row_lines)
            if verbose:
                print(output_name + " created")
    else:
        obj_name = name_fe
        print("No regions detected for {}".format(obj_name))
    print("extract flux elines table finish for {}".format(obj_name))

if __name__ == "__main__":
    description = "Extract the mean flux emission lines values for every ionized regions from segmentation map and fe file of CALIFA survey"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('SEG_MAP', help='Segmentation path file')
    parser.add_argument('FE_FILE', help='flux elines path file')
    parser.add_argument('OUTPUT', help='outpat path directory')
    parser.add_argument('--verbose', help="Verbose: 'True'"+
                        " or 'False'", default=False, metavar='verbose')
    args = parser.parse_args()
    seg_map = args.SEG_MAP
    fe_file = args.FE_FILE
    output = args.OUTPUT
    verbose = args.verbose
    extract_flux_elines_table(seg_map, fe_file, output, verbose)
