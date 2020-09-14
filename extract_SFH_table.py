#!/usr/bin/env python

#./extract_SFH_table.py /home/espinosa/tmp/seg_Ha_EW.NGC5947.fits.gz /home/espinosa/CALIFA_DATA/eCALIFA/sfh_files/NGC5947.SFH.cube.fits.gz /home/espinosa/tmp/

import argparse
import numpy as np
import pandas as pd
from datetime import date
from CALIFA_utils import read_seg_map, read_SFH_fits

def extract_SFH_table(seg_map, sfh_file, output, verbose=False,
                      obj_name_from_header=False):
    head, data = read_SFH_fits(sfh_file, header=True, verbose=verbose)
    hdr, seg_map = read_seg_map(seg_map, header=True, verbose=verbose)
    if head is None or data is None:
        print('Error with SFH file {}'.format(sfh_file))
        return None
    if obj_name_from_header:
        name_sfh = head['OBJECT']
    else:
        name_sfh = sfh_file.split('/')[-1][:-17]
    ns = int(np.max(seg_map))
    if ns > 0:
        name_seg = hdr['OBJECT']
        if name_sfh != name_seg:
            print('OBJECT in header files does not match')
        obj_name = name_sfh
        print('Starting extract SFH table for galaxy {}'.format(obj_name))
        nr = np.unique(seg_map)[1:]
        nz = data.shape[0]
        seg = np.copy(seg_map)
        means_cube = np.empty((nz, int(nr.max())))
        label_rows = []
        for n in nr:
            means_array_per_region = np.array([])
            #  print('current region',n)
            mask_region = seg == n
            label = '{}-{}'.format(obj_name, int(n))
            label_rows.append(label)
            for slide in np.arange(nz):
                #  print('\t current slide',slide)
                data[slide][data[slide] < -10] = 0
                data_sec = np.ma.array(data[slide], mask=~mask_region)
                #  if data_sec[data_sec > -1e308].shape[0] == 0:
                #    slide_mean = np.nan
                #  else:
                #    slide_mean = np.mean(data_sec[data_sec > -1e308])
                slide_mean = np.ma.mean(data_sec)
                means_array_per_region = np.append(means_array_per_region,
                                                   slide_mean)
            means_cube[:, int(n-1)] = means_array_per_region
        label_col = []
        header_csv = []
        for i in np.arange(nz):
            label = head['DESC_{}'.format(i)]
            header_csv.append(label)
            if label.startswith('Luminosity Fraction for age-met '):
                label = label.replace('Luminosity Fraction for ', '')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
            elif label.startswith('Luminosity Fraction for age '):
                label = label.replace('Luminosity Fraction for ', '')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
            elif label.startswith('Luminosity Fraction for met '):
                label = label.replace('Luminosity Fraction for ', '')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
            elif label.startswith('Error in the Lum. Fract. for age-met '):
                label = label.replace('Error in the Lum. Fract. for ',
                                      'e ')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
            elif label.startswith('Error in the Lum. Fract. for age '):
                label = label.replace('Error in the Lum. Fract. for ',
                                      'e ')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
            elif label.startswith('Error in the Lum. Fract. for met '):
                label = label.replace('Error in the Lum. Fract. for ',
                                      'e ')
                label = label.replace(' SSP', '').replace(' ', '_')
                label_col.append(label)
        df = pd.DataFrame(means_cube.T)
        df.index = label_rows
        df.columns = label_col
        if verbose:
            print('Output path: ' + output)
        output_name = "HII.{}.SFH.csv".format(obj_name)
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
            fp.write("# DATE: {}-{}-{}\n".format(day, month,
                                                 year))
            fp.write("# VERSION: 1.0\n")
            fp.write("# COLAPRV: S.F.Sanchez\n")
            fp.write("# PUBAPRV: None\n")           
            fp.write('#  COLUMN1:  HIIREGID           , string, ,  HII region ID\n')
            for i, label in enumerate(header_csv):
                fp.write('#  COLUMN{}:  {}, float, {} \n'.format(i+2,
                                                                 label_col[i],
                                                                 label))
            df.to_csv(fp, index_label='HIIREGID', header=False)
    else:
        obj_name = name_sfh
        print("No regions detected for {}".format(obj_name))
    print('Extract SFH table finish for {}'.format(obj_name))

if __name__ == "__main__":
    description = "Extract the mean SFH values for every ionized regions from segmentation map and fe file of CALIFA survey"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('SEG_MAP', help='Segmentation path file')
    parser.add_argument('SFH_FILE', help='SFH path file')
    parser.add_argument('OUTPUT', help='outpat path directory')
    parser.add_argument('--verbose', help="Level of verbose: 'True'"+
                        " or 'False'", default=False, metavar='verbose')
    args = parser.parse_args()
    seg_map = args.SEG_MAP
    sfh_file = args.SFH_FILE
    output = args.OUTPUT
    verbose = args.verbose
    extract_SFH_table(seg_map, sfh_file, output, verbose)
