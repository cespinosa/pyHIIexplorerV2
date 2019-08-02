#!/usr/bin/env python
import argparse
import os
import warnings
import numpy as np
from astropy.io import fits

class pyHIIexplorer:
    def __init__(self, Ha_map_path, max_dist):
        self.Ha_map_path = Ha_map_path
        self.max_area = np.pi * max_dist**2
        
    def HIIrecover(self):
        flag_file = os.path.isfile(self.Ha_map_path)
        if not flag_file:
            warn_messg = 'File not exists: {}'.format(Ha_map_path)
            warnings.warn(warn_messg)
            return None
        Ha_map, hdr = fits.getdata(Ha_map_path, header=True)
        Ha_map[np.isnan(Ha_map)] = 0
        nx = hdr['NAXIS1']
        ny = hdr['NAXIS2']
        if XC is None or YC is None:
            (XC, YC) = self.get_center()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identifies clumpy structures on a line emission map')
    parser.add_argument('HA_MAP', help='Ha emission map PATH')
    parser.add_argument('MAX_DIST', type=float, help='max distance tbw')
    args = parser.parse_args()
    Ha_map_path = args.HA_MAP
    max_dist = args.MAX_DIST
    pyHIIexplorer(Ha_map_path, max_dist).HIIrecover()
    print('Fin')
