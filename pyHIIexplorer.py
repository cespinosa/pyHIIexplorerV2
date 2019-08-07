#!/usr/bin/env python
import argparse
import os
import warnings
import numpy as np
import logging
from astropy.io import fits

class pyHIIexplorer:
    def __init__(self, Ha_map, nx, ny, max_dist, F_max, XC=None, YC=None,
                 log_level='info'):
        self.Ha_map = Ha_map
        self.max_area = np.pi * max_dist**2
        self.F_max = F_max
        self.nx = nx
        self.ny = ny
        self.XC = XC
        self.YC = YC
        self.logger = logging.getLogger('pyHIIexplorer')
        ch = logging.StreamHandler()
        if log_level == 'info':
            self.logger.setLevel(level=logging.INFO)
            ch.setLevel(logging.INFO)
        if log_level == 'debug':
            self.logger.setLevel(level=logging.DEBUG)
            ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def HIIrecover(self):
        Ha_map = self.Ha_map
        nx = self.nx
        ny = self.ny
        XC = self.XC
        YC = self.YC
        F_max = self.F_max
        Ha_map[np.isnan(Ha_map)] = 0
        seg_map = np.zeros_like(Ha_map)
        mask_map = np.ones_like(Ha_map)
        n_reg = 0
        map_data_now = np.copy(Ha_map)
        i_PEAKS = np.array([])
        j_PEAKS = np.array([])
        flux_PEAKS = np.array([])
        stop_flag = True
        n_good = self.check_good_px(Ha_map, F_l=0, F_m=1e30, nx=nx, ny=ny)
        self.logger.info('Good Pixels (F > {} and F < {} = {}'.format(0,
                                                                  1e30,
                                                                  n_good))
        self.logger.info('[{} x {}](ny x nx))'.format(ny, nx))
        while stop_flag:
            ip = -1
            jp = -1
            low_limit = F_max
            up_limit = 1e30
            flux_iter = np.nditer(map_data_now, flags=['c_index',
                                                       'multi_index'])
            for flux in flux_iter:
                j = flux_iter.multi_index[0]
                i = flux_iter.multi_index[1]
                if flux > low_limit:
        

    def check_good_px(self, map_data, F_l, F_m, nx, ny):
        n_good = 0
        for data in np.nditer(map_data):
            if data > F_l and data < 1e30:
                n_good += 1
        return n_good
            
        
        
    

if __name__ == "__main__":
    description = 'Identifies clumpy structures on a line emission map'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('HA_MAP', help='Ha emission map PATH')
    parser.add_argument('MAX_DIST', type=float, help='max distance tbw')
    args = parser.parse_args()
    Ha_map_path = args.HA_MAP
    max_dist = args.MAX_DIST
    pyHIIexplorer(Ha_map_path, max_dist).HIIrecover()
    print('Fin')
