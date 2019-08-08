#!/usr/bin/env python
import argparse
import os
import warnings
import numpy as np
import logging
from astropy.io import fits

class pyHIIexplorer:
    def __init__(self, Ha_map, nx, ny, max_dist, frac_peak, F_max, dist,
                 min_flux, obj_name, output_path, XC=None, YC=None,
                 log_level='info', PSF=2.3):
        self.Ha_map = Ha_map
        self.max_area = np.pi * max_dist**2
        self.F_max = F_max
        self.frac_peak = frac_peak
        self.dist = dist
        self.min_flux = min_flux
        self.max_dist = max_dist
        self.nx = nx
        self.ny = ny
        self.XC = XC
        self.YC = YC
        self.PSF = PSF
        self.obj_name = obj_name
        self.output_path = output_path
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
        frac_peak = self.frac_peak
        dist = self.dist
        max_dist = self.max_dist
        min_flux = self.min_flux
        output_path = self.output_path
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
                    is_max = self.is_max(map_data=Ha_map, i=i, j=j,
                                         nx=nx, ny=ny, ip=ip, jp=jp,
                                         frac_peak=frac_peak)
                    is_near = 0
                    dist_center = np.sqrt((i - XC)**2 + (j - YC)**2)
                    if flux < up_limit and is_max == 1 and is_near == 0 \
                       and dist_center > dist:
                        low_limit = flux
                        ip = i
                        jp = j
            i_PEAKS = np.append(i_PEAKS, ip)
            j_PEAKS = np.append(j_PEAKS, jp)
            flux_PEAKS = np.append(flux_PEAKS, Ha_map[jp, ip])
            n_reg += 1
            self.logger.debug('ip = {} jp = {} nReg = {}'.format(ip, jp, n_reg))
            if ip == -1 and jp == -1:
                stop_flag = False
            else:
                self.logger.debug('Adding spaxel to peak ({},{})'.format(ip,
                                                                         jp))
                stop_flag = self.add_points(data_map=Ha_map, ip=ip, jp=jp,
                                            i_PEAKS=i_PEAKS, j_PEAKS=j_PEAKS,
                                            flux_PEAKS=flux_PEAKS,
                                            nx=nx, ny=ny, seg_map=seg_map,
                                            mask_map=mask_map, n_reg=n_reg,
                                            max_dist=max_dist,
                                            min_flux=min_flux,
                                            frac_peak=frac_peak)
                map_data_now *= mask_map
            if not stop_flag:
                self.logger.debug('Adding spaxels done')
        mask_map = np.ones_like(seg_map)
        mask = seg_map > 0
        mask = ~mask
        mask_map = mask_map * mask

        hdu_mask_map = fits.PrimaryHDU(mask_map)
        hdu_seg_map = fits.PrimaryHDU(seg_map)
        hdul_mask_map = fits.HDUList([hdu_mask_map])
        hdul_seg_map = fits.HDUList([hdu_seg_map])
        mask_name = "mask_Ha_EW.{}.fits.gz".format(self.obj_name)
        seg_name = "seg_Ha_EW.{}.fits.gz".format(self.obj_name)
        mask_path = output_path + mask_name
        seg_path = output_path + seg_name
        hdul_mask_map.writeto(mask_path, overwrite=True)
        hdul_seg_map.writeto(seg_path, overwrite=True)
        self.logger.debug('Saving fit files in {}'.format(output_path))
        self.logger.info('HIIexplorer finish for {}'.format(self.obj_name))
        logging.shu

    def check_good_px(self, map_data, F_l, F_m, nx, ny):
        n_good = 0
        for data in np.nditer(map_data):
            if data > F_l and data < 1e30:
                n_good += 1
        return n_good

    def is_max(self, map_data, i, j, nx, ny, ip=-1, jp=-1, frac_peak=0.05):
        """
        Calulcate if a pixes is maximum
        Possible values:
        1: larger than nearby
        0: equal
        -1 smaller than nearby
        -3 Default value
        """
        is_max = -3
        square_data = self.mask_square(map_data=map_data, i=i, j=j,
                                       nx=nx, ny=ny, calc_k=False)
        med_a = np.median(square_data)
        std_a = np.std(square_data, ddof=1)
        if i > 0 and i < nx and j > 0 and j < ny:
            val = map_data[j, i]
            if ip == -1 and jp == -1:
                F_peak = 0.0
            else:
                F_peak = map_data[jp, ip]
            if val > frac_peak * F_peak:
                ncut = 2.0
                if val - med_a < 0:
                    is_max = -1
                if np.abs(med_a - val) <= ncut*std_a:
                    is_max = 0
                if val > (med_a + 0.5*std_a):
                    is_max = 1
        return is_max

    def mask_square(self, map_data, i, j, nx=None, ny=None, calc_k=True):
        i1 = i - 1
        i2 = i + 2
        j1 = j - 1
        j2 = j + 2
        if i1 < 0:
            i1 = 0
            i2 = i1 + 3
        if i2 > nx -1:
            i2 = nx - 1
            i1 = i2 - 3
        if j1 < 0:
            j1 = 0
            j2 = j1 + 3
        if j2 > ny - 1:
            j2 = ny - 1
            j1 = j2 - 3
        mask = np.ones(map_data.shape, np.bool)
        mask[j, :] = 0
        mask[:, i] = 0
        k = mask[j1:j2, i1:i2].sum()
        mask[0, :] = 0
        mask[:, 0] = 0
        mask[ny-1, :]=0
        mask[:, nx-1]=0
        sub_mask = mask[j1:j2, i1:i2]
        sub_data = map_data[j1:j2, i1:i2]
        sub_data = sub_data[sub_mask]
        if calc_k:
            return k, sub_data
        else:
            return sub_data

    def add_points(self, data_map, ip, jp, i_PEAKS, j_PEAKS, flux_PEAKS, nx, ny,
                   seg_map, mask_map, n_reg, max_dist, min_flux, frac_peak):
        i_s = ip
        j_s = jp
        i = ip
        j = jp
        n_now = 0
        delta = 2
        sig = 0
        seg_map[j, i] = n_reg
        mask_map[j, i] = 0
        s=0
        stop_flag = True
        n_pixels = 0
        while stop_flag:
            if sig == 0:
                for jj in np.arange(j_s + 1, j_s + delta):
                    s = self.check_is_point(data_map=data_map, i=i, j=jj,
                                            ip=ip, jp=jp, nx=nx, ny=ny,
                                            i_PEAKS=i_PEAKS, j_PEAKS=j_PEAKS,
                                            n_reg=n_reg,
                                            seg_map=seg_map, mask_map=mask_map,
                                            min_flux=min_flux,
                                            max_dist=max_dist,
                                            frac_peak=frac_peak)
                    n_pixels = n_pixels + s
                    n_now = n_now + 1
                j = jj
                for ii in np.arange(i_s + 1, i_s + delta):
                    s = self.check_is_point(data_map=data_map, i=ii, j=j,
                                            ip=ip, jp=jp, nx=nx, ny=ny,
                                            i_PEAKS=i_PEAKS, j_PEAKS=j_PEAKS,
                                            n_reg=n_reg,
                                            seg_map=seg_map, mask_map=mask_map,
                                            min_flux=min_flux,
                                            max_dist=max_dist,
                                            frac_peak=frac_peak)
                    n_pixels = n_pixels + s
                    n_now = n_now + 1
                i = ii
                sig = 1
            else:
                for jj in np.arange(j_s - 1, j_s - delta, -1):
                    s = self.check_is_point(data_map=data_map, i=i, j=jj,
                                            ip=ip, jp=jp, nx=nx, ny=ny,
                                            i_PEAKS=i_PEAKS, j_PEAKS=j_PEAKS,
                                            n_reg=n_reg,
                                            seg_map=seg_map, mask_map=mask_map,
                                            min_flux=min_flux,
                                            max_dist=max_dist,
                                            frac_peak=frac_peak)
                    n_pixels = n_pixels + s
                    n_now = n_now + 1
                j = jj
                for ii in np.arange(i_s - 1, i_s - delta, -1):
                    s = self.check_is_point(data_map=data_map, i=ii, j=j,
                                            ip=ip, jp=jp, nx=nx, ny=ny,
                                            i_PEAKS=i_PEAKS, j_PEAKS=j_PEAKS,
                                            n_reg=n_reg,
                                            seg_map=seg_map, mask_map=mask_map,
                                            min_flux=min_flux,
                                            max_dist=max_dist,
                                            frac_peak=frac_peak)
                    n_pixels = n_pixels + s
                    n_now = n_now + 1
                i = ii
                sig = 0
            delta += 1
            j_s = j
            i_s = i
            dist = np.sqrt((i - ip)**2 + (j - jp)**2)
            if dist > 1.5 * self.dist_max_now:
                n_pix_max = 0.2 * self.max_area
                if n_pixels > n_pix_max:
                    done = True
                else:
                    done = False
                stop_flag = False
        return done

    def check_is_point(self, data_map, i, j, ip, jp, nx, ny, i_PEAKS, j_PEAKS,
                       n_reg, seg_map, mask_map, min_flux,
                       max_dist, frac_peak):
        dist = np.sqrt((i-ip)**2 + (j-jp)**2)
        is_point = 0
        F_PEAK_NOW = data_map[jp, ip]
        self.dist_max_now = max_dist + 0.5 * np.log(F_PEAK_NOW/1000)
        if self.dist_max_now < self.PSF:
            self.dist_max_now = self.PSF
        self.max_area = np.pi * self.dist_max_now**2
        if i > 0 and i < nx and j > 0 and j < ny:
            val = data_map[j, i]
            F_PEAK_func = data_map[jp, ip]
            if dist < max_dist and val > min_flux:
                if dist < self.dist_max_now and val > frac_peak * F_PEAK_func:
                    nr = seg_map[j, i]
                    if nr == 0:
                        seg_map[j, i] = n_reg
                        is_point = 1
                    else:
                        d_now = np.sqrt((i - ip)**2 + (j - jp)**2)
                        d_r = np.sqrt((i - i_PEAKS[int(nr-1)])**2
                                      + (j - j_PEAKS[int(nr-1)])**2)
                        if d_now < d_r:
                            seg_map[j, i] = n_reg
                            is_point = 1
                mask_map[j, i] = 0
        return is_point
        
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
