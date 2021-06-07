# -*- coding: utf-8 -*-
"""
Created on Wed May 12 00:01:03 2021

@author: amori
"""

import os

import numpy as np
import matplotlib.pyplot as plt
import extinction
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

import context
from test_halpha import get_names

def calc_halpha_nii_3F(fnu):
    """ Calculates Halpha with 3 bands using eq. 3 of Vilella- Rojo (2015)."""
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i= fnu[2]
    fnu_F660_corr= 2* fnu_F660
    
    coef_file = os.path.join(context.tables_dir, "coeffs.fits")
    COEFS=fits.open(coef_file)

    t = Table.read(coef_file)
    alpha_F660=t['alpha_x'][0]
    beta_F660=t['beta_x'][0]
    #delta_F660=t['delta_x'][0]

    alpha_i=t['alpha_x'][1]
    beta_i=t['beta_x'][1]
    #delta_i=i=t['delta_x'][1]

    alpha_r=t['alpha_x'][2]
    beta_r=t['beta_x'][2]
    #delta_r=t['delta_x'][2]
    
    a =(alpha_r- alpha_i) / (alpha_F660 - alpha_i)
    numen = (fnu_r - fnu_i) - a * (fnu_F660 - fnu_i)
    denon = - a * beta_F660 + beta_r
    return numen / denon 

def calc_halpha_corrected(halpha_nii, g_i):
    """" Calculated NII corrected emission with eq. 21 or Vilella-Rojo+ (2015)"""
    y =halpha_nii
    vmin = np.nanpercentile(y, 10)
    vmax = np.nanpercentile(y, 90)
    plt.imshow(y, vmin=vmin, vmax=vmax, origin="lower")
    plt.show()
    log_halpha = np.where(g_i <= 0.5,
                          0.989 * np.power(10, halpha_nii)- 0.193,
                          0.954 * np.power(10, halpha_nii)- 0.193)
    return np.power(10, log_halpha)
# def save_the_disk(correcao,output_name):
#     save_the_disk
#
# def plot_corretion(log_halpha, g_i):
#     vmax = np.nanpercentile(log_halpha, 95)
#     vmin = np.nanpercentile(log_halpha, 80)
#     return plote

def process_galaxies():
    data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
    galaxies = os.listdir(data_dir)
    for galaxy in galaxies:
        wdir = os.path.join(data_dir, galaxy)
        os.chdir(wdir)
        imgnames = get_names(wdir, bands=["R","F660","I","G"])
        # image=fits.open(imgnames)
        # Loading data
        data = np.array([fits.getdata(name, ext=1) for name in imgnames])
        # Loading zero points
        zero_point = np.array(
            [fits.getval(name, "MAGZP", ext=1) for name in imgnames])
        fnu = data * np.power(10, -0.4 * (zero_point[:, None, None] + 48.6))
        halpha_nii = calc_halpha_nii_3F(fnu)
        ######################################################################
        # dust correction
        """extinction.calzetti00 with with A_V = 3.8 and R_V = 4.05 +-0.80"""
        wave = np.array([6266.6, 6614.0, 7683.8])
        dust_correction = extinction.calzetti00(wave, 3.8, 4.05)
        print ("dust_correction",dust_correction)
        ######################################################################
        #""" astropy (wcs)"""
        #w = WCS(data[1].header)
        #xdim, ydim = data.shape
        #x0 = xdim / 2.
        #y0 = ydim / 2.
        #sky = w.pixel_to_world(x0, y0)
        #xdim, ydim = data.shape
        #x0, y0 = w.world_to_pixel(sky)
        #print ("x0, y0=", x0, y0)
         
        magAB = -2.5 * np.log10(fnu) - 48.6
        g_i = magAB[3] - magAB[2]
        
        halpha = calc_halpha_corrected(halpha_nii, g_i)
        
        vmin = np.nanpercentile(halpha, 10)
        vmax = np.nanpercentile(halpha, 90)
        plt.imshow(halpha, vmin=vmin, vmax=vmax, origin="lower")
        plt.show()
        
        condition = np.where(g_i < 0.5, 1, 0)
        plt.imshow(condition)
        plt.show()
        # plt.imshow(g_i, origin="lower")
        # plt.show()
        # rcut = np.power(10, -0.4 * 0.5)
        # condition = np.where(r <= rcut, 1, 0)
        # fig = plt.figure(5)
        # plt.imshow(r, origin="lower")
        # plt.colorbar()
        # plt.show()


if __name__ == "__main__":
    process_galaxies()