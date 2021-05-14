# -*- coding: utf-8 -*-
"""
Created on Wed May 12 00:01:03 2021

@author: amori
"""

import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Table

import context
from test_halpha import get_names

"[N ii] correction"

coef_F660 = 125.3
coef_r = 1419.0


def three_filters(fnu):
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i= fnu[2]
    fnu_F660_corr= 2* fnu_F660
    
    coef_file = os.path.join(context.tables_dir, "coeffs.fits")
    COEFS=fits.open(coef_file)

    t = Table.read(coef_file)
    alpha_F660=t['alpha_x'][0]
    beta_F660=t['beta_x'][0]
    delta_F660=t['delta_x'][0]

    alpha_i=t['alpha_x'][1]
    beta_i=t['beta_x'][1]
    delta_i=i=t['delta_x'][1]

    alpha_r=t['alpha_x'][2]
    beta_r=t['beta_x'][2]
    delta_r=t['delta_x'][2]

    flux_three_bands = (((fnu_r- fnu_i)-((alpha_r-alpha_i)/(alpha_F660 -alpha_i))*
                         (fnu_F660 - fnu_i))/((beta_F660)*(alpha_i -alpha_r )-(beta_r)))
    
    return flux_three_bands

def corretion_nii(halpha_nii, g_i):
    log_halpha = np.where(g_i <= 0.5,
                          0.989 * np.log10(halpha_nii)-0.193,
                          0.954 * np.log10(halpha_nii)-0.193)
    return log_halpha

def save_the_disk(correcao,output_name):
    save_the_disk
    
def plot_corretion(log_halpha, g_i):
    vmax = np.nanpercentile(log_halpha, 95)
    vmin = np.nanpercentile(log_halpha, 80)
    return plote

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
        halpha_nii = three_filters(fnu)
        magAB = -2.5 * np.log10(fnu) - 48.6
        r = np.power(10, -0.4 * (fnu[3]/fnu[2]))
        rcut = np.power(10, -0.4 * 0.5)
        condition = np.where(r <= rcut, 1, 0)
        plt.imshow(condition, origin="lower")
        plt.colorbar()
        plt.show()


if __name__ == "__main__":
    process_galaxies()