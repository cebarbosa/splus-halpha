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


"[N ii] correction"

def get_names(wdir):
    filenames = [x for x in os.listdir(wdir) if x.endswith("_swp.fits")]
    names = []
    for band in ["R","I","F660","G"]:
        names.append([x for x in filenames if x.split("_")[2]==band][0])
    return names

coef_F660 = 125.3
coef_r = 1419.0

def read_galaxies(image):
    
    data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
    galaxies = os.listdir(data_dir)
    
    for galaxy in galaxies:
        wdir = os.path.join(data_dir, galaxy)
        os.chdir(wdir)
        imgnames = get_names(wdir)
        #image=fits.open(imgnames)
        # Loading data
        data = np.array([fits.getdata(name, ext=1) for name in imgnames])
        # Loading zero points
        zero_point = np.array([fits.getval(name, "MAGZP", ext=1) for name in imgnames])
        data_galaxies= fits.read(data)
        fnu =  data * np.power(10, -0.4 * (zero_point[:, None, None] + 48.6))
    return fnu


def three_filters(fnu):
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i= fnu[2]
    fnu_F660_corr= 2* fnu_F660
    
    coef_file = os.path.join(context.tables_dir, "coeffs.fits")
    COEFS=fits.open("file:///c:/Users/amori/Dropbox/splus-halpha (1)/tables/coeffs.fits")

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

def corretion_nii(fnu_F660_corr):
    log_halpha = np.where(("G" - "I") <=0.5,
                          0.989 * np.log10(fnu_F660_corr.value)-0.193,
                          0.954 * np.log10(fnu_F660_corr.value)-0.193)
    #output_name = calcular_correçao_nii
    #correcao= log_halpha.read()
    return log_halpha

def save_the_disk(correcao,output_name):
    save_the_disk
    
def plot_corretion(log_halpha):
    vmax = np.nanpercentile(log_halpha, 95)
    vmin = np.nanpercentile(log_halpha, 80)
    plt.imshow(log_halpha, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plote= plt.show()
    return plote
    
   
