# -*- coding: utf-8 -*-
"""
Created on Fri May  7 12:24:36 2021

@author: amori
"""
import os
import sys
import glob

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt

import context

def three_filters(flux_three_bands, fnu):  
    
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i= fnu[2]
    
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

    flux_three_bands = (((fnu_r- fnu_i)-((alpha_r-alpha_i)/(alpha_F660 -alpha_i))*(fnu_F660 - fnu_i))/((beta_F660)*(alpha_i -alpha_r )-(beta_r)))
    return flux_three_bands
    #log_halpha = np.where(g_i <=0.5,
                        #  0.989 * np.log10(fnu_F660_corr.value)-0.193,
                        #  0.954 * np.log10(fnu_F660_corr.value)-0.193)

   # fluxThreeBands = (((fnu_R- fnu_I)-((t['alpha_x'][2])-(t['alpha_x'][1]))/((t['alpha_x'][0] )- (t['alpha_x'][2]))*
   # (fnu_F660 - fnu_I))/(((t['beta_x'][0]))*((t['alpha_x'][1]) - (t['alpha_x'][2]) )-((t['beta_x'][2))))
    
def nii_correction(log_halpha,fnu_F660_corr):
    g_i=[]
    fnu_F660_corr= 2* fnu_F660
    log_halpha = np.where(g_i <=0.5,
                          0.989 * np.log10(fnu_F660_corr0.193,
                          0.954 * np.log10(fnu_F660_corr)-0.193)
    print log_halpha

    vmax = np.nanpercentile(flux_three_bands, 95)
    vmin = np.nanpercentile(flux_three_bands, 80)
    plt.imshow(flux_three_bands, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()