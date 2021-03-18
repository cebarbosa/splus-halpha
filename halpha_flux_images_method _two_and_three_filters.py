# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 16:58:48 2021

@author: 55119
"""

import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import splusdata
import context


def get_names(wdir):
    
# wdir=pasta+nome da galaxia

    """ Obtaining names of images  for one directory. """
    filenames = [x for x in os.listdir(wdir) if x.endswith("_swp.fits")]
    names = []
    for band in ["R","I","F660"]:
        names.append([x for x in filenames if x.split("_")[2]==band][0])
    return names


data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
galaxies = os.listdir(data_dir)

coef_F660 = 125.3
coef_r = 1419.0

for galaxy in galaxies:
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    imgnames = get_names(wdir)
    # Loading data
    data = np.array([fits.getdata(name, ext=1) for name in imgnames])
    # Loading zero points
    m0 = np.array([fits.getval(name, "MAGZP", ext=1) for name in imgnames])
    fnu =  data * np.power(10, -0.4 * (m0[:, None, None] + 48.6))
    
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i= fnu[2]
    flux_two_bands = np.dot(coef_F660,((fnu_F660 - fnu_r)/(1-(coef_F660/coef_r))))
    
    vmax = np.nanpercentile(flux_two_bands, 95)
    vmin = np.nanpercentile(flux_two_bands, 80)
    plt.imshow(flux_two_bands, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()
    
             #############   For 3 Filters    ####################
    coef_file = os.path.join(context.tables_dir, "coeffs.fits")

# For 3 Filters    
    
    COEFS=fits.open("file:///c:/Users/55119/Dropbox/splus-halpha (1)/tables/coeffs.fits")

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

    flux_three_bands = (((fnu_r- fnu_i)-((alpha_r-alpha_i)/(alpha_F660 - alpha_i))*(fnu_F660 - fnu_i))/((beta_F660)*(alpha_i -alpha_r )-(beta_r)))
    
   # fluxThreeBands = (((fnu_R- fnu_I)-((t['alpha_x'][2])-(t['alpha_x'][1]))/((t['alpha_x'][0] )- (t['alpha_x'][2]))*(fnu_F660 - fnu_I))/(((t['beta_x'][0]))*((t['alpha_x'][1]) - (t['alpha_x'][2]) )-((t['beta_x'][2))))
    vmax = np.nanpercentile(flux_three_bands, 95)
    vmin = np.nanpercentile(flux_three_bands, 80)
    plt.imshow(flux_three_bands, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()
    
 
    