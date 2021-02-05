# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 16:58:48 2021

@author: 55119
"""

import glob, os, sys, numpy as np, matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import numpy as np
import splusdata

import context

def get_names(wdir):
# wdir=pasta+nome da galaxia
    """ Obtaining names of images  for one directory. """
    filenames = [x for x in os.listdir(wdir) if x.endswith("_swp.fits")]
    names = []
    for band in ["R", "F660", "I"]:
        names.append([x for x in filenames if x.split("_")[2]==band][0])
    return names


data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
galaxies = os.listdir(data_dir)

coef_F660 = 125.3
coef_R = 1419.0

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
    fnu_R = fnu[0]
    fnu_I= fnu[2]
    fluxTwoBands = np.dot(coef_F660,((fnu_F660 - fnu_R)/(1-(coef_F660/coef_R))))
    
    vmax = np.nanpercentile(fluxTwoBands, 95)
    vmin = np.nanpercentile(fluxTwoBands, 80)
    plt.imshow(fluxTwoBands, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()
    
             #############   For 3 Filters    ####################
    
    COEFS=fits.open("file:///c:/Users/55119/Dropbox/splus-halpha (1)/tables/coeffs.fits")
   
    t = Table.read(COEFS)
    
    alpha_F660=t['alpha_x'][0]
    beta_F660=t['beta_x'][0]
    delta_F660=t['delta_x'][0]

    alpha_I=t['alpha_x'][1]
    beta_I=t['beta_x'][1]
    delta_I=t['delta_x'][1]

    alpha_R=t['alpha_x'][2]
    beta_R=t['beta_x'][2]
    delta_R=t['delta_x'][2]

    fluxThreeBands = (((fnu_R- fnu_I)-((alpha_R-alpha_I)/(alpha_F660 - alpha_I))*(fnu_F660 - fnu_I))/((beta_F660)*(alpha_I -alpha_R )-(beta_R)))
    
   # fluxThreeBands = (((fnu_R- fnu_I)-((t['alpha_x'][2])-(t['alpha_x'][1]))/((t['alpha_x'][0] )- (t['alpha_x'][2]))*(fnu_F660 - fnu_I))/(((t['beta_x'][0]))*((t['alpha_x'][1]) - (t['alpha_x'][2]) )-((t['beta_x'][2))))
    vmax = np.nanpercentile(fluxThreeBands, 95)
    vmin = np.nanpercentile(fluxThreeBands, 80)
    plt.imshow(fluxThreeBands, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()
    
    