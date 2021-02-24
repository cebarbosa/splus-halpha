# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:55:09 2021

@author: 55119
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import context
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.datasets import make_100gaussians_image
from photutils import CircularAperture, CircularAnnulus

def get_names(wdir):
# wdir=pasta+nome da galaxia
    """ Obtaining names of images  for one directory. """
    filenames = [x for x in os.listdir(wdir) if x.endswith("_swp.fits")]
    names = []
    for band in ["R"]:
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
    #Loading pixel coordinate of reference point
    x = np.array([fits.getval(name, "CRPIX1", ext=1) for name in imgnames])
    y = np.array([fits.getval(name, "CRPIX2", ext=1) for name in imgnames])
    
    #Creating Aperture Objects
                
    positions = []
    positions = np.append(x,y)
    
    aperture = CircularAperture(positions, r=3.)
   
    #Performing Aperture Photometry
    dat= np.ones((100, 100))
    phot_table = aperture_photometry(dat, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'
    print(phot_table)
    print(np.pi * 3. ** 2)
     
    #Multiple Apertures at Each Position
    radii = [3., 4., 5.]
    apertures = [CircularAperture(positions, r=r) for r in radii]
    phot_table = aperture_photometry(data, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'
    print(phot_table)    
   
    #Sigma-clipped median within a circular annulus
    data = make_100gaussians_image()
    positions = [(145.1, 168.3), (84.5, 224.1), (48.3, 200.3)]
    aperture = CircularAperture(positions, r=5)
    annulus_aperture = CircularAnnulus(positions, r_in=10, r_out=15)
    annulus_masks = annulus_aperture.to_mask(method='center')
    
    plt.imshow(annulus_masks[0], interpolation='nearest')
    plt.colorbar()