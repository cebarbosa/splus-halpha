# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 02:25:38 2021

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
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder

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

for galaxy in galaxies:
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    imgnames = get_names(wdir)
    # Loading data
    data = fits.getdata(imgnames[0], ext=1)
       
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    print(mean, median, std)
    #print((mean, median, std)) 
    
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) 
    sources = daofind(data - median) 
    #print(sources)
    xdim, ydim = data.shape
    x0 = xdim / 2.
    y0 = ydim / 2.
    
    r = np.sqrt((sources["xcentroid"]- x0)**2 + \
            (sources["ycentroid"] - y0)**2).data
   # print(r)
    idx = np.where(r >= 100)
    
    stars = sources[idx]
    
    img = plt.imshow(data - median, origin="lower")
    # for star in stars:
    #     plt.plot(star["xcentroid"], star["ycentroid"], "xr")
    # plt.colorbar(img)
    # plt.show()
    
    x = np.arange(xdim)
    y = np.arange(ydim)
    xx, yy = np.meshgrid(x, y)
    
    mask = np.zeros_like(data).astype(np.bool)
    rstars = 15
    for star in stars:
        r = np.sqrt((xx - star["xcentroid"])**2 + \
               (yy - star["ycentroid"])**2)
        idx = np.where(r < rstars)
        mask[idx] = True
    # plt.imshow(mask, origin="lower")
    # plt.colorbar()
    # plt.show()
    
    masked_data = data[:] - median
    plt.imshow(masked_data, origin="lower",
               vmax=np.percentile(masked_data, 99.9))
    plt.colorbar()
    plt.show()

 #Creating Aperture Objects
    positions =  np.array([0.5 * data.shape[0], 0.5 * data.shape[1]])
    radii = np.linspace(1, data.shape[0] / 2., 30)
    # apertures = [CircularAperture(positions, r=r) for r in radii]
    apertures = []
    for r in radii:
        aperture = CircularAperture(positions, r=r)
        apertures.append(aperture)
    # Performing Aperture Photometry
    phot_table = aperture_photometry(masked_data, apertures, mask=mask)
    # Lendo os valores da table
    phot = [-2.5 * np.log10(float(phot_table["aperture_sum_{}".format(i)]))
            for i in range(30)]
    plt.plot(radii, phot, "o")
    plt.show()
