# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:55:09 2021

@author: 55119
"""

import os
import shutil

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.datasets import make_100gaussians_image
from photutils import CircularAperture, CircularAnnulus
import context

data_dir = os.path.join(context.data_dir, "FCC_halpha")
galaxies = sorted(os.listdir(data_dir))
for galaxy in galaxies:
    print(galaxy)
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    try:
        cubename = [x for x in os.listdir(wdir) if x.endswith("pix.fits")][0]
    except IndexError:
        shutil.rmtree(wdir)
        continue
    # Loading data
    data = fits.getdata(cubename, ext=1)
    t = Table.read(cubename)
    t['FILTER'].tolist()
    idx = t['FILTER'].tolist().index('R')
    rband = data[idx, :, : ]
    #Creating Aperture Objects
    positions =  np.array([0.5 * rband.shape[0], 0.5 * rband.shape[1]])
    radii = np.linspace(1, rband.shape[0] / 2., 30)
    # apertures = [CircularAperture(positions, r=r) for r in radii]
    apertures = []
    plt.subplot(1,2,1)
    vmin = np.percentile(rband, 10)
    vmax = np.percentile(rband, 95)
    plt.imshow(rband, vmin=vmin, vmax=vmax)
    for r in radii:
        aperture = CircularAperture(positions, r=r)
        apertures.append(aperture)
        aperture.plot(color='r', lw=1)
    plt.subplot(1,2,2)
    # Performing Aperture Photometry
    phot_table = aperture_photometry(rband, apertures)
    # Lendo os valores da table
    phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in range(30)]
    table = Table([radii, phot], names=["r", "photsum"])
    table.write("photometry_R.fits", overwrite=True)
    plt.plot(radii, phot, "o")
    plt.savefig('CUBE_FOTOMETRIA_R.png')
    plt.show()
    
    # #Multiple Apertures at Each Position
    # radii = [3., 4., 5.]
    # apertures = [CircularAperture(positions, r=r) for r in radii]
    # phot_table = aperture_photometry(data, apertures)
    # for col in phot_table.colnames:
    #     phot_table[col].info.format = '%.8g'

    #Sigma-clipped median within a circular annulus
    # data = make_100gaussians_image()
    # positions = [(145.1, 168.3), (84.5, 224.1), (48.3, 200.3)]
    # aperture = CircularAperture(positions, r=5)
    # annulus_aperture = CircularAnnulus(positions, r_in=10, r_out=15)
    # annulus_masks = annulus_aperture.to_mask(method='center')
    
    # plt.imshow(annulus_masks[0], interpolation='nearest')
    # plt.colorbar()