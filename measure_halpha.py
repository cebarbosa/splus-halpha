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
    # TODO: Adaptar da banda r para images halpha
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    # Loading data
    halpha = fits.getdata(filename, ext=1)
    #Creating Aperture Objects
    positions =  np.array([0.5 * halpha.shape[0], 0.5 * halpha.shape[1]])
    radii = np.linspace(1, halpha.shape[0] / 2., 30)
    # apertures = [CircularAperture(positions, r=r) for r in radii]
    apertures = []
    plt.subplot(1,2,1)
    vmin = np.percentile(halpha, 10)
    vmax = np.percentile(halpha, 95)
    plt.imshow(halpha, vmin=vmin, vmax=vmax)
    for r in radii:
        aperture = CircularAperture(positions, r=r)
        apertures.append(aperture)
        aperture.plot(color='r', lw=1)
    plt.subplot(1,2,2)
    # TODO: FAzer máscara para estrelas
    # 1) Fazer query na região da imagem para achar as estrela próximas
    # 2) Fazer máscara das estrelas encontradas na query baseada no código
    # test_mascara.py
    # # ) Passar máscara como argumento do aperture_photometry
    ### Exemplo
    mask = np.zeros_like(halpha)
    # Performing Aperture Photometry
    phot_table = aperture_photometry(halpha, apertures, mask=mask)
    # Lendo os valores da table
    phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in range(30)]
    table = Table([radii, phot], names=["sma", "halpha"])
    table.write("photometry_halpha.fits", overwrite=True)
    plt.plot(radii, phot, "o")
    plt.savefig('CUBE_FOTOMETRIA_Halpha.png')
    plt.show()