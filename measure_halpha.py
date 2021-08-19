# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:55:09 2021

@author: 55119
"""

import os
import shutil
import getpass
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.datasets import make_100gaussians_image
from photutils import CircularAperture, CircularAnnulus
import splusdata 
import context
import make_halpha_fornax

data_dir = os.path.join(context.data_dir, "FCC_halpha")
galaxies = sorted(os.listdir(data_dir))
username = input("Login for SPLUS cloud:") # Change to your S-PLUS usernam
password = getpass.getpass(f"Password for {username}:")

for galaxy in galaxies:
    # TODO: Adaptar da banda r para images halpha
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    filename =[x for x in os.listdir(wdir) if x.endswith("halpha.fits")][0]
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

    conn = splusdata.connect(username, password)
    ps = 0.55 * u.arcsec / u.pix
    size = 256 * u.pix
    r = np.sqrt(2) * size * ps # arcsec
    r = r.to(u.degree)

    tablename = os.path.join(context.home_dir,
                     "tables/Literature_new_phot_structural_parameters_8arcsec_class_star.fits")
    table = Table.read(tablename)
    for i, t in enumerate(table):
        ra0 = t["ALPHA_J2000"]
        dec0 = t["DELTA_J2000"]
        qtable = conn.query(f"""SELECT det.id, det.ra, det.dec 
                            FROM dr2.detection_image as det 
                            JOIN dr2_vacs.star_galaxy_quasar as sgq ON (sgq.ID = det.ID)
                            WHERE (sgq.PROB_STAR>0.8) 
                            AND 1=CONTAINS( POINT('ICRS', det.ra, det.dec), CIRCLE('ICRS', {ra0}, {dec0}, {r.to(u.degree).value}) )""")
        print(qtable)
        input()
########Pixel Masking##############################
    
    #hdul = fits.open(filename)
    # hdr = WCS(hdul[1].header)
    # print("hdr=",hdr)
    # sky = hdr.pixel_to_world(ra, dec)
    # print(sky)
    # 2) Fazer máscara das estrelas encontradas na query baseada no código
    # test_mascara.py
    # # ) Passar máscara como argumento do aperture_photometry
    ### Exemplo
   # mask = np.zeros_like(halpha)
    # Performing Aperture Photometry
    phot_table = aperture_photometry(halpha, apertures, mask=mask)
    # Lendo os valores da table
    phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in range(30)]
    table = Table([radii, phot], names=["sma", "halpha"])
    table.write("photometry_halpha.fits", overwrite=True)
    plt.plot(radii, phot, "o")
    plt.savefig('CUBE_FOTOMETRIA_Halpha.png')
    plt.show()