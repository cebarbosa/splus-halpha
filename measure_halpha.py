# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:55:09 2021
@author: 55119
"""

import os
import shutil
import getpass

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
import splusdata

import context
import make_halpha_fornax

data_dir = os.path.join(context.home_dir, "FCC_halpha")
galaxies = sorted(os.listdir(data_dir))
username = input("Login for SPLUS cloud:") # Change to your S-PLUS usernam
password = getpass.getpass(f"Password for {username}:")

for galaxy in galaxies:
    # TODO: Adaptar da banda r para images halpha
    wdir = os.path.join(data_dir, galaxy)
    os.chdir(wdir)
    filename =[x for x in os.listdir(wdir) if x.endswith("halpha.fits")][0]
    f = fits.open(filename)
    # Loading data
    halpha = fits.getdata(filename, ext=1)
    w = WCS(f[1].header)
    
    vmin = np.percentile(halpha, 10)
    vmax = np.percentile(halpha, 95)
    img = plt.imshow(halpha, origin="lower", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.show()
    
    #Creating Aperture Objects
    positions =  np.array([0.5 * halpha.shape[0], 0.5 * halpha.shape[1]])
    radii = np.linspace(1, halpha.shape[0] / 2., 30)
    # apertures = [CircularAperture(positions, r=r) for r in radii]
    apertures = []
    plt.plot(1,2,1)
    
    for r in radii:
        aperture = CircularAperture(positions, r=r)
        apertures.append(aperture)
        aperture.plot(color='r', lw=1)
        plt.plot(1,2,2)
        
    # TODO: FAzer máscara para estrelas
    # 1) Fazer query na região da imagem para achar as estrela próximas

    conn = splusdata.connect(username, password)
    ps = 0.55 * u.arcsec / u.pix
    size = 256 * u.pix
    r = np.sqrt(2) / 2 * size * ps # arcsec
    r = (r.to(u.degree).value)

    tablename = os.path.join(context.home_dir,
                     "tables/Literature_new_phot_structural_parameters_8arcsec_class_star.fits")
    table = Table.read(tablename)
    
    for i, t in enumerate(table):
        ra0 = t["ALPHA_J2000"]
        dec0 = t["DELTA_J2000"]
        qtable = conn.query(f"""SELECT det.ID, det.ra, det.dec 
                 FROM idr3.detection_image as det  
                 JOIN idr3_vacs.star_galaxy_quasar as sgq ON (sgq.ID = det.ID)
                 WHERE (sgq.PROB_STAR>0.8) AND 1=CONTAINS( POINT('ICRS', det.ra, det.dec), CIRCLE('ICRS', {ra0}, {dec0}, {r}) )""")

        #for i in Result:
        ra = qtable["RA"].data * u.degree
        dec = qtable["DEC"].data * u.degree
        if len(ra) > 0:
            coord = SkyCoord(ra, dec)
            xpix, ypix = w.world_to_pixel(coord)
        positions = np.transpose((qtable["RA"], qtable["DEC"]))
        print(positions)
        apertures = CircularAperture(positions, r=r)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(halpha, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
    # #mean, median, std = sigma_clipped_stats(halpha, sigma=3.0)       
    # #daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) 
    # sources = daofind(halpha - median) 
    # daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) 
    # sources = daofind(halpha - median)
    
    #     mask = np.zeros_like(halpha)
    #     print(mask)
    #     input()
      
    # masked_data = halpha[:]
    # masked_data[mask] = median
    # plt.imshow(masked_data, origin="lower", vmax=vmax, vmin=vmin)
    # plt.colorbar()
    # plt.show()
########Pixel Masking##############################
    

    # 2) Fazer máscara das estrelas encontradas na query baseada no código
    # test_mascara.py
    # # ) Passar máscara como argumento do aperture_photometry
    ### Exemplo
   # mask = np.zeros_like(halpha)
    # Performing Aperture Photometry
    #phot_table = aperture_photometry(halpha, apertures, mask=mask)
    # Lendo os valores da table
   # phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in range(30)]
    # table = Table([radii, phot], names=["sma", "halpha"])
    # table.write("photometry_halpha.fits", overwrite=True)
    # plt.plot(radii, phot, "o")
    # plt.savefig('CUBE_FOTOMETRIA_Halpha.png')
    # plt.show()