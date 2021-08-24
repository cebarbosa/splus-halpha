# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:55:09 2021
@author: 55119
"""

import os
import getpass

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils import aperture_photometry
from photutils import CircularAperture
import splusdata

import context

def make_mask(data, xpix, ypix, rmask=15):
    """ Produces mask for h-alpha image photometry. """
    mask = np.zeros_like(halpha).astype(np.int)
    xdim, ydim = data.shape
    x = np.arange(xdim)
    y = np.arange(ydim)
    xx, yy = np.meshgrid(x, y)
    for x0, y0 in zip(xpix, ypix):
        r = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)
        mask[r<=rmask] = 1
    return mask

def make_halpha_image(halpha, galaxy):
    """ Produced PNG image for halpha. """
    vmin = np.percentile(halpha, 10)
    vmax = np.percentile(halpha, 99)
    plt.imshow(halpha, origin="lower", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.xlabel("X (pix)")
    plt.ylabel("Y (pix)")
    plt.savefig(f"{galaxy}_halpha.png", dpi=250)
    plt.close()

if __name__ == "__main__":
    data_dir = os.path.join(context.home_dir, "FCC_halpha")
    galaxies = sorted(os.listdir(data_dir))
    conn = None
    ############################################################################
    # Defining radius for query search
    ps = 0.55 * u.arcsec / u.pix
    size = 256 * u.pix
    r = np.sqrt(2) / 2 * size * ps # arcsec
    r = (r.to(u.degree).value)
    ############################################################################
    tablename = os.path.join(context.home_dir,
                     "tables/Literature_new_phot_structural_parameters_8arcsec_class_star.fits")
    table = Table.read(tablename)
    for i, t in enumerate(table):
        galaxy = f"Fornax{t['NUMBER']}"
        wdir = os.path.join(data_dir, galaxy)
        os.chdir(wdir)
        filename = [x for x in os.listdir(wdir) if x.endswith("halpha.fits")][0]
        f = fits.open(filename)
        # Loading data
        halpha = fits.getdata(filename, ext=1)
        w = WCS(f[1].header)
        ra0 = t["ALPHA_J2000"]
        dec0 = t["DELTA_J2000"]
        qtablefile = "query_idr3.fits"
        if not os.path.exists(qtablefile):
            if conn is None:
                username = input("Login for SPLUS cloud:") # Change to your S-PLUS usernam
                password = getpass.getpass(f"Password for {username}:")
                conn = splusdata.connect(username, password)
            qtable = conn.query(f"""SELECT det.ID, det.ra, det.dec 
                     FROM idr3.detection_image as det  
                     JOIN idr3_vacs.star_galaxy_quasar as sgq ON (sgq.ID = det.ID)
                     WHERE (sgq.PROB_STAR>0.8) AND 1=CONTAINS( POINT('ICRS', det.ra, det.dec), CIRCLE('ICRS', {ra0}, {dec0}, {r}) )""")
            qtable.write(qtablefile)
        else:
            qtable = Table.read(qtablefile)
        #for i in Result:
        ra = qtable["RA"].data * u.degree
        dec = qtable["DEC"].data * u.degree
        coord = SkyCoord(ra, dec)
        xpix, ypix = w.world_to_pixel(coord)
        mask = make_mask(halpha, xpix, ypix)
        halpha[mask==1] = 0
        make_halpha_image(halpha, galaxy)
        # Performing Aperture Photometry
        # Creating Aperture Objects
        positions = np.array([0.5 * halpha.shape[0], 0.5 * halpha.shape[1]])
        radii = np.linspace(1, halpha.shape[0] / 2., 30)
        # apertures = [CircularAperture(positions, r=r) for r in radii]
        apertures = []
        plt.subplot(1, 2, 1)
        vmin = np.percentile(halpha, 10)
        vmax = np.percentile(halpha, 95)
        plt.imshow(halpha, vmin=vmin, vmax=vmax)
        for r in radii:
            aperture = CircularAperture(positions, r=r)
            apertures.append(aperture)
            aperture.plot(color='r', lw=1)
        plt.subplot(1, 2, 2)
        phot_table = aperture_photometry(halpha, apertures, mask=mask)
        # Lendo os valores da tabela
        phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in
                range(30)]
        table = Table([radii, phot], names=["sma", "halpha"])
        table.write("photometry_halpha.fits", overwrite=True)
        plt.plot(radii, phot, "o")
        plt.savefig('halpha_photometry.png')
        plt.show()