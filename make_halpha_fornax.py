import os
import getpass 
import sys # For authentication
import splusdata
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt # To access the S-PLUS database
from splus_ifusci import SCube, make_RGB_with_overlay

import context
from halpha_3filters_corrections import dust_correction, nii_correction

if __name__ == "__main__":
    wdir = os.path.join(context.home_dir, "FCC_halpha")
    if not os.path.exists(wdir):
        os.mkdir(wdir)
    tablename = os.path.join(context.home_dir,
                         "tables/Literature_new_phot_structural_parameters_8arcsec_class_star.fits")
    table = Table.read(tablename)
    size = 256
    # Connect with S-PLUS
    username = input("Login: ")  # Change to your S-PLUS username
    password = getpass.getpass(f"Password:")
    conn = splusdata.connect(username, password)
    for i, t in enumerate(table):
        ra = t["ALPHA_J2000"]
        dec = t["DELTA_J2000"]
        coords = [ra, dec]
        galaxy = str(f"""Fornax{t["NUMBER"]}""")
        gal_dir = os.path.join(wdir, galaxy)
        if not os.path.exists(gal_dir):
            os.mkdir(gal_dir)
        halpha_img = os.path.join(gal_dir, f"{galaxy}_{size}x"
                                        f"{size}pix_halpha.fits")
        # Main routine to download the datacube.
        scube = SCube(galaxy, coords, size, conn=conn,
                      coord_unit=(u.degree, u.degree), wdir=gal_dir)
        scube.download_stamps()
        scube.make_cube()
        # Processing the data
        mag = scube.get_mag()
        magerr = scube.get_magerr()
        idx_g = scube.bands.index("G")
        idx_i = scube.bands.index("I")
        g_i = mag[idx_g] - mag[idx_i]
        halpha_obs, halpha_obs_err = scube.calc_halpha(store=True)
        halpha_nii = dust_correction(halpha_obs.value, g_i)
        halpha_nii_err = dust_correction(halpha_obs_err.value, g_i)
        # halpha = nii_correction(halpha_nii, g_i)
        # halpha_err = nii_correction(halpha_nii_err, g_i)
        # Saving fits
        output = scube.cubename.replace(".fits", "_halpha.fits")
        h = fits.getheader(os.path.join(scube.cutouts_dir, scube.cutnames[0]),
                           ext=1)
        h["EXTNAME"] = "DATA"
        hdu1 = fits.ImageHDU(halpha_nii, h)
        h["EXTNAME"] = "ERROR"
        hdu2 = fits.ImageHDU(halpha_nii_err, h)
        hdulist = fits.HDUList([fits.PrimaryHDU(), hdu1, hdu2])
        hdulist.writeto(output, overwrite=True)
        flam = scube.get_flam().value
        rgb_bands = ["I", "R", "G"]
        rgb = [flam[scube.bands.index(b)] for b in rgb_bands]
        outimg = os.path.join(gal_dir, f"{galaxy}_RGB.png")
        make_RGB_with_overlay(*rgb, outimg, overlay=halpha_obs)
