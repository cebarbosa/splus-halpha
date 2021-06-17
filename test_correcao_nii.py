# -*- coding: utf-8 -*-
""" 
Author : Carlos Eduardo Barbosa
Co-author: Jessica Silva Amorim

Methods to extract H-alpha and 
study the Star Formation Rate 
with SPLUS filters according 
to methods presented 
in Vilella-Rojo (2015) 

"""

import os


import matplotlib.pyplot as plt

import extinction
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

import context
from test_halpha import get_names


def calc_halpha_nii_3F(fnu):
    """ Calculates Halpha with 3 bands using eq. 3 of Vilella-Rojo+ (2015). """
    fnu_F660 = fnu[1]
    fnu_r = fnu[0]
    fnu_i = fnu[2]
    fnu_F660_corr = 2 * fnu_F660

    coef_file = os.path.join(context.tables_dir, "coeffs.fits")
    COEFS = fits.open(coef_file)

    t = Table.read(coef_file)
    alpha_F660 = t['alpha_x'][0]
    beta_F660 = t['beta_x'][0]
    delta_F660 = t['delta_x'][0]

    alpha_i = t['alpha_x'][1]
    beta_i = t['beta_x'][1]
    delta_i = i = t['delta_x'][1]

    alpha_r = t['alpha_x'][2]
    beta_r = t['beta_x'][2]
    delta_r = t['delta_x'][2]

    a = (alpha_r - alpha_i) / (alpha_F660 - alpha_i)
    numen = (fnu_r - fnu_i) - a * (fnu_F660 - fnu_i)
    denom = - a * beta_F660 + beta_r
    return numen / denom


def calc_halpha_without_nii(halpha_nii, g_i): 
    """ Calculated NII corrected emission with eq. 21 or Vilella-Rojo+ (2015)"""
    idx = np.where(halpha_nii < 0)
    halpha = np.power(10, np.where((g_i <= 0.5),
                          0.989 * np.log10(halpha_nii) - 0.193,
                          0.954 * np.log10(halpha_nii) - 0.753))
    halpha[idx] = halpha_nii[idx]
    corr = halpha / halpha_nii
    return halpha

def process_galaxies():
    data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
    galaxies = os.listdir(data_dir)
    for galaxy in galaxies:
        wdir = os.path.join(data_dir, galaxy)
        os.chdir(wdir)
        imgnames = get_names(wdir, bands=["R","F660","I","G"])
        # Loading data
        data = np.array([fits.getdata(name, ext = 1) for name in imgnames])
        # Loading zero points
        zero_point = np.array(
            [fits.getval(name, "MAGZP", ext = 1) for name in imgnames])
        fnu = data * np.power(10, - 0.4 * (zero_point[:, None, None] + 48.6))
        halpha_nii = calc_halpha_nii_3F(fnu)
        magAB = -2.5 * np.log10(fnu) - 48.6
        g_i = magAB[3] - magAB[2]
        wave = 6614.0 * u.Angstrom
        """" C = E(B-V) extinction law  with  eq. 20 or Vilella-Rojo+ (2015)"""
        ebv = np.array(0.206 * np.power(g_i, 1.68) - 0.0457)
        ebv[np.isnan(ebv)] = 0
        x = 1 / wave.to(u.micrometer).value
        wtran = (0.63 * u.micrometer)
        kappa = np.where(wave > wtran, 2.659 * (-1.857 + 1.040 * x),
                               2.659 * (-2.156 + 1.509 * x - 0.198 * x * x
                                       + 0.011 * (x * x * x)))
        """" Extinction.calzetti00 with A_V = R_V * c and R_V = 4.05 +-0.80 """
        dust_correction = np.power(10, -0.4 * ebv * kappa)
        halpha_nii_corr = halpha_nii * dust_correction
        halpha = calc_halpha_without_nii(halpha_nii_corr, g_i)
        
        vmax = np.percentile(halpha_nii, 95)
        vmin = np.percentile(halpha_nii, 10)
        plt.subplot(1,3,1)
        plt.imshow(halpha_nii, vmax=vmax, vmin=vmin, origin="lower")
        plt.colorbar()
        plt.subplot(1,3,2)
        plt.imshow(halpha_nii_corr, vmax=vmax, vmin=vmin, origin="lower")
        plt.colorbar()
        plt.subplot(1,3,3)
        
        vmax = np.percentile(halpha, 95)
        vmin = np.percentile(halpha, 10)
        plt.imshow(halpha, vmin=vmin, vmax=vmax, origin="lower")
        plt.colorbar()
        plt.show()


if __name__ == "__main__":
    process_galaxies()