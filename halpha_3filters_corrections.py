# -*- coding: utf-8 -*-
""" 
Author : Carlos Eduardo Barbosa
Co-author: Jessica Silva Amorim

Methods to extract H-alpha and study the Star Formation Rate with SPLUS filters
according to methods presented in Vilella-Rojo (2015)

"""
import os

import astropy.units as u
from astropy.io import fits
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
    dir_path = os.path.dirname(os.path.realpath(__file__))
    coef_file = os.path.join(dir_path, "coeffs.csv")
    t = Table.read(coef_file)
    alpha_F660 = t['alpha_x'][0]
    beta_F660 = t['beta_x'][0]
    alpha_i = t['alpha_x'][1]
    alpha_r = t['alpha_x'][2]
    beta_r = t['beta_x'][2]
    a = (alpha_r - alpha_i) / (alpha_F660 - alpha_i)
    numen = (fnu_r - fnu_i) - a * (fnu_F660 - fnu_i)
    denom = - a * beta_F660 + beta_r
    return numen / denom

def nii_correction(halpha_nii, g_i):
    """ Calculated NII corrected emission with eq. 21 or Vilella-Rojo+ (2015)"""
    idx = np.where(halpha_nii < 0)
    halpha = np.power(10, np.where((g_i <= 0.5),
                          0.989 * np.log10(halpha_nii) - 0.193,
                          0.954 * np.log10(halpha_nii) - 0.753))
    halpha[idx] = halpha_nii[idx]
    return halpha

def dust_correction(halpha, g_i):
    wave = 6614.0 * u.Angstrom
    """" C = E(B-V) extinction law  with  eq. 20 or Vilella-Rojo+ (2015)"""
    ebv = np.array(0.206 * np.power(g_i, 1.68) - 0.0457)
    #ebv = excesso de cor
    ebv[np.isnan(ebv)] = 0
    x = 1 / wave.to(u.micrometer).value
    wtran = (0.63 * u.micrometer)
    kappa = np.where(wave > wtran, 2.659 * (-1.857 + 1.040 * x),
                     2.659 * (-2.156 + 1.509 * x - 0.198 * x * x
                              + 0.011 * (x * x * x)))
    return halpha * np.power(10, -0.4 * ebv * kappa)
    dust = halpha * np.power(10, -0.4 * ebv * kappa)

    vmax = np.nanpercentile(dust, 95)
    vmin = np.nanpercentile(dust, 80)
    plt.imshow(dust, origin="lower", vmax=vmax, vmin=vmin)
    plt.colorbar()
    plt.show()
 
def test_halpha():
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
        halpha_nii_dust = calc_halpha_nii_3F(fnu)
        magAB = -2.5 * np.log10(fnu) - 48.6
        g_i = magAB[3] - magAB[2]
        halpha_nii = dust_correction(halpha_nii_dust, g_i)
        halpha = nii_correction(halpha_nii, g_i)
        vmax = np.percentile(halpha_nii_dust, 99)
        vmin = np.percentile(halpha_nii_dust, 10)
        # Making plot
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(6.36, 2.),
            gridspec_kw = {'wspace':0.05, 'hspace':0.05, 'left':0.02,
                           'right':0.88, "bottom":0.02, "top":0.98})
        axes[0].imshow(halpha_nii_dust, vmax=vmax, vmin=vmin, origin="lower")
        axes[0].legend(title=r"H$\alpha$+[NII]+poeira")
        axes[1].imshow(halpha_nii, vmax=vmax, vmin=vmin, origin="lower",
                       label="3F+")
        axes[1].legend(title=r"H$\alpha$+[NII]")
        im2 = axes[2].imshow(halpha, vmin=vmin, vmax=vmax, origin="lower")
        axes[2].legend(title=r"H$\alpha$")
        for ax in axes:
            ax.get_xaxis().set_ticklabels([])
            ax.get_yaxis().set_ticklabels([])
        # add space for colour bar
        cbar_ax = fig.add_axes([0.91, 0.06, 0.03, 0.85])
        fig.colorbar(im2, cax=cbar_ax)
        plt.savefig(f"correcao_{galaxy}.png", dpi=250)
        plt.show()
        break

if __name__ == "__main__":
    test_halpha()