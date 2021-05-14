
import os

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import context

def get_names(wdir, bands=["R", "F660", "I"]):
    """ Obtaining names of images  for one directory. """
    filenames = [x for x in os.listdir(wdir) if x.endswith("_swp.fits")]
    names = []
    for band in bands:
        names.append([x for x in filenames if x.split("_")[2]==band][0])
    return names

if __name__ == "__main__":
    data_dir = os.path.join(context.data_dir, "11HUGS/cutouts")
    galaxies = os.listdir(data_dir)
    coef_F660 = 125.3
    coef_R = 1419.0
    for galaxy in galaxies:
        wdir = os.path.join(data_dir, galaxy)
        os.chdir(wdir)
        imgnames = get_names(wdir)
        # Loading data
        data = np.array([fits.getdata(name, ext=1) for name in imgnames])
        # Loading zero points
        m0 = np.array([fits.getval(name, "MAGZP", ext=1) for name in imgnames])
        fnu =  data * np.power(10, -0.4 * (m0[:, None, None] + 48.6))
        R = fnu[0]
        F660 = fnu[1]
        fluxTwoBands = coef_F660 * (F660 - R) / (1-(coef_F660/coef_R))
        vmax = np.nanpercentile(fluxTwoBands, 95)
        vmin = np.nanpercentile(fluxTwoBands, 80)
        plt.imshow(fluxTwoBands, origin="lower", vmax=vmax, vmin=vmin)
        plt.colorbar()
        plt.show()


