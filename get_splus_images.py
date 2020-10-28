""" Retrieves images from S-PLUS data server. """
import os

import splusdata

import context

if __name__ == "__main__":
    galaxies = ["galaxy0"]
    ras = [313.229166666667]
    decs = [-0.7]
    sizes = [3000]
    for i, gal in enumerate(galaxies):
        for band in context.bands:
            img =splusdata.get_fits('test{}_{}'.format(i, band),
                                 ras[i], decs[i], sizes[i], band)
            print(img)
            input()