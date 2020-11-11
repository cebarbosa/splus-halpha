""" Obtaining calibrated datacubes for sample of galaxies. """
import os

from astropy.table import Table
import splusdata

from splus_cubes import make_cube

import context

if __name__ == "__main__":
    tablename = os.path.join(context.tables_dir, "kennicutt2008.fits")
    sample = Table.read(tablename)
    outdir = os.path.join(context.data_dir, "scubes-11HUGS")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)
    for galaxy in sample:
        for band in context.bands:
            splusdata.get_fits(galaxy["Name"], galaxy["RAJ2000"], galaxy[
                "DEJ2000"], 100, band, filename="test.fits")
            input(404)


