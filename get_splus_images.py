""" Retrieves images from S-PLUS data server. """

<<<<<<< HEAD
import os
=======
from astropy.table import Table
>>>>>>> e0e612189c11c14b4a0399fad07aedaff0e53a0e
import splusdata
import context
if __name__ == "__main__":
    tablename = os.path.join(context.tables_dir, "kennicutt2008.fits")
    sample = Table.read(tablename)
    size = 200
    for i, gal in enumerate(sample):
        name = gal["Name"]
        ra = gal["RAJ2000"]
        dec = gal["DEJ2000"]
        # for band in context.bands:
        #     img = splusdata.get_fits('test{}_{}'.format(i, band),
        #                          ra, dec, size, band)
        #     print(img)
