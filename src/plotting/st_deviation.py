"""
Generate .FITS data files containing magnitude standard deviations from master catalog tables.
Called from plot.py in /src/ directory.

Usage:
~/octarine/src$ python plot.py --std
"""

import os
import numpy as np
import re

from astropy.io import fits
from daomop import storage

MAG_STD_DATA_DIRECTORY = 'plotting/mag_std_data/'


def std(hpx):
    """
    Generate a data file containing the standard deviations of magnitudes in a given catalog with MATCHES > 2

    :param hpx: catalog HEALPIX
    """
    mag_filename = MAG_STD_DATA_DIRECTORY + str(hpx) + '_magnitudes.fits'
    std_filename = MAG_STD_DATA_DIRECTORY + str(hpx) + '_stds.fits'

    if not os.path.exists(std_filename):
        cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
        table = cat.table

        condition = (table['MATCHES'] > 2)
        hpxids = np.unique(table['HPXID'][condition][:100])
        print len(hpxids)
        mags = np.array([table[table['HPXID'] == hpxid]['MAG_AUTO'] for hpxid in hpxids])
        stds = [(np.std(x)) for x in mags]

        # if os.path.exists(mag_filename):
        #     os.remove(mag_filename)
        # mag_image = fits.PrimaryHDU(data=mags)
        # mag_image.writeto(mag_filename)

        std_image = fits.PrimaryHDU(data=stds)
        std_image.writeto(std_filename)


def main():
    if not os.path.exists(MAG_STD_DATA_DIRECTORY):
        os.mkdir(MAG_STD_DATA_DIRECTORY)

    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))

        if hpx not in reg and hpx > 0:  # alter 'hpx > 0' to use specific files
            reg.append(hpx)
            std(hpx)
