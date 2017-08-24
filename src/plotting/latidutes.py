import os
import re
import numpy as np
from astropy.io import fits
from astropy.coordinates import HeliocentricTrueEcliptic, SkyCoord
from astropy import wcs, units

from .plot_viewer import FITS_DATA_DIR, HELIO_DATA_DIR


def latitude_factory(hpx, hpx_files):
    """
    Generate expected counts per latitude data files. Opens up magnitude data files from plotting/mag_data/
    and generates the expected amount of objects per latitude degree.

    :param hpx: HEALPIX value
    :param hpx_files: List of all HEALPIX values from present files in /mag_data/
    """
    # hpx_files is sorted, allowing a simple iteration over the list
    for filename in hpx_files:
        if hpx in filename and '_mag_' in filename:
            x = re.match('(?P<number>\d{3,5})_(?P<qrun>\d{2}[A-z]{2}\d{2})', filename)
            qrun = x.group('qrun')
            if len(hpx) == 3:
                hpx = '0' + hpx

            hdu = fits.open(FITS_DATA_DIR + filename)[0]
            w = wcs.WCS(hdu.header)

            helio_data = np.zeros(90)
            helio_filename = HELIO_DATA_DIR + hpx + '_' + qrun + '_latitude_counts.fits'
            print helio_filename, os.path.exists(helio_filename)

            if os.path.exists(helio_filename):
                continue

            for i in range(hdu.data.shape[0]):
                for j in range(hdu.data.shape[1]):
                    # print i, j
                    ra, dec = w.all_pix2world(i, j, 0)
                    c = SkyCoord(ra, dec, 100 * units.au, unit=('degree', 'degree'))
                    coord = c.transform_to(HeliocentricTrueEcliptic)

                    # if hdu.data[i, j] != 0:
                    count = 0.000625 * 10 ** (0.9 * (hdu.data[i, j] - 23.4))
                    helio_data[int(coord.lat.deg)] += count

            # print helio_data
            mag_image = fits.PrimaryHDU(data=helio_data)
            mag_image.writeto(helio_filename)


def main():
    if not os.path.exists(HELIO_DATA_DIR):
        os.mkdir(HELIO_DATA_DIR)

    directory = sorted(os.listdir(FITS_DATA_DIR))
    hpx_values = []
    for name in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('(?P<number>\d{3,5})_', name)
        if int(x.group('number')) not in hpx_values:
            hpx_values.append(int(x.group('number')))

    for i in sorted(hpx_values):
        latitude_factory(str(i), sorted(directory))
