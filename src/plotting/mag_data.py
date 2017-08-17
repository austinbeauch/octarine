import os
import numpy as np
import re
from astropy.io import fits
from astropy import wcs

from daomop import storage

PIXEL_WIDTH = 90.0/3600.0
FITS_DIRECTORY = 'fits_data/'


class IDX(object):
    def __init__(self, range_min, range_max, width):
        self.min = range_min
        self.max = range_max
        self.width = width
        self.nbin = int((range_max - range_min) / width) + 1
        self.centre = (self.max+self.min) / 2

    def __len__(self):
        return self.nbin


def fits_factory(hpx):
    print hpx

    cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
    table = cat.table

    table['frame'] = [dataset_name.split('p')[0] for dataset_name in table['dataset_name']]
    dataset_names = np.unique(table['frame'])
    table_entries = [table[table['frame'] == name] for name in dataset_names]

    mag_lims = {}

    for table_entry in table_entries:
        previous_mag = None
        for mag in np.arange(table_entry['MAG_AUTO'].min(), table_entry['MAG_AUTO'].max(), 0.20):
            condition = np.all((table_entry['MAG_AUTO'] > mag, table_entry['MAG_AUTO'] <= mag + 0.2), axis=0)
            magerr = np.median(table_entry['MAGERR_AUTO'][condition])
            if magerr > 0.15:
                if previous_mag is None:
                    previous_mag = mag
                else:
                    frame = table_entry['frame'][0].split('p')[0]
                    mag_lims[frame] = mag
                    break
            else:
                previous_mag = None

    ra_idx = IDX(table['X_WORLD'].min(), table['X_WORLD'].max(), PIXEL_WIDTH)
    dec_idx = IDX(table['Y_WORLD'].min(), table['Y_WORLD'].max(), PIXEL_WIDTH)

    mag_data = np.zeros((len(dec_idx), len(ra_idx)))

    w = wcs.WCS(naxis=2)
    w.wcs.cd = [[PIXEL_WIDTH, 0], [0, PIXEL_WIDTH]]
    crpix1 = mag_data.shape[1]/2
    crpix2 = mag_data.shape[0]/2
    w.wcs.crpix = [crpix1, crpix2]
    crval = [ra_idx.centre, dec_idx.centre]
    w.wcs.crval = crval
    w.wcs.ctype = ['ra--tan', 'dec-tan']
    w.wcs.cunit = ['deg', 'deg']
    header = w.to_header()

    print mag_data.shape[1], mag_data.shape[0]
    for qrun in np.unique(table['QRUNID']):

        # set file names and reset data arrays for each qrunid
        mag_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_mag_data.fits'
        overlap_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_overlap_image.fits'
        density_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_density_image.fits'
        mag_data = np.zeros((len(dec_idx), len(ra_idx)))
        overlap_data = np.zeros((len(dec_idx), len(ra_idx)))
        stellar_density_data = np.zeros((len(dec_idx), len(ra_idx)))
        print qrun, mag_image_filename

        if not os.path.exists(density_image_filename):
            count = 0
            for xx in range(mag_data.shape[1]):
                print count
                count += 1

                for yy in range(mag_data.shape[0]):
                    ra, dec = w.all_pix2world(xx, yy, 0)

                    ra_cond = np.all((table['X_WORLD'] >= ra,
                                      table['X_WORLD'] < ra + PIXEL_WIDTH,
                                      table['OVERLAPS'] > 1,
                                      table['QRUNID'] == qrun),
                                     axis=0)

                    dec_cond = np.all((ra_cond,
                                       table['Y_WORLD'] >= dec,
                                       table['Y_WORLD'] < dec + PIXEL_WIDTH),
                                      axis=0)

                    dataset_names = np.unique(table[dec_cond]['frame'])

                    if len(dataset_names) != 0:
                        density = float(len(table[dec_cond])) / len(dataset_names)
                        stellar_density_data[yy, xx] = density

                        if len(dataset_names) > 2:
                            l = []
                            try:
                                for frame in dataset_names:
                                    l.append(float(mag_lims[frame]))  # getting third faintest magnitude
                                third_faintest = sorted(l)[-3]
                                print 'index: ({} , {})'.format(yy, xx)
                                mag_data[yy, xx] = third_faintest
                                overlap_data[yy, xx] = len(dataset_names)
                            except KeyError:
                                print KeyError
                                continue

            for row in mag_data:
                if row.any() != 0:
                    # only write these files if there's data that has been inserted
                    mag_image = fits.PrimaryHDU(data=mag_data, header=header)
                    mag_image.writeto(mag_image_filename)

                    overlap_image = fits.PrimaryHDU(data=overlap_data, header=header)
                    overlap_image.writeto(overlap_image_filename)
                    break

            # write the stellar density image regardless so it's not regathered again in case of empty mag_data
            density_image = fits.PrimaryHDU(data=stellar_density_data, header=header)
            density_image.writeto(density_image_filename)


def main():
    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))

        if hpx not in reg and hpx > 0:  # alter 'hpx > 0' to use specific files
            reg.append(hpx)
            fits_factory(hpx)


if __name__ == "__main__":
    main()
