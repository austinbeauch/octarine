import os
import numpy as np
import re
from astropy.io import fits
from astropy import wcs

from daomop import storage

PIXEL_WIDTH = 90.0/3600.0


def find_nearest(array, value):
    return array[(np.abs(array-value)).argmin()]


def arduino_map(x, in_min, in_max, out_min, out_max):
    return (x - in_min) * (out_max - out_min) // (in_max - in_min) + out_min


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

    qrunids = np.unique(table['QRUNID'])

    mag_lims = {}
    print "mag lims ", len(table_entries)
    for table_entry in table_entries:
        # print table_entry['frame'][0], len(table_entry['dataset_name'])
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
    overlap_data = np.zeros((len(dec_idx), len(ra_idx)))
    stellar_density_data = np.zeros((len(dec_idx), len(ra_idx)))

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
    for qrun in qrunids:
        mag_image_filename = 'fits_data/' + str(hpx) + '_' + qrun + '_mag_data.fits'
        overlap_image_filename = 'fits_data/' + str(hpx) + '_' + qrun + '_overlap_image.fits'
        density_image_filename = 'fits_data/' + str(hpx) + '_' + qrun + '_density_image.fits'
        print qrun, mag_image_filename, overlap_image_filename
        if not os.path.exists(mag_image_filename):

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
                            for frame in dataset_names:
                                try:
                                    l.append(float(mag_lims[frame]))  # getting third faintest magnitude
                                except KeyError:
                                    print "Key {} not in {}.".format(frame, mag_lims)
                                    continue

                            third_faintest = sorted(l)[-3]
                            print 'index: ({} , {})'.format(yy, xx)
                            mag_data[yy, xx] = third_faintest
                            overlap_data[yy, xx] = len(dataset_names)

            for row in mag_data:
                if row.any() != 0:
                    mag_image = fits.PrimaryHDU(data=mag_data, header=header)
                    mag_image.writeto(mag_image_filename)

                    overlap_image = fits.PrimaryHDU(data=overlap_data, header=header)
                    overlap_image.writeto(overlap_image_filename)

                    density_image = fits.PrimaryHDU(data=stellar_density_data, header=header)
                    density_image.writeto(density_image_filename)
                    break


def main():
    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))

        if hpx not in reg and hpx == 2811:
            reg.append(hpx)
            try:
                fits_factory(hpx)

            except Exception as ex:
                print ex
                raise ex

    # condition = (table['OVERLAPS'] > 2)
    # hpxids = np.unique(table['HPXID'][condition][:2000])
    # x = table['MAG_AUTO']
    # y = table['MAGERR_AUTO']
    # idx = np.random.choice(np.arange(len(x)), 1000)
    # x, y = x[idx], y[idx]
    # xx, yy = np.meshgrid(x, y)
    # plt.plot(xx, yy, '.')
    # plt.show()
    # magni2ds = {}
    # for data_set in table_entries:
    #     nearest_magerr = find_nearest(data_set['MAGERR_AUTO'], 0.2)
    #     condition = (table['MAGERR_AUTO'] == nearest_magerr)
    #     nearest_magnitude = float(table['MAG_AUTO'][condition])
    #     magni2ds[np.unique(data_set['dataset_name']).tostring().rstrip('\x00')] = [nearest_magerr, nearest_magnitude]
    # cat = storage.HPXCatalog(1161, catalog_dir='catalogs/master', dest_directory='master')
    # table = cat.table
    # mags = [table[table['HPXID'] == hpxid]['MAG_AUTO'] for hpxid in hpxids]
    # stds = [(numpy.std(x)) for x in mags]
    # plt.plot(table['MAGERR_AUTO'][:2000], stds, 'x', alpha=0.7)
    # plt.ylabel("std")
    # plt.xlabel("mag error")
    # plt.show()
    # mean_mag = []
    # for i in mags:
    #     mean_mag.append(numpy.mean(i))
    # plt.subplot(121)
    # plt.plot(mean_mag, stds, '.', alpha=0.7)
    # plt.ylim(ymin=-0.05, ymax=3/2)
    # # plt.xlim(xmin=0)
    # plt.ylabel("magnitude standard deviation")
    # plt.xlabel("average magnitude")
    # plt.subplot(122)
    # plt.plot(mean_mag, table['MAGERR_AUTO'][:2000], '.', alpha=0.7)
    # plt.ylabel("magnitude error")
    # plt.xlabel("average magnitude")
    # plt.ylim(ymin=-0.05, ymax=3/2)
    # plt.show()


if __name__ == "__main__":
    main()
