import os
import numpy as np
import re
from astropy.io import fits
from astropy import wcs
from matplotlib import pyplot as plt
from daomop import storage


def find_nearest(array, value):
    return array[(np.abs(array-value)).argmin()]


def arduino_map(x, in_min, in_max, out_min, out_max):
    return (x - in_min) * (out_max - out_min) // (in_max - in_min) + out_min


def std(hpx):
    cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
    table = cat.table

    table['frame'] = [dataset_name.split('p')[0] for dataset_name in table['dataset_name']]
    dataset_names = np.unique(table['frame'])
    table_entries = [table[table['frame'] == name] for name in dataset_names]
    condition = (table['OVERLAPS'] > 2)
    hpxids = np.unique(table['HPXID'][condition][:2000])

    # magni2ds = {}
    # for data_set in table_entries:
    #     nearest_magerr = find_nearest(data_set['MAGERR_AUTO'], 0.2)
    #     condition = (table['MAGERR_AUTO'] == nearest_magerr)
    #     nearest_magnitude = float(table['MAG_AUTO'][condition])
    #     magni2ds[np.unique(data_set['dataset_name']).tostring().rstrip('\x00')] = [nearest_magerr, nearest_magnitude]

    # cat = storage.HPXCatalog(1161, catalog_dir='catalogs/master', dest_directory='master')
    # table = cat.table

    mags = [table[table['HPXID'] == hpxid]['MAG_AUTO'] for hpxid in hpxids]
    stds = [(np.std(x)) for x in mags]
    plt.plot(table['MAGERR_AUTO'][:2000], stds, 'x', alpha=0.7)
    plt.ylabel("std")
    plt.xlabel("mag error")
    plt.show()
    mean_mag = []
    for i in mags:
        mean_mag.append(np.mean(i))
    plt.subplot(121)
    plt.plot(mean_mag, stds, '.', alpha=0.7)
    plt.ylim(ymin=-0.05, ymax=3/2)
    # plt.xlim(xmin=0)
    plt.ylabel("magnitude standard deviation")
    plt.xlabel("average magnitude")
    plt.subplot(122)
    plt.plot(mean_mag, table['MAGERR_AUTO'][:2000], '.', alpha=0.7)
    plt.ylabel("magnitude error")
    plt.xlabel("average magnitude")
    plt.ylim(ymin=-0.05, ymax=3/2)
    plt.show()


def main():
    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))

        if hpx not in reg and hpx == 2811:  # alter 'hpx > 0' to use specific files
            reg.append(hpx)
            std(hpx)


if __name__ == "__main__":
    main()
