import os
import numpy as np
import re
from matplotlib import pyplot as plt
from daomop import storage


def std(hpx):
    cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
    table = cat.table

    condition = (table['MATCHES'] > 2)
    hpxids = np.unique(table['HPXID'][condition][:1000])

    mags = [table[table['HPXID'] == hpxid]['MAG_AUTO'] for hpxid in hpxids]
    stds = [(np.std(x)) for x in mags]

    plt.plot(table['MAGERR_AUTO'][:len(hpxids)], stds, 'x', alpha=0.7)
    plt.ylabel("std")
    plt.xlabel("mag error")
    plt.show()

    mean_mag = []
    for i in mags:
        mean_mag.append(np.mean(i))

    plt.subplot(121)
    plt.plot(mean_mag, stds, '.', alpha=0.7)
    plt.ylim(ymin=-0.05, ymax=3/2)
    plt.ylabel("magnitude standard deviation")
    plt.xlabel("average magnitude")

    plt.subplot(122)
    plt.plot(mean_mag, table['MAGERR_AUTO'][:len(hpxids)], '.', alpha=0.7)
    plt.ylabel("magnitude error")
    plt.xlabel("average magnitude")
    plt.ylim(ymin=-0.05, ymax=3/2)
    plt.show()

    plt.plot()
    plt.hist([table['MAGERR_AUTO'][index] / i for index, i in enumerate(stds) if i != 0],
             50,
             range=(-0.2, 5))

    plt.xlabel("Magnitude Error / Magnitude St. Deviation")
    plt.ylabel("Frequency")
    plt.title("Histogram of MAGERR_AUTO / Magnitude Standard Deviation")
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
