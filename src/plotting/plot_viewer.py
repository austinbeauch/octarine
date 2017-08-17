import os
import sys
import re
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs


def load_images(hpx, hpx_files):
    qrun = None
    for filename in hpx_files:
        x = re.match('(?P<number>\d{3,5})_17AQ(?P<qrun>\d{2})', filename)
        if qrun is None:
            qrun = '17AQ' + x.group('qrun')

        if hpx in filename and qrun in filename and 'density' in filename:
            hdu = fits.open('./fits_data/' + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(313, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='binary')
            plt.title(x.group('number') + ' ' + qrun + ' Stellar Densities')
            plt.xlabel('RA')
            plt.ylabel('Dec')

        elif hpx in filename and qrun in filename and '_mag_' in filename:
            hdu = fits.open('./fits_data/' + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(311, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='viridis', vmin=22)
            plt.title(x.group('number') + ' ' + qrun + ' Magnitudes')
            plt.xlabel('RA')
            plt.ylabel('Dec')

        elif hpx in filename and qrun in filename and 'overlap' in filename:
            hdu = fits.open('./fits_data/' + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(312, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='binary')
            plt.title(x.group('number') + ' ' + qrun + ' Overlaps')
            plt.xlabel('RA')
            plt.ylabel('Dec')
            plt.show()  # overlap file comes last in the directory, wait to show plots until it's reached
            qrun = None
        else:
            qrun = None


def main():
    directory = sorted(os.listdir('./fits_data/'))
    directory_names = []
    for name in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('(?P<number>\d{3,5})_', name)
        if int(x.group('number')) not in directory_names:
            directory_names.append(int(x.group('number')))

    while True:
        print "Enter one of: ",
        for i in sorted(directory_names):
            print i,
        print
        print "Hit <enter> key to view all, or type 'exit' to quit."
        hpx = raw_input("Enter HPX to display: ")

        hpx_files = []
        for filename in directory:
            if hpx in filename:
                hpx_files.append(filename)

        if hpx == 'exit':
            sys.exit(0)

        if hpx == '':
            for i in sorted(directory_names):
                load_images(str(i), hpx_files)

        else:
            load_images(hpx, hpx_files)


if __name__ == '__main__':
    main()
