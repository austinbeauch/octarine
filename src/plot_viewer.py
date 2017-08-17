import os
import sys
import re
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs


def main():
    directory = sorted(os.listdir('./fits_data/'))
    directory_names = []
    for name in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('(?P<number>\d{3,5})_', name)
        if x.group('number') not in directory_names:
            directory_names.append(x.group('number'))

    while True:
        print "Enter one of: ",
        for i in directory_names:
            print i,
        print
        print "Type 'exit' to quit"
        hpx = raw_input("Enter HPX to display: ")

        if hpx == 'exit':
            sys.exit(0)
        else:
            qrun = None
            for filename in directory:
                x = re.match('(?P<number>\d{3,5})_17AQ(?P<qrun>\d{2})', filename)
                if qrun is None:
                    qrun = '17AQ' + x.group('qrun')

                if hpx and qrun and 'density' in filename:
                    hdu = fits.open('./fits_data/' + filename)[0]
                    w = wcs.WCS(hdu.header)
                    plt.subplot(313, projection=w)
                    plt.imshow(hdu.data, origin='lower', cmap='binary')
                    plt.title(x.group('number') + ' ' + qrun + ' Stellar Densities')
                    plt.xlabel('RA')
                    plt.ylabel('Dec')

                elif hpx and qrun and '_mag_' in filename:
                    hdu = fits.open('./fits_data/' + filename)[0]
                    w = wcs.WCS(hdu.header)
                    plt.subplot(311, projection=w)
                    plt.imshow(hdu.data, origin='lower', cmap='binary')
                    plt.title(x.group('number') + ' ' + qrun + ' Magnitudes')
                    plt.xlabel('RA')
                    plt.ylabel('Dec')

                elif hpx and qrun and 'overlap' in filename:
                    hdu = fits.open('./fits_data/' + filename)[0]
                    w = wcs.WCS(hdu.header)
                    plt.subplot(312, projection=w)
                    plt.imshow(hdu.data, origin='lower', cmap='binary')
                    plt.title(x.group('number') + ' ' + qrun + ' Overlaps')
                    plt.xlabel('RA')
                    plt.ylabel('Dec')
                    plt.show()
                    qrun = None


if __name__ == '__main__':
    main()
