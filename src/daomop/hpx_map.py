import healpy as hp
import numpy as np
from astropy.table import Table
from matplotlib import pyplot
import os


def build_map(filename):
    """
    plot the RA/DEC values of all objects in a given file, colour points by value of expnum.
    """

    t = Table.read(filename)

    expnums = np.array([x.split('p')[0] for x in t['dataset_name']])
    colours = ['r', 'g', 'b', 'y']
    count = 0
    for expnum in np.unique(expnums):
        colour = colours[count % len(colours)]
        cond = np.all((expnums == expnum, t['MATCHES'] <= n_matches, t['MAGERR_AUTO'] < 0.25), axis=0)
        pyplot.plot(t['X_WORLD'][cond], t['Y_WORLD'][cond], ',{}'.format(colour), ms=1, alpha=.25)
        count += 1

    nside = 32
    corners = hp.boundaries(nside, 3065)

    pyplot.xlabel("RA")
    pyplot.ylabel("DEC")
    pyplot.title(filename)

    # Now convert to phi,theta representation:
    phi_theta = hp.vec2ang(np.transpose(corners), lonlat=True)
    x = list(phi_theta[0]).append(phi_theta[0][0])
    y = list(phi_theta[1]).append(phi_theta[1][0])
    pyplot.plot(x, y, '-k')
    map_file_name = os.path.splitext(os.path.basename(filename))[0]+".png"

    pyplot.savefig(map_file_name)

if __name__ == '__main__':
    build_map(sys.argv[1])
