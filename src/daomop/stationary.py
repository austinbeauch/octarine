"""Mark the stationary sources in a given source catalog by matching with other source catalogs"""
import sys
import errno
import storage
import util
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import vstack
import numpy
import argparse
import logging

task="stationary"
dependency = None

def add_column(hdu, column_name, column_format, data):
    """

    @param hdu:  The FITS Binary table to add the column to
    @param column_name: The name of the Column to add
    @param column_format: The data format for the column
    @param data: The data for the column
    @return: fits.BinTableHDU
    """
    orig_cols = hdu.data.columns
    new_cols = fits.ColDefs([fits.Column(name=column_name, format=column_format, array=data),])
    new_hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols, name=hdu.name)
    return new_hdu


def run(expnum, ccd, prefix, version, dry_run, force):
    """
    Retrieve the catalog from VOSspace, find the matching dataset_name/ccd combos and match against those.

    @param expnum: exposure number to retrieve for match
    @param ccd: chip to retrieve for matching
    @return:
    """
    message = storage.SUCCESS

    if storage.get_status(task, prefix, expnum, version=version, ccd=ccd) and not force:
        logging.info("{} completed successfully for {} {} {} {}".format(task, prefix, expnum, version, ccd))
        return

    with storage.LoggingManager(task, prefix, expnum, ccd, version, dry_run):
        try:
            if dependency is not None and not storage.get_status(dependency, prefix, 
                                                                 expnum, "p", ccd=ccd):
                raise IOError("{} not yet run for {}".format(dependency, expnum))

            # get catalog from the vospace storage area
            logging.info("Getting fits image from VOSpace")

            logging.info("Running match on %s %d" % (expnum, ccd))
            catalog = match(expnum, ccd)
            split_to_hpx(catalog)

            if dry_run:
                return

            # place the results into VOSpace
            logging.info(message)
        except Exception as e:
            message = str(e)
            logging.error(message)

        storage.set_status(task, prefix, expnum, version, ccd=ccd, status=message)

NSIDE = 32
def split_to_hpx(catalog, nside=NSIDE):

    number_of_pix = 12 * nside**2
    field_size = len(str(number_of_pix))
    table = catalog.table
    dataset_name = "{}{}{}".format(catalog.observation.dataset_name, catalog.version, catalog.ccd)
    image = storage.Image(catalog.observation, ccd=catalog.ccd, version=catalog.version)
    catalog.table['dataset_name'] = len(catalog.table)*[dataset_name]
    catalog.table['mid_mjdate'] = image.header['MJDATE'] + image.header['EXPTIME']/24./3600.0

    ra_dec = SkyCoord(catalog.table['X_WORLD'], 
                      catalog.table['Y_WORLD'],
                      unit=('degree', 'degree'))
    healpix = util.skycoord_to_healpix(ra_dec, nside=nside)
    for pix in numpy.unique(healpix):
        hpx_dataset_name = ("{"+":0{:d}".format(field_size)+"}").format(pix)
        skycoord = util.healpix_to_skycoord(pix, nside=nside)
        hpx_dataset_name = "HPX_{}_RA_{:4.1f}_DEC_{:+4.1f}".format(hpx_dataset_name, skycoord.ra.degree, skycoord.dec.degree)
        observation = storage.Observation(hpx_dataset_name, dbimages="vos:cfis/solar_system/catalogs/")
        healpix_catalog = storage.FitsTable(observation, 
                                            ccd=None, 
                                            subdir="", 
                                            version='_cat', 
                                            prefix=None, 
                                            ext='.fits')
        try:
            healpix_catalog.get()
            healpix_catalog.table = healpix_catalog.table[healpix_catalog.table['dataset_name']!=dataset_name]
            healpix_catalog.table = vstack([healpix_catalog.table, catalog.table[healpix==pix]])
        except Exception as ex:
            if ex.errno == errno.ENOENT:
                healpix_catalog.hdulist=fits.HDUList()
                healpix_catalog.hdulist.append(catalog.hdulist[0])
                healpix_catalog.table = catalog.table[healpix==pix]
            else:
                raise(ex)
        healpix_catalog.write()
        healpix_catalog.put()

def match(expnum, ccd):

    observation = storage.Observation(expnum)
    image = storage.Image(observation, ccd=ccd)
    match_list = image.overlaps()
    catalog = storage.FitsTable(observation, ccd=ccd, ext='.cat.fits')
    # reshape the position vectors from the catalogues for use in match_lists
    p1 = numpy.transpose((catalog.table['X_WORLD'],
                          catalog.table['Y_WORLD']))
    matches = numpy.zeros(len(catalog.table['X_WORLD']))
    overlaps = numpy.zeros(len(catalog.table['X_WORLD']))
    for match_set in match_list:
        logging.info("trying to match against catalog {}p{:02d}.cat.fits".format(match_set[0], match_set[1]))
        try:
            match_catalog = storage.FitsTable(storage.Observation(match_set[0]), ccd=match_set[1], ext='.cat.fits')
            match_image = storage.Image(storage.Observation(match_set[0]), ccd=match_set[1])
        except Exception as ex:
            logging.error(ex)
            logging.error(type(ex))
            continue
        except OSError as ioe:
            logging.debug(str(ioe))
            continue
            
        overlaps += [match_image.polygon.isInside(row['X_WORLD'], row['Y_WORLD']) for row in catalog.table]
        # reshape the position vectors from the catalogues for use in match_lists
        p2 = numpy.transpose((match_catalog.table['X_WORLD'],
                              match_catalog.table['Y_WORLD']))
        idx1, idx2 = util.match_lists(p1, p2, tolerance=1.0/3600.0)
        matches[idx2.data[~idx2.mask]] += 1

    catalog.table['MATCHES'] = matches
    catalog.table['OVERLAPS'] = overlaps

    return catalog


def main():
    parser = argparse.ArgumentParser(
        description='Create a matches column in a source catalog to determine if a source is a stationary object.')

    parser.add_argument('--ccd', '-c',
                        action='store',
                        type=int,
                        dest='ccd',
                        default=None,
                        help='which ccd to process, default is all')
    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("dataset_name",
                        type=int,
                        nargs='+',
                        help="dataset_name(s) to process")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="DRY RUN, don't copy results to VOSpace, implies --force")
    parser.add_argument("--verbose", "-v",
                        action="store_true")
    parser.add_argument("--force", default=False,
                        action="store_true")
    parser.add_argument("--debug", "-d",
                        action="store_true")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()

    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages
    prefix = ''
    version = 'p'

    exit_code = 0
    for expnum in args.dataset_name:
        if args.ccd is None:
           if int(expnum) < 1785619:
               # Last exposures with 36 CCD Megaprime
               ccdlist = range(0,36)
           else:
               # First exposrues with 40 CCD Megaprime
               ccdlist = range(0, 40)
        else:
           ccdlist = [args.ccd]
        for ccd in ccdlist:
            run(expnum, ccd, prefix, version, args.dry_run, args.force)
    return exit_code

if __name__ == '__main__':
    sys.exit(main())
