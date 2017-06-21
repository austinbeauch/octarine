"""Mark the stationary sources in a given source catalog by matching with other source catalogs"""
import sys
import storage
import util
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import vstack
import numpy
import argparse
import logging
import traceback
from cadcutils.exceptions import NotFoundException

task = "stationary"
dependency = None


def run(pixel, expnum, ccd, prefix, version, dry_run, force, start_date, end_date):
    """
    Retrieve the catalog from VOSspace, find the matching dataset_name/ccd combos and match against those.

    :param pixel: Which HPX Pixel should we build a catalog for.
    :param ccd: chip to retrieve for matching
    :param expnum: exposure number to retrieve for match
    :param start_date: MJD at start of period to match catalogs in
    :param end_date: MJD of end of period to match catalogs in
    :param force:
    :param dry_run:
    :param version:
    :param prefix:
    """
    message = storage.SUCCESS

    if storage.get_status(task, prefix, expnum, version=version, ccd=ccd) and not force:
        logging.info("{} completed successfully for {} {} {} {}".format(task, prefix, expnum, version, ccd))
        return

    with storage.LoggingManager(task, str(expnum), expnum, ccd, version, dry_run):
        try:
            if dependency is not None and not storage.get_status(dependency, prefix, 
                                                                 expnum, "p", ccd=ccd):
                raise IOError("{} not yet run for {}".format(dependency, expnum))

            # get catalog from the vospace storage area
            logging.info("Getting fits image from VOSpace")

            logging.info("Running match on %s %d" % (expnum, ccd))
            catalog = match(pixel, expnum, ccd, start_date, end_date)
            split_to_hpx(pixel, catalog)

            if dry_run:
                return

            # place the results into VOSpace
            logging.info(message)
        except Exception as e:
            logging.debug(traceback.format_exc())
            logging.debug(type(e))
            message = str(e)
            logging.error(message)

        storage.set_status(task, prefix, expnum, version, ccd=ccd, status=message)


def split_to_hpx(pixel, catalog, catalog_dir=None):
    dataset_name = catalog.observation.dataset_name
    pix = pixel
    try:
        healpix_catalog = storage.HPXCatalog(pixel=pix, catalog_dir=catalog_dir)
        healpix_catalog.get()
        healpix_catalog.table = healpix_catalog.table[healpix_catalog.table['dataset_name'] != dataset_name]
        healpix_catalog.table = vstack([healpix_catalog.table, catalog.table[catalog.table['HEALPIX'] == pix]])
    except NotFoundException:
        healpix_catalog = storage.HPXCatalog(pixel=pix)
        healpix_catalog.hdulist = fits.HDUList()
        healpix_catalog.hdulist.append(catalog.hdulist[0])
        healpix_catalog.table = catalog.table[catalog.table['HEALPIX'] == pix]
    healpix_catalog.write()
    healpix_catalog.put()


def match(pixel, expnum, ccd, start_date, end_date):

    observation = storage.Observation(expnum)

    catalog = storage.FitsTable(observation, ccd=ccd, ext='.cat.fits')
    dataset_name = "{}{}{}".format(catalog.observation.dataset_name, catalog.version, catalog.ccd)
    image = storage.FitsImage(catalog.observation, ccd=catalog.ccd, version=catalog.version)
    catalog.table['dataset_name'] = len(catalog.table)*[dataset_name]
    catalog.table['mid_mjdate'] = image.header['MJDATE'] + image.header['EXPTIME']/24./3600.0
    catalog.table['exptime'] = image.header['EXPTIME']

    # First match against the HPX catalogs (if they exist)
    ra_dec = SkyCoord(catalog.table['X_WORLD'],
                      catalog.table['Y_WORLD'],
                      unit=('degree', 'degree'))
    catalog.table['HEALPIX'] = util.skycoord_to_healpix(ra_dec)

    npts = numpy.sum([catalog.table['MAGERR_AUTO'] < 0.002])
    if npts < 10:
        flux_radius_lim = 1.8
    else:
        flux_radius_lim = numpy.median(catalog.table['FLUX_RADIUS'][catalog.table['MAGERR_AUTO'] < 0.002])

    datasec = storage.datasec_to_list(image.header['DATASEC'])
    trim_condition = numpy.all((catalog.table['X_IMAGE'] > datasec[0],
                                catalog.table['X_IMAGE'] < datasec[1],
                                catalog.table['Y_IMAGE'] > datasec[2],
                                catalog.table['Y_IMAGE'] < datasec[3],
                                catalog.table['MAG_PSF'] < 99,
                                catalog.table['FLUX_RADIUS'] > flux_radius_lim), axis=0)

    catalog.table = catalog.table[trim_condition]
    match_list = image.polygon.cone_search(runids=storage.RUNIDS,
                                           minimum_time=2.0/24.0,
                                           mjdate=image.header.get('MJDATE', None),
                                           start_date=start_date,
                                           end_date=end_date)

    # First match against the HPX catalogs (if they exist)
    # reshape the position vectors from the catalogues for use in match_lists
    p1 = numpy.transpose((catalog.table['X_WORLD'],
                          catalog.table['Y_WORLD']))

    # Build the HPXID column by matching against the HPX catalogs that might exit.
    catalog.table['HPXID'] = -1
    healpix = pixel

    hpx_cat = storage.HPXCatalog(pixel=healpix, catalog_dir='catalogs_master')
    hpx_cat_len = 0

    try:
        hpx_cat.get()
        p2 = numpy.transpose((hpx_cat.table['X_WORLD'],
                              hpx_cat.table['Y_WORLD']))
        idx1, idx2 = util.match_lists(p1, p2, tolerance=0.5 / 3600.0)
        catalog.table['HPXID'][idx2.data[~idx2.mask]] = hpx_cat.table['HPXID'][~idx2.mask]
        hpx_cat_len = len(hpx_cat.table)
    except NotFoundException:
        pass

    # for all non-matched sources in this healpix we increment the counter.
    cond = numpy.all((catalog.table['HPXID'] < 0,
                      catalog.table['HEALPIX'] == healpix), axis=0)
    catalog.table['HPXID'][cond] = [hpx_cat_len + numpy.arange(cond.sum()), ]
    catalog.table['MATCHES'] = 0
    catalog.table['OVERLAPS'] = 0
    # Now append these new source (cond) to the end of the master catalog.
    split_to_hpx(pixel, catalog, catalog_dir='catalogs_master')

    for match_set in match_list:
        logging.info("trying to match against catalog {}p{:02d}.cat.fits".format(match_set[0], match_set[1]))
        try:
            match_catalog = storage.FitsTable(storage.Observation(match_set[0]), ccd=match_set[1], ext='.cat.fits')
            match_image = storage.FitsImage(storage.Observation(match_set[0]), ccd=match_set[1])
            datasec = storage.datasec_to_list(match_image.header['DATASEC'])

            npts = numpy.sum([match_catalog.table['MAGERR_AUTO'] < 0.002])
            if npts < 10:
                flux_radius_lim = 1.8
            else:
                flux_radius_lim = numpy.median(
                    match_catalog.table['FLUX_RADIUS'][match_catalog.table['MAGERR_AUTO'] < 0.002])
                
            trim_condition = numpy.all((match_catalog.table['X_IMAGE'] > datasec[0],
                                        match_catalog.table['X_IMAGE'] < datasec[1],
                                        match_catalog.table['Y_IMAGE'] > datasec[2],
                                        match_catalog.table['Y_IMAGE'] < datasec[3],
                                        match_catalog.table['MAG_PSF'] < 99,
                                        match_catalog.table['FLUX_RADIUS'] > flux_radius_lim), axis=0)

            match_catalog.table = match_catalog.table[trim_condition]

            # reshape the position vectors from the catalogues for use in match_lists
            p2 = numpy.transpose((match_catalog.table['X_WORLD'],
                                  match_catalog.table['Y_WORLD']))
            idx1, idx2 = util.match_lists(p1, p2, tolerance=0.5/3600.0)
            catalog.table['MATCHES'][idx2.data[~idx2.mask]] += 1
            catalog.table['OVERLAPS'] += \
                [match_image.polygon.isInside(row['X_WORLD'], row['Y_WORLD']) for row in catalog.table]
        except NotFoundException:
            pass

    return catalog


def main():
    parser = argparse.ArgumentParser(
        description='Create a matches column in a source catalog to determine if a source is a stationary object.')

    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("--catalogs",
                        action="store",
                        default="catalogs",
                        help='dbimages subdirectory where catalogs will be stored.')
    parser.add_argument("healpix",
                        type=int,
                        help="healpix to process")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="DRY RUN, don't copy results to VOSpace, implies --force")
    parser.add_argument("--verbose", "-v",
                        action="store_true")
    parser.add_argument("--force", default=False,
                        action="store_true")
    parser.add_argument("--debug", "-d",
                        action="store_true")
    parser.add_argument("start-date", help="MJD of start of darkrun to build HPX catalog for.")
    parser.add_argument("end-date", help="MJD of end of darkrun to build HPX catalog for.")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()

    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages
    storage.CATALOG = args.catalogs
    prefix = ''
    version = 'p'

    exit_code = 0
    overlaps = storage.MyPolygon.from_healpix(args.healpix).cone_search(runids=storage.RUNIDS)
    for overlap in overlaps:
        expnum = overlap[0]
        ccd = overlap[1]
        run(args.healpix, expnum, ccd, prefix, version, args.dry_run, args.force, args.start_date, args.end_datae)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
