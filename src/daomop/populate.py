import argparse
import sys
import logging
import util

import storage
from storage import archive_url
from storage import Observation
from storage import Image
from storage import make_path
from storage import make_link
from storage import cfis_uri
from storage import isfile


def run(dataset_name):
    """Given a dataset_name created the desired dbimages directories
    and links to the raw and processed data files stored at CADC and vospace.

    @param dataset_name: the name of the CFHT dataset to make a link to.
    """

    observation = Observation(dataset_name)

    version = 'o'
    ext = '.fits.fz'
    artifact = Image(observation, version=version, ext=ext)
    source = archive_url(dataset_name, version)
    logging.debug("Making link between {} and {}".format(source, artifact.uri))
    make_path(artifact.uri)
    make_link(source, artifact.uri)

    logging.debug("Linking to RAW header")
    source = archive_url(dataset_name, version=version, fhead='true')
    make_link(source, artifact.header_artifact.uri)

    logging.debug("Making link between PROC'd image and dbimages")
    version = 'p'
    ext = '.fits.fz'
    artifact = Image(observation, version=version, ext=ext)
    # Source is either CFIS processing or archive URL
    source = cfis_uri(dataset_name)
    if not isfile(source):
        source = archive_url(dataset_name, version)
    make_link(source, artifact.uri)

    logging.debug("Link up the header")
    # Can be the CFIS produced header or the one in the archive
    source = archive_url(dataset_name, version, ext='.head', archive='CFHTSG')
    if not isfile(source):
        source = archive_url(dataset_name, version, ext='', fhead='true')
    make_link(source, artifact.header_artifact.uri)

    return True


def main():
    parser = argparse.ArgumentParser(
        description='Create the vospace entries required for pipeline processing')

    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("expnum",
                        type=int,
                        nargs='+',
                        help="expnum(s) to create directories for")
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

    exit_code = 0
    for expnum in args.expnum:
        run(expnum)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())