"""OSSOS VOSpace storage convenience package."""
import logging
import errno
import os
import urllib
import re
import tempfile
import Polygon
import numpy
import requests
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table
from astropy.io import fits, ascii
from astropy.time import Time

import util
import vospace
from wcs import WCS

# Try and turn off warnings, only works for some releases of requests.
try:
    requests.packages.urllib3.disable_warnings()
    requests_log = logging.getLogger("requests.packages.urllib3")
    requests_log.setLevel(logging.ERROR)
except:
    pass

VOS_PROTOCOL = 'vos:'
MAX_RETRY = 10
MAXCOUNT = 30000
_TARGET = "TARGET"
SSOIS_SERVER = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/fixedssos.pl"
DATA_WEB_SERVICE = 'https://www.canfar.phys.uvic.ca/data/pub/'
VOSPACE_WEB_SERVICE = 'https://www.canfar.phys.uvic.ca/vospace/nodes/'
TAP_WEB_SERVICE = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap/sync'
TAG_URI_BASE = 'ivo://canfar.uvic.ca/daomop'
OBJECT_COUNT = "object_count"
ZEROPOINT_KEYWORD = "PHOTZP"
SUCCESS = 'success'
APCOR_EXT = "apcor"
ZEROPOINT_USED_EXT = "zeropoint.used"
PSF_EXT = "psf.fits"
HEADER_EXT = ".head"
IMAGE_EXT = '.fits.fz'
TEXT_EXT = ".txt"
PROCESSED_VERSION = 'p'
RAW_VERSION = 'o'
RUNIDS = ['%P30', '%P31']
DBIMAGES = 'vos:cfis/solar_system/dbimages'
CATALOG = 'catalogs'
PITCAIRN = 'vos:cfis/pitcairn'
FLATS_VOSPACE = 'vos:sgwyn/flats'
ARCHIVE = 'CFHT'
DEFAULT_FORMAT = 'fits'
NSIDE = 32


class MyRequests(object):

    def __init__(self):

        self.requests = requests

    def get(self, *args, **kwargs):
        resp = self.requests.get(*args, **kwargs)
        resp.raise_for_status()
        return resp


requests = MyRequests()


class Task(object):
    """
    A task within the OSSOS pipeline work-flow.
    """

    def __init__(self, executable, dependency=None):
        self.executable = executable
        self.name = os.path.splitext(self.executable)[0]
        self._target = None
        self._status = None
        self._dependency = None
        self.dependency = dependency

    def __str__(self):
        return self.name

    @property
    def tag(self):
        """
        Get the string representation of the tag used to annotate the status in VOSpace.
        @return: str
        """
        return "{}{}_{}{:02d}".format(self.target.prefix,
                                      self,
                                      self.target.version,
                                      self.target.ccd)

    @property
    def target(self):
        """

        @return: The target that this task is set to run on.
        @rtype: Target
        """
        return self._target

    @target.setter
    def target(self, target):
        assert isinstance(Artifact, target)
        self._target = target

    @property
    def dependency(self):
        """
        @rtype: Task
        """
        raise NotImplementedError()

    @dependency.setter
    def dependency(self, dependency):
        if dependency is None:
            self._dependency = dependency
        else:
            assert isinstance(Task, dependency)

    @property
    def status(self):
        """

        @return: The status of running this task on the given target.
        @rtype: str
        """
        return get_tag(self.target.expnum, self.tag)

    @status.setter
    def status(self, status):
        status += Time.now().iso
        set_tag(self.target.expnum, self.tag, status)

    @property
    def finished(self):
        """
        @rtype: bool
        """
        return self.status.startswith(SUCCESS)

    @property
    def ready(self):
        if self.dependency is None:
            return True
        else:
            return self.dependency.finished


class MyPolygon(Polygon.Polygon):

    def __init__(self, *args, **kwargs):
        if 'footprint' in kwargs:
            self.footprint = kwargs['footprint']
            del(kwargs['footprint'])

        super(MyPolygon, self).__init__(*args, **kwargs)


    @classmethod
    def from_footprint(cls, footprint):
        """
        Build a Polygon object using a wcs.calc_footprint output as input.
        """
        return cls(numpy.concatenate((footprint, numpy.array([footprint[0]])), axis=0), footprint=footprint)

    @classmethod
    def from_healpix(cls, healpix, nside=None):
        corners = util.healpix_to_corners(healpix, nside)
        return cls.from_footprint(corners)

    def cone_search(self, runids=None, mjdate=None, minimum_time=None):

        query = (" SELECT Observation.observationID as collectionID "
                 " FROM caom2.Observation AS Observation "
                 " JOIN caom2.Plane AS Plane "
                 " ON Observation.obsID = Plane.obsID "
                 " WHERE Observation.collection = 'CFHT'  "
                 " AND Plane.calibrationLevel = 1 "
                 " AND Plane.energy_bandpassName LIKE 'r.%'  ")

        # Restrict the overlap search to particular runids
        if runids is not None and len(runids) > 0:
            query += " AND ( "
            sep = ""
            for runid in runids:
                query += sep + " Observation.proposal_id LIKE '{}' ".format(runid)
                sep = " OR "
            query += " ) "

        polygon_str = "POLYGON('ICRS GEOCENTER', " + ", ".join([str(x) for x in self.footprint.ravel()]) + ")"

        query += " AND INTERSECTS( {}, Plane.position_bounds ) = 1 ".format(polygon_str)
        if mjdate is not None:
            query += " AND ( Plane.time_bounds_lower < {} ".format(mjdate - minimum_time)
            query += " OR  Plane.time_bounds_upper > {} ) ".format(mjdate + minimum_time)

        table = tap_query(query)
        overlaps = []
        for collectionID in table['collectionID']:
            try:
                headers = Header(Observation(collectionID))
                for header in headers.headers:
                    if header.get('EXTVER', None) is None:
                        continue
                    headers.ccd = int(header['EXTVER'])
                    logging.debug("Checking {} {} ".format(collectionID, headers.ccd))
                    if headers.polygon.overlaps(self):
                        logging.debug("{} {} OVERLAPS ".format(collectionID, headers.ccd))
                        overlaps.append([collectionID, headers.ccd])
            except Exception as ex:
                logging.debug("ERROR processing {}: {}".format(collectionID, ex))
                continue
        logging.debug("Found these overlapping CCDs\n" + str(overlaps))
        return overlaps


class Observation(object):

    def __init__(self, dataset_name, dbimages=None):
        if dbimages is None:
            dbimages = DBIMAGES
        self.dataset_name = dataset_name
        self.dbimages = dbimages

    @property
    def ccd_list(self):
        """
        A list of all the CCDs that this CFHT MegaPrime Exposure should have.
        :return:
        """
        if int(self.dataset_name) < 1785619:
            # Last exposures with 36 CCD MegaPrime
            return range(0, 36)
        return range(0, 40)

    def __str__(self):
        return "{}".format(self.dataset_name)


class Artifact(object):

    def __init__(self, observation, version=PROCESSED_VERSION,
                 ccd=None, subdir=None, ext=IMAGE_EXT, prefix=None):
        """
        :type version: str
        :type observation: Observation
        :type ccd: int
        :type subdir: str
        :type ext: str
        :type prefix: str
        :rtype: Artifact
        """
        self.observation = observation
        self.ccd = ccd
        self._subdir = subdir
        self._ext = ext
        self._version = version
        self._prefix = prefix
        self._hdulist = None

    @property
    def subdir(self):
        """
        :rtype: basestring
        """
        if self._subdir is None:
            self._subdir = self.observation.dataset_name
        return self._subdir

    @property
    def ext(self):
        """
        :rtype: basestring
        """
        if self._ext is None:
            self._ext = ""
        return self._ext

    @property
    def version(self):
        """
        :rtype: basestring
        """

        if self._version is None:
            self._version = ""
        return self._version

    @property
    def prefix(self):
        """
        :rtype: basestring
        """
        if self._prefix is None:
            self._prefix = ""
        return self._prefix

    @property
    def uri(self):
        """
        Build the uri for an OSSOS image stored in the dbimages
        containerNode.

        :rtype : basestring
        dataset_name: CFHT exposure number
        ccd: CCD in the mosaic [0-35]
        version: one of p,s,o etc.
        dbimages: dbimages containerNode.
        """
        uri = "{}/{}".format(self.observation.dbimages, self.subdir)

        if self.ccd is None:
            return "{}/{}{}{}{}".format(uri, self.prefix, self.observation, self.version, self.ext)

        return "{}/ccd{:02d}/{}{}{}{:02d}{}".format(uri, int(self.ccd),
                                                    self.prefix, self.observation, self.version,
                                                    self.ccd, self.ext)

    def get(self):
        """Get the artifact from VOSpace."""
        if not os.access(self.filename, os.F_OK):
            logging.info("Retrieving {} from VOSpace".format(self.uri))
            return copy(self.uri, self.filename)
        return 0

    @property
    def filename(self):
        return os.path.basename(self.uri)

    def put(self):
        """Put the artifact to VOSpace."""
        # first ensure the path exists
        logging.info("Checking that path {} exists".format(self.uri))
        make_path(self.uri)
        logging.info("Copying {} to {}".format(self.filename, self.uri))
        return copy(self.filename, self.uri)

    def delete(self):
        """Delete a file from VOSpace"""
        delete(self.uri)


class FitsArtifact(Artifact):

    @property
    def hdulist(self):
        if self._hdulist is not None:
            return self._hdulist
        if not os.access(self.filename, os.R_OK):
            self.get()
        self._hdulist = fits.open(self.filename)
        return self._hdulist

    @hdulist.setter
    def hdulist(self, hdulist):
        self._hdulist = hdulist


class FitsTable(FitsArtifact):

    def __init__(self, *args, **kwargs):
        super(FitsTable, self).__init__(*args, **kwargs)
        self._table = None

    @property
    def table(self):
        """

        :return: Table object built from this FitsArtifact
        :rtype: Table
        """
        if self._table is None:
            if not os.access(self.filename, os.R_OK):
                self.get()
            self._table = Table.read(self.filename)
        return self._table

    @table.setter
    def table(self, table):
        self._table = table

    def write(self):
        """
        Write the current data in the table to disk.
        """
        if self.table.meta['EXTNAME'] in self.hdulist:
            del(self.hdulist[self.table.meta['EXTNAME']])
        self.hdulist.append(fits.table_to_hdu(self.table))
        self.hdulist.writeto(self.filename, clobber=True)


# TODO: create new function in parallel with mp_ephem package
class FitsImage(FitsArtifact):

    def __init__(self, *args, **kwargs):
        super(FitsImage, self).__init__(*args, **kwargs)
        self._ccd = None
        self.ccd = kwargs.get('ccd', None)
        self._header = None
        self._wcs = None
        self._flat_field_name = None
        self._flat_field = None
        self._footprint = None
        self._polygon = None
    
    @property
    def ccd(self):
        return self._ccd

    @ccd.setter
    def ccd(self, ccd):
        self._wcs = self._polygon = self._header = None
        self._ccd = ccd

    @property
    def filename(self):
        if self.ccd is None:
            return os.path.basename(self.uri)
        return "{}{}{}{:02d}.fits".format(self.prefix,
                                          self.observation,
                                          self.version,
                                          self.ccd)

    @property
    def header(self):
        """
        :return:         The Header of the FITS Extension hold this image.
        :rtype: fits.Header
        """
        if self.ccd is None:
            ext = 0
        else:
            ext = self.ccd + 1
        if self._header is None:
            self._header = Header(self.observation, version=self.version).headers[ext]
        return self._header

    @property
    def wcs(self):
        if self._wcs is None:
            if self.ccd is not None:
                self._wcs = WCS(self.header)
            else:
                self._wcs = [None, ]
                for header in self.header[1:]:
                    self._wcs.append(WCS(header))
        return self._wcs

    @property
    def footprint(self):
        datasec = util.get_pixel_bounds_from_datasec_keyword(self.header['DATASEC'])
        corners = numpy.array([[datasec[0][0], datasec[1][0]],
                            [datasec[0][0], datasec[1][1]],
                            [datasec[0][1], datasec[1][1]],
                            [datasec[0][1], datasec[1][0]]], dtype = numpy.float64)
        self._footprint = self.wcs.wcs_pix2world(corners, 1)
        self._footprint = numpy.concatenate((self._footprint, numpy.array([self._footprint[0]])), axis=0)
        return self._footprint

    @property
    def polygon(self):
        if self._polygon is None:
           self._polygon = MyPolygon.from_footprint(self.footprint)
        return self._polygon

    @property
    def uri(self):
        """
        Build the uri for an OSSOS image stored in the dbimages
        containerNode.

        :rtype : basestring
        dataset_name: CFHT exposure number
        ccd: CCD in the mosaic [0-35]
        version: one of p,s,o etc.
        dbimages: dbimages containerNode.
        """
        uri = "{}/{}".format(self.observation.dbimages, self.subdir)

        if self.ccd is None:
            return "{}/{}{}{}{}".format(uri, self.prefix, self.observation, self.version, self.ext)

        return "{}/{}{}{}{}[{}]".format(uri,
                                        self.prefix, self.observation, self.version, self.ext,
                                        self.ccd+1)

    @property
    def flat_field(self):
        """

        :rtype: FitsImage
        """
        if self._flat_field is None:
            self._flat_field = FitsImage(Observation(self.flat_field_name, dbimages=FLATS_VOSPACE),
                                         subdir="",
                                         ext=".fits",
                                         version="",
                                         ccd=self.ccd)
        return self._flat_field

    @property
    def flat_field_name(self):

        if self._flat_field_name is not None:
            return self._flat_field_name
        if self.ccd is not None:
            self._flat_field_name = self.header.get('FLAT', None)
        else:
            for header in self.header:
                self._flat_field_name = (header.get("FLAT", None))
        if self._flat_field_name is None:
            self._flat_field_name = "weight.fits"
        self._flat_field_name = self._flat_field_name.rstrip(".fits")
        return self._flat_field_name

    @property
    def fwhm(self):
        uri = Artifact(self.observation, ccd=self.ccd, version=self.version, ext='fwhm', prefix=self.prefix).uri
        filename = os.path.basename(uri)
        copy(uri, filename)
        return open(filename).read()


    def cutout(self, cutout, return_file=False):
        """
        Get a sub-section of the image as a cutout.
        :param cutout
        :param return_file
        :return:
        """
        fpt = tempfile.NamedTemporaryFile(suffix='.fits')
        copy(self.uri + cutout, fpt.name)
        fpt.seek(0)
        hdu_list = fits.open(fpt, scale_back=False)
        hdu_list.verify('silentfix+ignore')
        hdu_list[0].header['DATASEC'] = reset_datasec(cutout,
                                                      hdu_list[0].header['DATASEC'],
                                                      hdu_list[0].header['NAXIS1'],
                                                      hdu_list[0].header['NAXIS2'])
        if not hdu_list:
            raise OSError(errno.EFAULT, "Failed to retrieve cutout of image", self.uri)

        self._header = [None]
        for hdu in hdu_list:
            self._header.append(hdu.header)

        if return_file:
            cutout_list = datasec_to_list(cutout)
            if len(cutout_list) == 5:
                cutout_list = cutout_list[1:]
            cutout_filename = os.path.splitext(self.filename)[0] + "_{}_{}_{}_{}.fits".format(cutout_list[0],
                                                                                              cutout_list[1],
                                                                                              cutout_list[2],
                                                                                              cutout_list[3])
            hdu_list.writeto(cutout_filename)
            return cutout_filename

        return hdu_list

    @property
    def zeropoint(self):
        return float(self.header.get(ZEROPOINT_KEYWORD, 30.0))


class HPXCatalog(FitsTable):

    def __init__(self, pixel, version="_cat", ext=".fits", subdir="", nside=None, **kwargs):
        self.pixel = pixel
        if nside is None:
            nside = util.HEALPIX_NSIDE
        self.nside = nside
        dbimages = os.path.join(os.path.dirname(DBIMAGES), CATALOG)
        super(HPXCatalog, self).__init__(Observation(self.dataset_name, dbimages=dbimages),
                                         version=version, ext=ext, subdir=subdir, **kwargs)

    @property
    def skycoord(self):
        return util.healpix_to_skycoord(self.pixel, nside=self.nside)

    @property
    def dataset_name(self):
        number_of_pix = 12 * self.nside ** 2
        field_size = len(str(number_of_pix))
        dataset_name = ("{" + ":0{:d}".format(field_size) + "}").format(self.pixel)
        return "HPX_{}_RA_{:4.1f}_DEC_{:+4.1f}".format(dataset_name,
                                                       self.skycoord.ra.degree,
                                                       self.skycoord.dec.degree)


def set_tags_on_uri(uri, keys, values=None):
    node = vospace.client.get_node(uri)
    if values is None:
        values = []
        for idx in range(len(keys)):
            values.append(None)
    assert (len(values) == len(keys))
    for idx in range(len(keys)):
        key = keys[idx]
        value = values[idx]
        tag = tag_uri(key)
        node.props[tag] = value
    vospace.client.add_props(node)
    return vospace.client.get_node(uri, force=True)


def _set_tags(expnum, keys, values=None):
    uri = os.path.join(DBIMAGES, str(expnum))
    node = vospace.client.get_node(uri, force=True)
    if values is None:
        values = []
        for idx in range(len(keys)):
            values.append(None)
    assert (len(values) == len(keys))
    for idx in range(len(keys)):
        key = keys[idx]
        tag = tag_uri(key)
        value = values[idx]
        node.props[tag] = value
    vospace.client.add_props(node)
    return vospace.client.get_node(uri, force=True)


def set_tags(expnum, props):
    """Assign the key/value pairs in props as tags on on the given dataset_name.

    @param expnum: str
    @param props: dict
    @return: success
    """
    # now set all the props
    return _set_tags(expnum, props.keys(), props.values())


def set_tag(expnum, key, value):
    """Assign a key/value pair tag to the given dataset_name containerNode.

    @param expnum:  str
    @param key: str
    @param value: str
    @return: success
    """

    return set_tags(expnum, {key: value})


def tag_uri(key):
    """Build the uri for a given tag key. 

    @param key: what is the key that we need a stadanrd uri tag for? eg 'mkpsf_00'
    @type key: str
    """
    if TAG_URI_BASE in key:
        return key
    return TAG_URI_BASE + "#" + key.strip()


def get_tag(expnum, key):
    """given a key, return the vospace tag value.

    @param expnum: Number of the CFHT exposure that a tag value is needed for
    @param key: The process tag (such as mkpsf_00) that is being looked up.
    @return: the value of the tag
    @rtype: str
    """

    uri = tag_uri(key)
    force = uri not in get_tags(expnum)
    value = get_tags(expnum, force=force).get(uri, None)
    return value


def get_process_tag(program, ccd, version=PROCESSED_VERSION):
    """
    make a process tag have a suffix indicating which ccd its for.
    @param program: Name of the process that a tag is built for.
    @param ccd: the CCD number that this process ran on.
    @param version: The version of the exposure (s, p, o) that the process ran on.
    @return: The string that represents the processing tag.
    """
    return "%s_%s%s" % (program, str(version), str(ccd).zfill(2))


def get_tags(expnum, force=False):
    """

    @param expnum:
    @param force:
    @return: dict
    @rtype: dict
    """
    uri = os.path.join(DBIMAGES, str(expnum))
    return vospace.client.get_node(uri, force=force).props


def get_status(task, prefix, expnum, version, ccd, return_message=False):
    """
    Report back status of the given program by looking up the associated VOSpace annotation.

    @param task:  name of the process or task that will be checked.
    @param prefix: prefix of the file that was processed (often fk or None)
    @param expnum: which exposure number (or base filename)
    @param version: which version of that exposure (p, s, o)
    @param ccd: which CCD within the exposure.
    @param return_message: Return what did the TAG said or just /True/False/ for Success/Failure?
    @return: the status of the processing based on the annotation value.
    """
    key = get_process_tag(prefix+task, ccd, version)
    status = get_tag(expnum, key)
    logging.debug('%s: %s' % (key, status))
    if return_message:
        return status
    else:
        return status == SUCCESS


def set_status(task, prefix, expnum, version, ccd, status):
    """
    set the processing status of the given program.
    @param task: name of the processing task
    @param prefix: was there a prefix on the exposure number processed?
    @param expnum: exposure number processed.
    @param version: which version of the exposure? (p, s, o)
    @param ccd: the number of the CCD processing.
    @param status: What status to record:  "SUCCESS" we hope.
    @return: Success?
    """
    return set_tag(expnum, get_process_tag(prefix+task, ccd, version), status)


def decompose_content_decomposition(content_decomposition):
    """

    :param content_decomposition:
    :return:
    """
    # check for '-*' in the cutout string and replace is naxis:1
    content_decomposition = re.findall('(\d+)__(\d*)_(\d*)_(\d*)_(\d*)', content_decomposition)
    if len(content_decomposition) == 0:
        content_decomposition = [(0, 1, -1, 1, -1)]
    return content_decomposition


def datasec_to_list(datasec):
    """
    convert an IRAF style PIXEL DATA section as to a list of integers.
    @param datasec: str
    @return: list
    """

    return [int(x) for x in re.findall(r"([-+]?[*\d]+?)[:,\]]+", datasec)]


def reset_datasec(cutout, datasec, naxis1, naxis2):
    """
    reset the datasec to account for a possible cutout.

    @param cutout:
    @param datasec:
    @param naxis1: size of the original image in the 'x' direction
    @param naxis2: size of the oringal image in the 'y' direction
    @return:
    """

    if cutout is None or cutout == "[*,*]":
        return datasec
    try:
        datasec = datasec_to_list(datasec)
    except:
        return datasec

    # check for '-*' in the cutout string and replace is naxis:1
    cutout = cutout.replace(" ", "")
    cutout = cutout.replace("[-*,", "{}:1,".format(naxis1))
    cutout = cutout.replace(",-*]", ",{}:1]".format(naxis2))
    cutout = cutout.replace("[*,", "[1:{},".format(naxis1))
    cutout = cutout.replace(",*]", ",1:{}]".format(naxis1))

    try:
        cutout = datasec_to_list(cutout)
    except:
        logging.debug("Failed to processes the cutout pattern: {}".format(cutout))
        return datasec

    if len(cutout) == 5:
        # cutout likely starts with extension, remove
        cutout = cutout[1:]

    # -ve integer offsets indicate offset from the end of array.
    for idx in [0, 1]:
        if cutout[idx] < 0:
            cutout[idx] = naxis1 - cutout[idx] + 1
    for idx in [2, 3]:
        if cutout[idx] < 0:
            cutout[idx] = naxis2 - cutout[idx] + 1

    flip = cutout[0] > cutout[1]
    flop = cutout[2] > cutout[3]

    logging.debug("Working with cutout: {}".format(cutout))

    if flip:
        cutout = [naxis1 - cutout[0] + 1, naxis1 - cutout[1] + 1, cutout[2], cutout[3]]
        datasec = [naxis1 - datasec[1] + 1, naxis1 - datasec[0] + 1, datasec[2], datasec[3]]

    if flop:
        cutout = [cutout[0], cutout[1], naxis2 - cutout[2] + 1, naxis2 - cutout[3] + 1]
        datasec = [datasec[0], datasec[1], naxis2 - datasec[3] + 1, naxis2 - datasec[2] + 1]

    datasec = [max(datasec[0] - cutout[0] + 1, 1),
               min(datasec[1] - cutout[0] + 1, naxis1),
               max(datasec[2] - cutout[2] + 1, 1),
               min(datasec[3] - cutout[2] + 1, naxis2)]

    return "[{}:{},{}:{}]".format(datasec[0], datasec[1], datasec[2], datasec[3])


class Header(FitsImage):

    def __init__(self, *args, **kwargs):
        if kwargs.get('ext', None) is None:
            kwargs['ext'] = HEADER_EXT
        super(Header, self).__init__(*args, **kwargs)
        self._headers = None

    @property
    def headers(self):
        """
        The a list of the headers assocaited with the MEF.
        :return: List of image headers
        :rtype: list
        """
        if self._headers is not None:
            return self._headers
        if not os.access(self.filename, os.R_OK):
            self.get()

        header_str_list = re.split('END      \n', open(self.filename, 'r').read())

        # make the first entry in the list a Null
        _headers = []
        for header_str in header_str_list[:-1]:
            header = fits.Header.fromstring(header_str, sep='\n')
            if len(_headers) == 0 and not header.get('SIMPLE', False):
                _headers.append({"SOURCE": self.uri})
            _headers.append(header)
        self._headers = _headers
        return self._headers

    @property
    def header(self):
        """
        The FITS header.
        :return: image header
        :rtype: Header
        """
        return self.headers[self.ccd+1]
                            

def make_path(uri):
    """Build a path, with recursion. Don't catch errors."""
    mkdir(os.path.dirname(uri))


def mkdir(dirname):
    """make directory tree in vospace.

    @param dirname: name of the directory to make
    """
    dir_list = []

    while not vospace.client.isdir(dirname):
        logging.debug("Queuing {} for mkdir.".format(dirname))
        dir_list.append(dirname)
        dirname = os.path.dirname(dirname)
    while len(dir_list) > 0:
        logging.info("Creating directory: %s" % (dir_list[-1]))
        try:
            vospace.client.mkdir(dir_list.pop())
        except IOError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise e


def delete(uri):
    vospace.client.delete(uri)


def make_link(source, destination):
    logging.debug("Linking {} to {}".format(source, destination))
    return vospace.client.link(source, destination)


def listdir(directory, force=False):
    return vospace.client.listdir(directory, force=force)


def list_dbimages(dbimages=None):
    if dbimages is None:
        dbimages = DBIMAGES
    return listdir(dbimages)


def exists(uri, force=False):
    try:
        return vospace.client.get_node(uri, force=force) is not None
    except EnvironmentError as e:
        logging.error(str(e))  # not critical enough to raise
        # Sometimes the error code returned is the OS version, sometimes the HTTP version
        if e.errno in [404, os.errno.ENOENT]:
            return False


def move(old_uri, new_uri):
    vospace.client.move(old_uri, new_uri)


def has_property(node_uri, property_name, ossos_base=True):
    """
    Checks if a node in VOSpace has the specified property.

    @param node_uri:
    @param property_name:
    @param ossos_base:
    @return:
    """
    if get_property(node_uri, property_name, ossos_base) is None:
        return False
    else:
        return True


def get_property(node_uri, property_name, ossos_base=True):
    """
    Retrieves the value associated with a property on a node in VOSpace.

    @param node_uri:
    @param property_name:
    @param ossos_base:
    @return:
    """
    # Must use force or we could have a cached copy of the node from before
    # properties of interest were set/updated.
    node = vospace.client.get_node(node_uri, force=True)
    property_uri = tag_uri(property_name) if ossos_base else property_name

    if property_uri not in node.props:
        return None

    return node.props[property_uri]


def set_property(node_uri, property_name, property_value, ossos_base=True):
    """
    Sets the value of a property on a node in VOSpace.  If the property
    already has a value then it is first cleared and then set.

    @param node_uri:
    @param property_name:
    @param property_value:
    @param ossos_base:
    @return:
    """
    node = vospace.client.get_node(node_uri)
    property_uri = tag_uri(property_name) if ossos_base else property_name

    # If there is an existing value, clear it first
    if property_uri in node.props:
        node.props[property_uri] = None
        vospace.client.add_props(node)

    node.props[property_uri] = property_value
    vospace.client.add_props(node)


def log_filename(prefix, task, version, ccd):
    return "{}{}_{}{}{}".format(prefix, task, version, ccd, TEXT_EXT)


def log_location(expnum, ccd):
    return os.path.dirname(Artifact(Observation(expnum), ccd=ccd).uri)


class LoggingManager(object):

    def __init__(self, task, prefix, expnum, ccd, version, dry_run=False):
        self.logging = logging.getLogger('')
        self.log_format = logging.Formatter('%(asctime)s - %(module)s.%(funcName)s %(lineno)d: %(message)s')
        self.filename = log_filename(prefix, task, ccd=ccd, version=version)
        self.location = log_location(expnum, ccd)
        self.dry_run = dry_run

    def __enter__(self):
        if not self.dry_run:
            self.vo_handler = util.VOFileHandler("/".join([self.location, self.filename]))
            self.vo_handler.setFormatter(self.log_format)
            self.logging.addHandler(self.vo_handler)
        self.file_handler = logging.FileHandler(filename=self.filename)
        self.file_handler.setFormatter(self.log_format)
        self.logging.addHandler(self.file_handler)
        return self

    def __exit__(self, *args):
        if not self.dry_run:
            self.logging.removeHandler(self.vo_handler)
            self.vo_handler.close()
            del self.vo_handler
        self.logging.removeHandler(self.file_handler)
        self.file_handler.close()
        del self.file_handler


def archive_url(dataset_name, version, ext=IMAGE_EXT, archive=ARCHIVE, **kwargs):
    url = "{}/{}/{}{}{}".format(DATA_WEB_SERVICE, archive, dataset_name, version, ext)
    if len(kwargs) > 0:
        data = urllib.urlencode(kwargs)
        url += "?"+data
    return url


def pitcairn_uri(dataset_name, ext=None):
    if ext is None:
        ext = IMAGE_EXT
    return "{}/{}{}{}".format(PITCAIRN, dataset_name, PROCESSED_VERSION, ext)


def isfile(uri):
    """
    Check to see if URI points to a VOSpace or URL file.

    :param uri:  The URI to check existence on
    :type uri: str
    :return:
    """
    logging.debug("Checking if {} exists.".format(uri))
    if uri.startswith(VOS_PROTOCOL):
        return vospace.client.isfile(uri)
    try:
        response = vospace.client.conn.session.head(uri)
        response.raise_for_status()
        return True
    except Exception as ex:
        logging.debug(str(ex))
        return False


def copy(source, destination):
    """Copy a file to/from VOSpace. With upto 10 retries on errors"""
    count = 1
    while True:
        try:
            return vospace.client.copy(source, destination)
        except OSError as ex:
            if ex.errno == errno.ENOENT:
                raise ex
            count += 1
            logging.debug(str(ex))
            if count > MAX_RETRY:
                raise ex


def list_exposures(proposal_title='cfis'):
    """
    List all exposures that are part of the project
    :return:
    """

    query = (" SELECT Observation.sequenceNumber AS expnum, "
             " COORD1(CENTROID(Plane.position_bounds)) AS RA, "
             " COORD2(CENTROID(Plane.position_bounds)) AS DE "
             " FROM caom2.Observation AS Observation "
             " JOIN caom2.Plane AS Plane "
             " ON Observation.obsID = Plane.obsID "
             " WHERE Observation.collection = 'CFHT'  "
             " AND Plane.calibrationLevel = 1 "
             " AND lower(Observation.proposal_title) LIKE '%"+proposal_title+"%' "
             " AND Plane.energy_bandpassName LIKE 'r.%'  ")

    return tap_query(query)


def list_healpix():
    exposure_table = list_exposures()
    skycoords = SkyCoord(exposure_table['RA']*units.degree,
                         exposure_table['DE']*units.degree)
    healpix = []
    for coord in skycoords:
        healpix.append(util.skycoord_to_healpix(coord))

    return numpy.unique(healpix)


def tap_query(query):
    """
    Send query to a TAP service and return an astropy table. (could use astroquery instead)
    :param query: SQL to send to the TAP service.
    :type query: str
    :return: Table of results
    :rtype: Table
    """

    data = dict(QUERY=query,
                REQUEST="doQuery",
                LANG="ADQL",
                FORMAT="tsv")

    logging.debug("QUERY: {}".format(data["QUERY"]))
    try:
        result = requests.get(TAP_WEB_SERVICE, params=data, verify=False)
        result.raise_for_status()
    except Exception as ex:
        logging.error(str(ex))

    logging.debug("Doing TAP Query using url: %s" % (str(result.url)))

    table_reader = ascii.get_reader(Reader=ascii.Basic)
    table_reader.header.splitter.delimiter = '\t'
    table_reader.data.splitter.delimiter = '\t'
    table = table_reader.read(result.text)
    return table
