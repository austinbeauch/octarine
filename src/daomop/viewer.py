from ginga.web.pgw import ipg
import storage
from astropy.wcs import WCS
import wcs
import mp_ephem
import logging

WCS.all_pix2world = wcs.WCS.all_pix2world

logging.basicConfig(level=logging.INFO)


class Viewer(object):

    def __init__(self, target, dbimages="vos:jkavelaars/TNORecon/dbimages"):
        """
        Initialization of a local, web-client based server for displaying images. The viewer should automatically pop
         up in a new tab.

        :param target: Filename to read for creating ObsRecord objects
        :param dbimages: Directory location for images
        """

        # standard setup commands; creating viewing window
        self.server = ipg.make_server(host='localhost', port=9914, use_opencv=False)
        self.server.start(no_ioloop=True)
        self.viewer = self.server.get_viewer('v1')
        self.viewer.enable_autocuts('on')
        self.viewer.set_autocut_params('zscale')
        self.viewer.open()

        # creating drawing canvas; initializing polygon types
        self.canvas = self.viewer.add_canvas()
        self.circle = self.canvas.get_draw_class('circle')
        self.ellipse = self.canvas.get_draw_class('ellipse')

        # creating key-press event handling
        self.canvas.add_callback('key-press', self._key_press, 'key', self.viewer)

        # reading target file; setting database location; creating objects
        storage.DBIMAGES = dbimages
        self.obs = mp_ephem.EphemerisReader().read(target)
        self._loaded_images = [None] * len(self.obs)
        logging.info("Computing preliminary orbit.")
        self.orb = mp_ephem.BKOrbit(self.obs)
        logging.info(self.orb.summarize())
        self.obs_number = 0
        self.image = None
        self.image_values = None

        # eventually store different HDU's for image switching,
        # probably will be replaced with a single index counter variable in _key_press
        self._current_image = None

    def load(self, obs_number=0):
        """
        With the viewing window already created, Creates a FitsImage object and loads its cutout into the window.
        Calls _mark_aperture to draw a small circle around the celestial object pertaining to that cutout.

        :param obs_number: index of which line in the file gets loaded/displayed in the viewer
        """
        self.obs_number = obs_number
        self._load()

        # Attempt to monkey patch all_pix2world, doesn't seem to be getting called before keyword error is thrown.
        # Warnings probably being thrown before pix2world is called.

    def _mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object
         being observed.
        """
        ra, dec = self.obs[self.obs_number].coordinate.ra, self.obs[self.obs_number].coordinate.dec
        x, y = wcs.WCS(self._loaded_images[self.obs_number].header).sky2xy(ra, dec, usepv=True)
        self.canvas.deleteAllObjects()
        self.canvas.add(self.circle(x, y, radius=10, color='red'))

    def _predict_obs_record(self):
        """
        Creates a prediction from the ObsRecord object using mp_ephem.BKOrbit's predict method.
        Creates a circle on the canvas which is where the celestial object is predicted to be, along with its
         associated uncertainty, also from BKOrbit's predict method.
        """
        self.orb.predict(self.obs[self.obs_number].date)

        # x, y = self.orb.coordinate.ra.deg, self.orb.coordinate.dec.deg
        # tmp = self.image_values.radectopix(x, y)  # not working, hopefully will once wcs conversion is fixed

        tmp = (294, 318)  # placeholder values
        x, y = tmp[0], tmp[1]

        # prediction circle
        self.canvas.add(self.circle(x, y, radius=10, color='green'))

        # uncertainty ellipse
        # Not sure about changing the ra/dec uncertainty from acrseconds into pix values.
        # Current values in acrsec/degrees are on the magnitude of 10^-5, which doesn't do much as far as creating
        #  an ellipse goes.
        self.canvas.add(self.ellipse(x,
                                     y,
                                     self.orb.dra.to('deg').value,
                                     self.orb.ddec.to('deg').value, color='red'))

    # button press methods NEED all these arguments to work even if they are not used.
    # Might be a weird Ginga thing.
    def _key_press(self, canvas, keyname, opn, viewer):
        """
        Method called once a keyboard stoke has been detected. Using two un-bound keys, f & g, to cycle different
         cutout hdu's from the ObsRecord.
        Parameters canvas, opn, and viewer are all needed fo the method to be called even though they are not
         directly used in the method. Method appears to be in connection with the add_callback method call from
         __init__, but isn't used directly.

        :param canvas: Ginga DrawingCanvas Object
        :param keyname: Name of the key that has been pressed
        :param opn: str "key"
        :param viewer: Ginga EnhancedCanvasView object
        """
        logging.debug("Got key: {} from canvas: {} with opn: {} from viewer: ".format(canvas, keyname, opn, viewer))
        if keyname == 'f':  # index 0
            self.viewer.set_onscreen_message("pressed f")  # TODO: remove onscreen messages once finished implementation

            # if/else statements to speed up image switching. Only needs to grab the hdu from the database
            #  once, after that it saves the hdu so it can be quickly loaded into the viewer again.
            # Currently only switching between 2 images. Need more information as to how other .ast files are set up
            #  before implementing an iterative solution to image switching.
            # Is there a quantity variable storing the amount of lines that can be iterated?
            # In the case that there is an amount variable, set up a simple counter to act as the index,
            #  increment or decrement appropriately, loading each image into the viewer, checking to make sure it's not
            #  going out of range. If the user is about to go below 0 or above the max range, just do nothing.
            #  Probably don't need any warning, but it would be possible to have an onscreen message pop up.
            # Retrieving the image cutouts definitely takes more time than needed, not sure if anything in storage.py
            #  can be sped up or if it's a server issue.
            self.obs_number += 1
        elif keyname == 'g':
            self.obs_number -= 1
        else:
            logging.warning("Unknown keystroke {}".format(keyname))
            return

        self.obs_number %= len(self.obs)
        self._load()

    def _load(self):
        """
        Checks if an HDU has been loaded already and retrieves if needed and then displays that HDU.
        """
        # load the image if not already available, for now we'll put this in here.
        self.viewer.clear()
        if self._loaded_images[self.obs_number] is None:
            self.get_hdu(self.obs_number)
        self.viewer.load_hdu(self._loaded_images[self.obs_number])
        self._mark_aperture()
        self.viewer.set_onscreen_message("Loaded: {}".format(self.obs[self.obs_number].comment.frame))

    def get_hdu(self, obs_number):
        """
        Retrieve a fits image associated with a given obs_record and return the HDU off the assocaited cutout.

        :param obs_number: the number of the record to retrieve.
        :type obs_number: int
        :return:
        """
        obs_record = self.obs[obs_number]
        image = storage.FitsImage.from_frame(obs_record.comment.frame)
        hdu = image.ra_dec_cutout(obs_record.coordinate)[-1]  # always index 2?
        self._loaded_images[obs_number] = hdu
        logging.info("Retrieved {} and stored {}".format(obs_number, self._loaded_images))
        return

    def __del__(self):
        self.server.stop()

    # def __enter__(self):
    #     return self
    #
    # def __exit__(self, exc_type, exc_val, exc_tb):
    #     self.server.stop()
