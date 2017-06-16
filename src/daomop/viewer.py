from ginga.web.pgw import ipg, Widgets, Viewers
from astropy.wcs import WCS
import storage
import mp_ephem
import logging


logging.basicConfig(level=logging.INFO)


class GoodGui(ipg.EnhancedCanvasView):
    def build_gui(self, container):
        vbox = Widgets.VBox()
        vbox.set_border_width(2)
        vbox.set_spacing(1)

        w = Viewers.GingaViewerWidget(viewer=self)
        vbox.add_widget(w, stretch=1)

        self.pixel_base = 1.0

        self.readout = Widgets.Label("")
        vbox.add_widget(self.readout, stretch=0)

        # self.set_callback('none-move', self.motion_cb)
        self.set_callback('cursor-changed', self.motion_cb)
        accept = Widgets.Button("Accept")
        accept.add_callback('activated', lambda x: self.onscreen_message("Accept", 1))

        reject = Widgets.Button("Reject")
        reject.add_callback('activated', lambda x: self.onscreen_message("Reject", 1))

        wclear = Widgets.Button("Clear")
        wclear.add_callback('activated', lambda x: self.clear())

        vbox.add_widget(accept, stretch=1)
        vbox.add_widget(reject, stretch=1)
        vbox.add_widget(wclear, stretch=1)

        # need to put this in an hbox with an expanding label or the
        # browser wants to resize the canvas, distorting it
        hbox = Widgets.HBox()
        hbox.add_widget(vbox, stretch=0)
        hbox.add_widget(Widgets.Label(''), stretch=1)

        container.set_widget(hbox)


class Viewer(object):

    def _write(self):

        with open('filename', 'w') as fobj:
            for ob in self.obs:
                fobj.write(ob.to_string())

    def __init__(self, target, dbimages="vos:jkavelaars/TNORecon/dbimages"):
        """
        Initialization of a local, web-client based server for displaying images. The viewer should automatically pop
         up in a new tab.

        :param target: Filename to read for creating ObsRecord objects
        :param dbimages: Directory location for images
        """

        # standard setup commands; creating viewing window
        self.server = ipg.make_server(host='localhost', port=9914, use_opencv=False, viewer_class=GoodGui)
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
        self._current_image = None

    def load(self, obs_number=0):
        """
        With the viewing window already created, Creates a FitsImage object and loads its cutout into the window.

        :param obs_number: index of which line in the file gets loaded/displayed in the viewer
        """
        self.obs_number = obs_number
        self._load()

    def _mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object
         being observed.
        """
        ra, dec = self.obs[self.obs_number].coordinate.ra, self.obs[self.obs_number].coordinate.dec
        x, y = WCS(self._loaded_images[self.obs_number].header).all_world2pix(ra, dec, 0)
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
        if keyname == 'f':
            # if/else statements for blinking. Only needs to grab the hdu from the database
            #  once, after that it saves the hdu so it can be quickly loaded into the viewer again.
            self.obs_number -= 1
        elif keyname == 'g':
            self.obs_number += 1
        elif keyname == 'z':
            pass                # soon to be accept/reject shortcuts
        elif keyname == 'c':
            self.accept()
        else:
            logging.warning("Unknown keystroke {}".format(keyname))
            return

        self.obs_number %= len(self.obs)
        self._load()

    def _load(self):
        """
        Checks if an HDU has been loaded already and retrieves if needed and then displays that HDU.
        Calls _mark_aperture to draw a small circle around the celestial object.
        """
        # load the image if not already available, for now we'll put this in here.
        self.viewer.clear()
        if self._loaded_images[self.obs_number] is None:
            self.get_hdu(self.obs_number)

        # pv_to_sip(self._loaded_images[self.obs_number].header)
        # logging.info(self._loaded_images[self.obs_number].header)

        self.viewer.load_hdu(self._loaded_images[self.obs_number])
        self._mark_aperture()
        self.viewer.onscreen_message("Loaded: {}".format(self.obs[self.obs_number].comment.frame), delay=3)

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

    def accept(self):
        logging.warning("Accept keypress")
        return

    def __del__(self):
        self.server.stop()

    # def __enter__(self):
    #     return self
    #
    # def __exit__(self, exc_type, exc_val, exc_tb):
    #     self.server.stop()
