import logging
from math import atan2, degrees

import candidate
from ginga.web.pgw import ipg, Widgets, Viewers
from astropy.wcs import WCS

dbimages = "vos:jkavelaars/TNORecon/dbimages"
logging.basicConfig(level=logging.INFO)


class ValidateGui(ipg.EnhancedCanvasView):

    def build_gui(self, container):
        """
        Building the GUI to be displayed in an HTML5 canvas. Currently consists of three buttons and a text box.
        Tested and working in Mozilla Firefox web browser.
        :param container: ginga.web.pgw.Widgets.TopLevel object
        """
        self.candidates = None
        self.load = None
        self.write_record = None
        self.skip = None

        vbox = Widgets.VBox()
        vbox.set_border_width(2)
        vbox.set_spacing(1)

        w = Viewers.GingaViewerWidget(viewer=self)
        vbox.add_widget(w, stretch=1)

        self.pixel_base = 1.0
        self.readout = Widgets.Label("")

        vbox.add_widget(self.readout, stretch=0)

        self.set_callback('cursor-changed', self.motion_cb)

        load_candidates = Widgets.TextEntry()
        load_candidates.add_callback('activated', lambda x: self.load_candiates(x))

        accept = Widgets.Button("Accept")
        accept.add_callback('activated', lambda x: self.accept())

        reject = Widgets.Button("Reject")
        reject.add_callback('activated', lambda x: self.reject())

        wclear = Widgets.Button("Skip")
        wclear.add_callback('activated', lambda x: self.jump())

        hbox = Widgets.HBox()
        h2box = Widgets.HBox()
        h3box = Widgets.HBox()

        h2box.add_widget(accept, stretch=0)
        h2box.add_widget(reject, stretch=0)
        h2box.add_widget(wclear, stretch=0)
        h2box.set_spacing(5)
        h2box.set_margins(0, 0, 25, 0)
        h3box.add_widget(load_candidates, stretch=1)

        vbox.add_widget(h2box, stretch=1)
        vbox.add_widget(h3box, stretch=1)

        # need to put this in an hbox with an expanding label or the
        # browser wants to resize the canvas, distorting it
        hbox.add_widget(vbox, stretch=0)
        hbox.add_widget(Widgets.Label(''), stretch=1)

        container.set_widget(hbox)

    def jump(self):
        self.onscreen_message("Skipping", 1)

        # this didn't work for whatever reason
        # self.candidates.next()
        # if self.load is not None:
        #     self.load()

        if self.skip is not None:
            self.skip()

    def reject(self):
        self.onscreen_message("Rejected", delay=1)
        if self.write_record is not None:
            self.write_record(rejected=True)

    def accept(self):
        self.onscreen_message("Accepted", delay=1)
        if self.write_record is not None:
            self.write_record()

    def load_candiates(self, event):
        self.onscreen_message("Entered: {}".format(event.text), 1)
        logging.info("Accepted candidate entry: {}".format(event.text))

        self.candidates = candidate.CandidateSet(int(event.text))
        if self.load is not None:
            self.load()


class WebServerFactory(object):
    """
    The Server that the validate app will be run via.
    """

    def __enter__(self):
        self.web_server = ipg.make_server(host=self.host,
                                          port=self.port,
                                          use_opencv=False,
                                          viewer_class=self.viewer_class)
        self.web_server.start(no_ioloop=True)
        return self.web_server

    def __init__(self, viewer_class=ValidateGui, port=9914, host='localhost'):
        # standard setup commands; creating viewing window
        self.viewer_class = viewer_class
        self.port = port
        self.host = host

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.web_server.stop()

    def __del__(self):
        self.web_server.stop()


class ImageViewer(object):

    def __init__(self, downloader, *args, **kwargs):
        """
        Initialization of a local, web-client based server for displaying images. The viewer should automatically pop
         up in a new tab.

        :param candidate: A BK Orbfit object that contains the information about a possible candidate.
        :type candidate: BKOrbit
        :param viewer: The Viewer that will be used to display images.
        :type viewer: ipg.EnhancedCanvasViewer
        """
        self.pixel_base = 1
        self.readout = Widgets.Label("")

        # standard setup commands; creating viewing window
        self.web_server = ipg.make_server(host="localhost",
                                          port=9914,
                                          use_opencv=False,
                                          viewer_class=ValidateGui)
        self.web_server.start(no_ioloop=True)
        self.viewer = self.web_server.get_viewer("ID")
        self.viewer.enable_autocuts('on')
        self.viewer.set_autocut_params('zscale')
        self.viewer.open()

        self.viewer.load = self.load
        self.viewer.write_record = self.write_record
        self.viewer.skip = self.skip

        self.downloader = downloader

        # creating drawing canvas; initializing polygon types
        self.canvas = self.viewer.add_canvas()
        self.circle = self.canvas.get_draw_class('circle')
        self.ellipse = self.canvas.get_draw_class('ellipse')

        # creating key-press event handling
        self.canvas.add_callback('key-press', self._key_press, 'key', self.viewer)

        self.obs_number = 0
        self.candidate = None
        self.pan = None
        self.zoom = None
        self._center = None

    def write_record(self, rejected=False):
        try:
            with open(self.candidate.observations[0].provisional_name+".ast", 'w+') as fobj:
                for ob in self.candidate.observations:
                    if rejected:
                        ob.null_observation = True
                    fobj.write(ob.to_string()+'\n')
            logging.info("Written to file {}".format(self.candidate.observations[0].provisional_name+".ast"))
        except IOError as ex:
            logging.error("Unable to write to file.")
            raise ex
        self.skip()

    def skip(self):
        """
        Skipping to the next candidate set without writing to a file.
        """
        # so this is where I realize there must be a better way
        self.candidate = self.viewer.candidates.next()
        self.load()

    def load(self, obs_number=0):
        """
        With the viewing window already created, Creates a FitsImage object and loads its cutout into the window.
        Define the center of the first image to be the reference point for aligning the other two images.
        :param obs_number: index of which line in the file gets loaded/displayed in the viewer
        """
        self._center = None
        self.obs_number = obs_number
        self._load()
        self._center = WCS(self.header).all_pix2world(self.viewer.get_data_size()[0] / 2,
                                                      self.viewer.get_data_size()[1] / 2, 0)

    def _load(self):
        """
        Checks if an HDU has been loaded already and retrieves if needed and then displays that HDU.
        Calls _mark_aperture to draw a small circle around the celestial object.
        Calls _rotate to orient the image North up East left
        """
        # load the image if not already available, for now we'll put this in here.
        print self.obs_number
        if self.viewer.candidates is None:
            logging.debug("No candidates loaded.")
            return

        if self.candidate is None:
            self.candidate = self.viewer.candidates.next()

        self.viewer.clear()
        self.viewer.load_hdu(self.downloader.get(self.candidate.observations[self.obs_number]))

        if (self.pan and self.zoom) is not None:
            self.set_position()

        self._mark_aperture()
        self._rotate()

        if self.center is not None:
            self._align()

        self.viewer.onscreen_message("Loaded: {}".format(self.candidate.observations[self.obs_number].comment.frame),
                                     delay=1)

    def _mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object
         being observed.
        """
        ra = self.candidate.observations[self.obs_number].coordinate.ra
        dec = self.candidate.observations[self.obs_number].coordinate.dec
        x, y = WCS(self.header).all_world2pix(ra, dec, 0)
        self.canvas.deleteAllObjects()
        self.canvas.add(self.circle(x, y, radius=10, color='red'))

    def _rotate(self):
        """
        Rotates the current viewer image to be oriented North up East left. This is done by taking outward vectors from
         the origin and using their WCS values to determine the original orientation of the image. Images are then
         flipped/rotated accordingly to be North up East left.
        """
        wcs = WCS(self.header)
        self.viewer.transform(False, False, False)
        x = wcs.all_pix2world([[0, 0], [1, 1], [1, 0]], 0)
        ra1 = x[0][0]
        ra2 = x[1][0]
        ra3 = x[2][0]
        dec1 = x[0][1]
        dec2 = x[1][1]
        dec3 = x[2][1]

        # phi = acos((ra2-ra1)/((ra2-ra1)**2+(dec2-dec1)**2)**0.5)

        delta_x = ra2 - ra1
        delta_y = dec2 - dec1

        flip_x = 1
        flip_y = 1
        if not delta_x < 0:
            flip_x = -1
            if not delta_y > 0:
                flip_y = -1
                self.viewer.transform(True, True, False)  # def transform(self, flip_x, flip_y, swap_xy):
            else:
                self.viewer.transform(True, False, False)
        elif not delta_y > 0:
            flip_y = -1
            self.viewer.transform(False, True, False)

        delta_delta = (dec3 - dec1) * flip_y
        delta_ra = (ra3 - ra1) * flip_x * -1

        theta = degrees(atan2(delta_delta, delta_ra))

        self.viewer.rotate(theta)

    def _align(self):
        """
        Aligns images via panning so their backgrounds stay consistent. Images requiring a pan greater than 1/2 the
         viewing window will be ignored.
        """
        try:
            x, y = WCS(self.header).all_world2pix(self.center[0], self.center[1], 0)

            if not(0 < x < self.viewer.get_data_size()[0] and 0 < y < self.viewer.get_data_size()[1]):
                logging.info("Pan out of range: ({}, {}) greater than half the viewing window.".format(x, y))
                return

            self.viewer.set_pan(x, y)

        except Exception as ex:
            logging.warning("Could not convert ra/dec to pixel values")
            raise ex

    # def _predict_obs_record(self):
    #     """
    #     Creates a prediction from the ObsRecord object using mp_ephem.BKOrbit's predict method.
    #     Creates a circle on the canvas which is where the celestial object is predicted to be, along with its
    #      associated uncertainty, also from BKOrbit's predict method.
    #     """
    #     self.candidate.predict(self.candidate[self.obs_number].date)
    #
    #     # x, y = self.orb.coordinate.ra.deg, self.orb.coordinate.dec.deg
    #     # tmp = self.image_values.radectopix(x, y)
    #
    #     tmp = (294, 318)  # placeholder values
    #     x, y = tmp[0], tmp[1]
    #
    #     # prediction circle
    #     self.canvas.add(self.circle(x, y, radius=10, color='green'))
    #
    #     # uncertainty ellipse
    #     # Not sure about changing the ra/dec uncertainty from acrseconds into pix values.
    #     # Current values in acrsec/degrees are on the magnitude of 10^-5, which doesn't do much as far as creating
    #     #  an ellipse goes.
    #     self.canvas.add(self.ellipse(x,
    #                                  y,
    #                                  self.candidate.dra.to('deg').value,
    #                                  self.candidate.ddec.to('deg').value, color='red'))

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
            self.get_position()
            self.obs_number -= 1

        elif keyname == 'g':
            self.get_position()
            self.obs_number += 1

        else:
            logging.debug("Unknown keystroke {}".format(keyname))  # keystrokes might be for ginga
            return

        self.obs_number %= len(self.candidate.observations)
        self._load()

    @property
    def center(self):
        """
        Returns the center of the image in ra/dec coordinates
        """
        if self._center is not None:
            return self._center

    @property
    def loaded_hdu(self):
        """
        Return current HDU
        """
        return self.downloader.get(self.candidate.observations[self.obs_number])

    @property
    def header(self):
        """
        Return current HDU's header
        """
        return self.loaded_hdu.header

    def get_position(self):
        """
        Stores the current image's pan location and zoom amount
        """
        self.pan = self.viewer.get_pan()
        self.zoom = self.viewer.get_zoom()

    def set_position(self):
        """
        Sets the current image's pan and zoom to what was saved from get_position
        """
        # self.viewer.set_pan(self.pan[0], self.pan[1]) TODO: make pans relative depending on orientation
        self.viewer.zoom_to(self.zoom)
