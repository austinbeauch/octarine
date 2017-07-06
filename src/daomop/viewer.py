import sys
import logging
from multiprocessing.dummy import Pool, Lock
from math import atan2, degrees
from multiprocessing.pool import ApplyResult

import candidate
import downloader
import storage
from ginga import AstroImage
from ginga.web.pgw import ipg, Widgets, Viewers
from ginga.misc import log
from astropy.wcs import WCS

logging.basicConfig(level=logging.INFO, format="%(module)s.%(funcName)s:%(lineno)s %(message)s")
DISPLAY_KEYWORDS = ['EXPNUM', 'DATE-OBS', 'UTC-OBS', 'EXPTIME', 'FILTER']
LEGEND = 'Keyboard Shortcuts: \n' \
         'f: image backwards \n' \
         'g: image forwards \n' \
         't: contrast mode'
PROCESSES = 5


class ValidateGui(ipg.EnhancedCanvasView):

    def __init__(self, logger, window, bindings=None):
        super(ValidateGui, self).__init__(logger=logger, bindings=bindings)

        self.console_box = Widgets.TextArea(editable=False)

        # self.console_streamer = logging.StreamHandler(stream=self.console_box.append_text)

        self.downloader = downloader.Downloader()
        self.pool = Pool(processes=PROCESSES)
        self.lock = Lock()
        self.image_list = {}
        self.astro_images = {}

        self.logger = logger
        self.top = window

        self.enable_autocuts('on')
        self.set_autocut_params('zscale')

        # creating drawing canvas; initializing polygon types
        self.canvas = self.add_canvas()
        self.circle = self.canvas.get_draw_class('circle')

        # creating key-press event handling
        self.canvas.add_callback('key-press', self._key_press, 'key', self)

        self.obs_number = 0
        self.candidate = None
        self.candidates = None
        self.zoom = None
        self._center = None
        self.event = None

        # GUI elements
        self.pixel_base = 1.0
        self.readout = Widgets.Label("")
        self.header_box = Widgets.TextArea(editable=False)

        self.legend = Widgets.TextArea(wrap=True)
        self.legend.set_text(LEGEND)
        self.build_gui(self.top)

    def build_gui(self, container):
        """
        Building the GUI to be displayed in an HTML5 canvas.
        Tested and working in Mozilla Firefox and Google Chrome web browsers.

        :param container: ginga.web.pgw.Widgets.TopLevel object
        """
        bindings = self.get_bindings()
        bindings.enable_all(True)

        # keyboard mode indicator, upper right corner
        self.show_mode_indicator(True, corner='ur')

        viewer_vbox = Widgets.VBox()  # box containing the viewer
        viewer_vbox.set_border_width(2)
        viewer_vbox.set_spacing(1)
        viewer_widget = Viewers.GingaViewerWidget(viewer=self)
        viewer_vbox.add_widget(viewer_widget, stretch=1)
        viewer_vbox.add_widget(self.readout, stretch=0)  # text directly below the viewer for coordinate display

        self.set_callback('cursor-changed', self.motion_cb)

        load_candidates = Widgets.TextEntry()
        load_candidates.add_callback('activated', lambda x: self.load_candidates(x))

        accept = Widgets.Button("Accept")
        accept.add_callback('activated', lambda x: self.accept_reject())

        reject = Widgets.Button("Reject")
        reject.add_callback('activated', lambda x: self.accept_reject(rejected=True))

        next_set = Widgets.Button("Next")
        next_set.add_callback('activated', lambda x: self.next())

        previous_set = Widgets.Button("Previous")
        previous_set.add_callback('activated', lambda x: self.previous())

        quit_button = Widgets.Button("Quit")
        quit_button.add_callback('activated', lambda x: self.exit())

        reload_button = Widgets.Button("Reload")
        reload_button.add_callback('activated', lambda x: self.reload_candidates())

        # accept/reject/next buttons
        buttons_hbox = Widgets.HBox()
        buttons_hbox.add_widget(accept)
        buttons_hbox.add_widget(reject)
        buttons_hbox.add_widget(previous_set)
        buttons_hbox.add_widget(next_set)
        buttons_hbox.add_widget(load_candidates)
        buttons_hbox.set_spacing(7)
        viewer_vbox.add_widget(buttons_hbox)  # add buttons below the viewer

        # quit/reload buttons
        quit_box = Widgets.HBox()
        quit_box.add_widget(quit_button)
        quit_box.add_widget(reload_button)

        viewer_header_hbox = Widgets.HBox()  # box containing the viewer/buttons and rightmost text area
        viewer_header_hbox.add_widget(viewer_vbox)
        viewer_header_hbox.add_widget(Widgets.Label(''))
        viewer_header_hbox.add_widget(self.header_box)
        viewer_header_hbox.add_widget(self.legend)

        full_vbox = Widgets.VBox()  # vbox container for all elements
        full_vbox.add_widget(viewer_header_hbox)

        full_vbox.add_widget(self.console_box)
        full_vbox.add_widget(quit_box)
        self.console_box.set_text('Logging output:\n')
        container.set_widget(full_vbox)

    def next(self):
        """
        Load the next set of images into the viewer
        """
        if self.candidates is not None:
            self.obs_number = 0
            self.candidate = self.candidates.next()
            self.load()

    def previous(self):
        """
        Load the previous set of images into the viewer
        """
        if self.candidates is not None:
            self.obs_number = 0
            self.candidate = self.candidates.previous()
            self.load()

    def accept_reject(self, rejected=False):
        """
        Accept or reject current observation depending on button press. Write to file and load next set into the viewer

        :param rejected: whether the candidate set has been accepted or rejected
        """
        if self.candidates is not None:
            self.write_record(rejected=rejected)
            self.next()

    def exit(self):
        """
        Shuts down the application
        """
        self.readout.set_text("Shutting down application.")
        self.logger.info("Attempting to shut down the application...")
        if self.top is not None:
            self.top.close()
        sys.exit()

    def load_candidates(self, event):
        """
        Initial candidates loaded into the viewer. Starts up a threadpool to download images simultaneously.

        :param event: Catalogue number containing dataset
        """
        if hasattr(event, 'text'):
            self.event = int(event.text)

        self.readout.set_text("Accepted candidate entry: {}".format(self.event))
        self.candidates = candidate.CandidateSet(self.event)

        with self.lock:
            for bk_orbit in self.candidates:
                for obs_record in bk_orbit.observations:
                    key = self.downloader.image_key(obs_record)
                    self.image_list[key] = self.pool.apply_async(self.downloader.get, (obs_record,))

        self.candidates = candidate.CandidateSet(self.event)
        self.load()

    def reload_candidates(self):
        """
        Performs a hard reload on all images for the case of loading errors.
        Closes current worker pool and reopens a new one.
        """
        if self.event is not None:
            self.console_box.append_text('Reloading all canditates...\n')
            self.pool.terminate()
            self.pool = Pool(processes=PROCESSES)
            self.load_candidates(self.event)
            self.next()

    def load(self, obs_number=0):
        """
        With the viewing window already created, Creates a FitsImage object and loads its cutout into the window and
         displays select header values (see: DISPLAY_KEYWORDS).
        Define the center of the first image to be the reference point for aligning the other two images in the set.

        :param obs_number: index of which line in the file gets loaded/displayed in the viewer
        """
        self._center = None
        self.obs_number = obs_number
        self._load()
        self._center = WCS(self.header).all_pix2world(self.get_data_size()[0] / 2,
                                                      self.get_data_size()[1] / 2, 0)

    def _load(self):
        """
        Loads an image into the viewer, applying appropriate transformations for proper display.
        Checks if an HDU has been loaded already and retrieves if needed and then displays that HDU.
        Uses multiprocessing techniques for simultaneous downloads and dictionaries to keep track of which images
         have been already loaded for faster image switching.
        """
        # load the image if not already available, for now we'll put this in here.
        if self.candidates is None:
            logging.debug("No candidates loaded.")
            return

        if self.candidate is None:
            self.next()

        while True:
            # noinspection PyBroadException
            try:
                key = self.downloader.image_key(self.candidate.observations[self.obs_number])

                with self.lock:
                    hdu = (isinstance(self.image_list[key], ApplyResult) and self.image_list[key].get()
                           or self.image_list[key])
                    self.image_list[key] = hdu

                if key not in self.astro_images:
                    image = AstroImage.AstroImage(logger=self.logger)
                    if len(self.image_list[key]) > 1:
                        image.load_hdu(self.image_list[key][2])
                    else:
                        image.load_hdu(self.image_list[key][0])

                    self.astro_images[key] = image

                self.set_image(self.astro_images[key])
                break

            except Exception as ex:
                logging.error(str(ex))
                self.console_box.append_text(str(ex) + '\n')
                logging.debug("Skipping candidate {} due to load failure".format(self.obs_number))
                self.next()

        if self.zoom is not None:
            self.zoom_to(self.zoom)

        self._rotate()

        if self.center is not None:
            self._align()

        self._mark_aperture()

        self.header_box.set_text(self.info)
        self.onscreen_message("Loaded: {}".format(self.candidate.observations[self.obs_number].comment.frame),
                              delay=0.5)

    def write_record(self, rejected=False):
        """
        Writing observation lines to a new file.

        :param rejected: Whether or not the candidate set contains a valid celestial object
        :type rejected: bool
        """
        try:
            art = storage.ASTRecord(self.candidate.observations[0].provisional_name)
            with open(art.filename, 'w+') as fobj:
                for ob in self.candidate.observations:
                    if rejected:
                        ob.null_observation = True
                    fobj.write(ob.to_string() + '\n')
            art.put()
            logging.info("Written to file {}".format(self.candidate.observations[0].provisional_name+".ast"))
        except IOError as ex:
            logging.error("Unable to write to file.")
            raise ex

    def _mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object detected.
        """
        # the image cutout is considered the first object on the canvas, this deletes everything over top of it
        self.canvas.delete_objects(self.canvas.get_objects()[1:])

        ra = self.candidate.observations[self.obs_number].coordinate.ra
        dec = self.candidate.observations[self.obs_number].coordinate.dec
        x, y = WCS(self.header).all_world2pix(ra, dec, 0)
        self.canvas.add(self.circle(x, y, radius=10, color='red'))

    def _rotate(self):
        """
        Rotates the current viewer image to be oriented North up East left. This is done by taking outward vectors from
         the origin and using their WCS values to determine the original orientation of the image. Images are then
         flipped/rotated accordingly to be North up East left.
        """
        wcs = WCS(self.header)
        self.transform(False, False, False)
        x = wcs.all_pix2world([[0, 0], [1, 1], [1, 0]], 0)
        ra1 = x[0][0]
        ra2 = x[1][0]
        ra3 = x[2][0]
        dec1 = x[0][1]
        dec2 = x[1][1]
        dec3 = x[2][1]

        delta_x = ra2 - ra1
        delta_y = dec2 - dec1

        flip_x = 1
        flip_y = 1
        if not delta_x < 0:
            flip_x = -1
            if not delta_y > 0:
                flip_y = -1
                self.transform(True, True, False)  # def transform(self, flip_x, flip_y, swap_xy):
            else:
                self.transform(True, False, False)
        elif not delta_y > 0:
            flip_y = -1
            self.transform(False, True, False)

        delta_delta = (dec3 - dec1) * flip_y
        delta_ra = (ra1 - ra3) * flip_x

        theta = degrees(atan2(delta_delta, delta_ra))

        self.rotate(theta)

    def _align(self):
        """
        Aligns images via panning so their backgrounds stay consistent. Images requiring a pan greater than 1/2 the
         viewing window will be ignored.
        """
        x, y = WCS(self.header).all_world2pix(self.center[0], self.center[1], 0)

        if not(0 < x < self.get_data_size()[0] and 0 < y < self.get_data_size()[1]):
            logging.debug("Pan out of range: ({}, {}) is greater than half the viewing window.".format(x, y))
            self.console_box.append_text("Pan out of range: ({}, {}) is greater than half the viewing window."
                                         .format(x, y) + '\n')
            return

        self.set_pan(x, y)

    def _key_press(self, canvas, keyname, opn, viewer):
        """
        Method called once a keyboard stoke has been detected. Using two un-bound keys, f & g, to cycle different
         cutout hdu's from the ObsRecord.
        Parameters canvas, opn, and viewer are all needed for the method to be called even though they are not
         directly used.

        :param canvas: Ginga DrawingCanvas Object
        :param keyname: Name of the key that has been pressed
        :param opn: str "key"
        :param viewer: Ginga EnhancedCanvasView object
        """
        self.logger.info("Got key: {} from canvas: {} with opn: {} from viewer: {}".format(canvas,
                                                                                           keyname,
                                                                                           opn,
                                                                                           viewer))
        if keyname == 'f':
            self.obs_number -= 1

        elif keyname == 'g':
            self.obs_number += 1

        self.zoom = self.get_zoom()
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
        hdu_list = self.downloader.get(self.candidate.observations[self.obs_number])

        if len(hdu_list) > 1:
            return hdu_list[1]  # Return index 1 if multi extension file? Possibly. Could use index 1 as reference.

        return hdu_list[0]

    @property
    def header(self):
        """
        Return current HDU's header
        """
        return self.loaded_hdu.header

    @property
    def info(self):
        return "\n".join([x + " = " + str(self.header.get(x, "UNKNOWN")) for x in DISPLAY_KEYWORDS])


class ImageViewer(object):

    def __init__(self):
        """
        Initialization of a local, web-client based server for displaying images. The viewer should automatically pop
         up in a new tab.
        """
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


def main(params):

    logger = log.get_logger("daomop", options=params)

    if params.use_opencv:
        from ginga import trcalc
        try:
            trcalc.use('opencv')
        except Exception as ex:
            logger.warning("Error using OpenCL: {}".format(ex))

    if params.use_opencl:
        from ginga import trcalc
        try:
            trcalc.use('opencl')
        except Exception as ex:
            logger.warning("Error using OpenCL: {}".format(ex))

    app = Widgets.Application(logger=logger, host=params.host, port=params.port)

    #  create top level window
    window = app.make_window("Validate 0.2.0", wid='Validate')

    # our own viewer object, customized with methods (see above)
    ValidateGui(logger, window)

    try:
        app.start()

    except KeyboardInterrupt:
        logger.info("Terminating viewer...")
        window.close()
