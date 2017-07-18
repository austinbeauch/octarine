import os
import re
import logging
from multiprocessing.dummy import Pool, Lock
from math import atan2, degrees
from multiprocessing.pool import ApplyResult

from astropy.time import Time
from mp_ephem import ObsRecord

from . import candidate
from . import downloader
from . import storage
from ginga import AstroImage
from ginga.web.pgw import ipg, Widgets, Viewers
from ginga.misc import log
from astropy.wcs import WCS

logging.basicConfig(level=logging.INFO, format="%(module)s.%(funcName)s:%(lineno)s %(message)s")
DISPLAY_KEYWORDS = ['EXPNUM', 'DATE-OBS', 'UTC-OBS', 'EXPTIME', 'FILTER']
LEGEND = 'Keyboard Shortcuts: \n' \
         'f: image backwards \n' \
         'g: image forwards \n' \
         'q: pan mode \n(click and drag on canvas)\n' \
         't: contrast mode \n(right click on canvas after pressing "t" to reset contrast)\n' \
         'esc: reset keyboard mode\n'
PROCESSES = 5


class ConsoleBoxStream(object):
    """
    A class that writes to a console box as a stream.
    """
    def __init__(self, console_box):
        self.console_box = console_box

    def write(self, content):
        """

        :param content: content to write to the console stream.
        :return:
        """
        self.console_box.append_text(content)

    def flush(self):
        pass


class ValidateGui(ipg.EnhancedCanvasView):

    def __init__(self, logger, window, bindings=None):
        """

        :param logger: a logger object to send messages to
        :type logger: logging.Logger
        :param window: The main window of the application
        :param bindings: Any bindings previously set on this window.
        """
        super(ValidateGui, self).__init__(logger=logger, bindings=bindings)

        self.console_box = Widgets.TextArea(editable=False)

        self.downloader = downloader.Downloader()
        self.pool = Pool(processes=PROCESSES)
        self.lock = Lock()
        self.image_list = {}
        self.astro_images = {}

        self.logger = logger
        console_handler = logging.StreamHandler(stream=ConsoleBoxStream(self.console_box))
        self.logger.addHandler(console_handler)
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
        self.pixel = None
        self.storage_list = None
        self.override = None
        self.qrun_id = None
        self.length_check = False

        # GUI elements
        self.pixel_base = 1.0
        self.readout = Widgets.Label("")
        self.header_box = Widgets.TextArea(editable=False)
        self.accept = Widgets.Button("Accept")
        self.reject = Widgets.Button("Reject")
        self.next_set = Widgets.Button("Next Set >")
        self.previous_set = Widgets.Button("< Previous Set")
        self.load_json = Widgets.Button("Load")

        self.legend = Widgets.TextArea(wrap=True)
        self.legend.set_text(LEGEND)
        self.build_gui(self.top)
        self.comparison_images = {}
        self.null_observation = {}
        self.next_image = None

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

        candidate_set = Widgets.TextEntry()
        candidate_set.add_callback('activated', lambda x: self.set_pixel(event=x))
        candidate_set.set_length(6)

        candidate_override = Widgets.TextEntry()
        candidate_override.add_callback('activated', lambda x: self.override_set(event=x))
        candidate_override.set_length(10)

        catalog = Widgets.TextEntrySet(text='17AQ10')
        catalog.add_callback('activated', lambda x: self.set_qrun_id(x))
        catalog.set_length(5)

        self.accept.add_callback('activated', lambda x: self.accept_reject())
        self.reject.add_callback('activated', lambda x: self.accept_reject(rejected=True))
        self.load_json.add_callback('activated', lambda x: self.load_candidates())
        self.next_set.add_callback('activated', lambda x: self.next())
        self.previous_set.add_callback('activated', lambda x: self.previous())

        # quit_button = Widgets.Button("Quit")
        # quit_button.add_callback('activated', lambda x: self.exit())

        reload_button = Widgets.Button("Reload")
        reload_button.add_callback('activated', lambda x: self.reload_candidates())

        # accept/reject/next buttons
        buttons_hbox = Widgets.HBox()
        buttons_hbox.add_widget(self.previous_set)
        buttons_hbox.add_widget(self.accept)
        buttons_hbox.add_widget(self.reject)
        buttons_hbox.add_widget(self.next_set)
        buttons_hbox.add_widget(self.load_json)
        self.load_json.set_enabled(False)
        buttons_hbox.add_widget(reload_button)
        buttons_hbox.set_spacing(10)

        # catalog directory text box
        catalog_box = Widgets.HBox()
        catalog_label = Widgets.Label(text="Set QRUNID:", style='color:red')
        catalog_box.add_widget(catalog_label)
        catalog_box.add_widget(catalog)
        catalog_box.set_margins(15, 0, 10, 0)  # top, right, bottom, left

        candidates_hbox = Widgets.HBox()
        candidate_label = Widgets.Label(text="(Optional) Enter candidate set: ")
        candidates_hbox.add_widget(candidate_label)
        candidates_hbox.add_widget(candidate_set)
        candidates_hbox.set_margins(15, 0, 15, 0)  # top, right, bottom, left

        override_hbox = Widgets.HBox()
        override_label = Widgets.Label(text="(Optional) Override provisional name for viewing: ")
        override_hbox.add_widget(override_label)
        override_hbox.add_widget(candidate_override)

        # button and text entry vbox
        buttons_vbox = Widgets.VBox()
        buttons_vbox.add_widget(buttons_hbox)
        buttons_vbox.add_widget(catalog_box)
        buttons_vbox.add_widget(candidates_hbox)
        buttons_vbox.add_widget(override_hbox)

        viewer_vbox.add_widget(buttons_vbox)  # add buttons below the viewer

        viewer_header_hbox = Widgets.HBox()  # box containing the viewer/buttons and rightmost text area
        viewer_header_hbox.add_widget(viewer_vbox)
        viewer_header_hbox.add_widget(Widgets.Label(''))
        hbox = Widgets.HBox()
        hbox.add_widget(self.header_box)
        hbox.add_widget(self.legend)
        viewer_header_hbox.add_widget(hbox)

        full_vbox = Widgets.VBox()  # vbox container for all elements
        full_vbox.add_widget(viewer_header_hbox)

        full_vbox.add_widget(self.console_box)
        self.console_box.set_text('Logging output:\n')
        self.header_box.set_text("Header:")
        container.set_widget(full_vbox)

    def next(self):
        """
        Load the next set of images into the viewer
        """
        if self.candidates is not None:
            # noinspection PyBroadException
            try:
                self.next_set.set_enabled(False)
                self.previous_set.set_enabled(False)
                self.accept.set_enabled(False)
                self.reject.set_enabled(False)
                self.obs_number = 0
                self.candidate = self.candidates.next()

                # Checks if candidate has already been examined with a file written on VOSpace.
                # length_check is necessary because it means the sub directory exists, if it doesn't an error will be
                # thrown when looking in the directory list.
                if self.length_check and self.candidate[0].provisional_name + '.ast' in storage.listdir(
                        os.path.join(os.path.dirname(storage.DBIMAGES), storage.CATALOG, self.qrun_id,
                                     self.candidates.catalog.catalog.dataset_name), force=True):

                    if self.override == self.candidate[0].provisional_name:
                        self.console_box.append_text("Candidate {} being overridden for viewing.\n"
                                                     .format(self.candidate[0].provisional_name))

                    else:
                        self.console_box.append_text("Candidate {} has been investigated.\n"
                                                     .format(self.candidate[0].provisional_name))
                        self.next()
                        return

                self.load()
                self.next_set.set_enabled(True)
                self.previous_set.set_enabled(True)
                self.accept.set_enabled(True)
                self.reject.set_enabled(True)
            except Exception as ex:
                self.console_box.append_text('Loading next candidate set failed.\n')
                if isinstance(ex, StopIteration):
                    self.console_box.append_text('StopIteration error: End of candidate set.\n'
                                                 'Hit "Load" button to move onto the next set.\n')
                    self.previous_set.set_enabled(True)
                    self.load_json.set_enabled(True)

    def previous(self):
        """
        Load the previous set of images into the viewer
        """
        if self.candidates is not None:
            self.next_set.set_enabled(False)
            self.previous_set.set_enabled(False)
            self.accept.set_enabled(False)
            self.reject.set_enabled(False)
            self.obs_number = 0
            self.candidate = self.candidates.previous()
            if self.candidate is not None:
                self.load()
            self.previous_set.set_enabled(True)
            self.next_set.set_enabled(True)
            self.accept.set_enabled(True)
            self.reject.set_enabled(True)

    def accept_reject(self, rejected=False):
        """
        Accept or reject current observation depending on button press. Write to file and load next set into the viewer

        :param rejected: whether the candidate set has been accepted or rejected
        """
        self.console_box.append_text("Rejected.\n") if rejected else self.console_box.append_text("Accepted.\n")
        if self.candidates is not None:
            self.write_record(rejected=rejected)
            self.next()

    def set_qrun_id(self, qrun_id):
        """
        :param qrun_id: QRUNID in a header file
        """
        self.qrun_id = qrun_id.text
        self.storage_list = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                                         storage.CATALOG,
                                                         self.qrun_id), force=True)
        self.load_json.set_enabled(True)
        self.console_box.append_text("QRUNID set to {}. \n".format(self.qrun_id))

    def lookup(self):
        """
        Looks up which pixel value for candidate files in VOSpace.

        :return pixel value for the candidate files; 0 if no candidate files have been found
        """
        count = 0
        self.length_check = False
        for filename in self.storage_list:

            # ex: filename = HPX_00887_RA_203.6_DEC_+58.9_bk.json,
            #     filename[:-len(storage.MOVING_TARGET_VERSION)] = HPX_00887_RA_203.6_DEC_+58.9
            # sub_directory will be the directory where a candidate's .ast files are written
            sub_directory = filename[:-len(storage.MOVING_TARGET_VERSION + '.json')]
            count += 1

            # if the file extension is in the filename, then it is a file containing candidate information
            if storage.MOVING_TARGET_VERSION in filename:
                x = re.match('(?P<hpx>HPX_)(?P<pixel>\d{5})(?P<leftover>_.*)', filename)

                if self.pixel is not None and int(x.group('pixel')) < self.pixel:
                    continue  # skipping over json files until the specified catalog has been reached

                # if the sub directory exists, we will have to check that all the candidates have been investigated
                elif sub_directory in self.storage_list:
                    self.length_check = True

                # TODO: go back to server for storage_list in case two people are actively writing from unique servers
                # cutting down the storage list for further iterating
                self.storage_list = self.storage_list[count:]
                return int(x.group('pixel'))
        return 0

    def set_pixel(self, event):
        """
        Sets the pixel for the current Candidate set.

        :param event: Pixel value
        """
        if hasattr(event, 'text'):
            self.pixel = int(event.text)
            self.console_box.append_text("Set pixel as {}\n".format(self.pixel))
            if self.qrun_id is not None:
                self.load_json.set_enabled(True)

    def override_set(self, event):
        """
        Look at the cutout even if it has already been investigated. Primarily used for double checking
         accepted candidates.
        """
        if hasattr(event, 'text'):
            self.override = str(event.text)
            self.console_box.append_text("Will override {}.\n".format(self.override))

    def load_candidates(self, pixel=None):
        """
        Initial candidates loaded into the viewer. Starts up a thread pool to download images simultaneously.

        :param pixel: Catalogue number containing dataset
        """
        self.load_json.set_enabled(False)

        if pixel is None:
            self.pixel = self.lookup()

        if self.pixel == 0:  # recursive base case (when there are no more open candidate sets in the VOSpace directory)
            self.console_box.append_text("No more candidate sets for this QRUNID.\n")
            raise StopIteration

        self.console_box.append_text("Accepted candidate entry: {}\n".format(self.pixel))

        try:
            self.candidates = candidate.CandidateSet(self.pixel, catalog_dir=self.qrun_id)

            if self.length_check:
                sub_directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                                             storage.CATALOG,
                                                             self.qrun_id,
                                                             self.candidates.catalog.catalog.dataset_name),
                                                force=True)
                if self.override is not None:
                    x = self.override+'.ast'
                    if x in sub_directory:
                        self.console_box.append_text("Overriding {}.\n".format(x))
                else:
                    count = 0
                    # counting the total amount of candidates that are in self.candidates
                    for _ in self.candidates:
                        count += 1

                    # re-set self.candidates since the for loop removes all its candidates in a dequeuing fashion
                    self.candidates = candidate.CandidateSet(self.pixel, catalog_dir=self.qrun_id)
                    # the amount of files in the accompanying subdirectory for the .json candidate file
                    directory_length = len(sub_directory)

                    if count == directory_length:
                        self.console_box.append_text("Candidate set {} fully examined.\n".format(self.pixel))
                        self.load_candidates()
                        return

                    elif count > directory_length:
                        self.console_box.append_text("Candidate set {} not fully examined.\n".format(self.pixel))

                    else:
                        logging.error("Value error: count {} or directory_length {} is out of range."
                                      .format(count, directory_length))
                        raise ValueError

        except Exception as ex:
            self.console_box.append_text("Failed to load candidates: {} \n".format(str(ex)))
            if isinstance(ex, StopIteration):
                self.console_box.append_text('StopIteration error. Candidate set might be empty.\n')
                self.load_candidates()  # recursive call to find next candidate which needs validation
                return  # prevents further execution of this method, handle better by shortening this method?
            else:
                raise ex

        self.logger.warning("Launching image prefetching. Please be patient.\n")

        with self.lock:
            for obs_records in self.candidates:
                previous_record = None
                for obs_record in obs_records:
                    assert isinstance(obs_record, ObsRecord)
                    key = self.downloader.image_key(obs_record)
                    if key not in self.image_list:
                        self.image_list[key] = self.pool.apply_async(self.downloader.get, (obs_record,))

                    # Check if we should load a comparison for the previous image.
                    if previous_record is not None:
                        offset = obs_record.coordinate.separation(previous_record.coordinate)
                        if offset > storage.CUTOUT_RADIUS:
                            # Insert a blank image in the list
                            previous_key = self.downloader.image_key(previous_record)
                            comparison = storage.get_comparison_image(previous_record.coordinate,
                                                                      previous_record.date.mjd)
                            frame = "{}{}".format(comparison[0]['observationID'], 'p00')
                            comparison_obs_record = ObsRecord(null_observation=True,
                                                              provisional_name=previous_record.provisional_name,
                                                              date=Time(comparison[0]['mjdate'], format='mjd',
                                                                        precision=5).mpc,
                                                              ra=previous_record.coordinate.ra.degree,
                                                              dec=previous_record.coordinate.dec.degree,
                                                              frame=frame,
                                                              comment=previous_key)
                            key = self.downloader.image_key(comparison_obs_record)
                            self.null_observation[key] = comparison_obs_record
                            self.comparison_images[previous_key] = key
                            if key not in self.image_list:
                                self.image_list[key] = self.pool.apply_async(self.downloader.get,
                                                                             (comparison_obs_record,))

                    previous_record = obs_record

        self.candidates = candidate.CandidateSet(self.pixel, catalog_dir=self.qrun_id)
        self.candidate = None  # reset on candidate to clear it of any leftover from previous sets
        self.load()

    def reload_candidates(self):
        """
        Performs a hard reload on all images for the case of loading errors.
        Closes current worker pool and reopens a new one.
        """
        if self.pixel is not None:
            self.console_box.append_text('Reloading all candidates...\n')
            self.pool.terminate()
            self.pool = Pool(processes=PROCESSES)
            self.next_set.set_enabled(True)
            self.previous_set.set_enabled(True)
            self.accept.set_enabled(True)
            self.reject.set_enabled(True)
            self.load_candidates(self.pixel)
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
            self.console_box.append_text("No candidates loaded.\n")
            return

        # loads first candidate
        if self.candidate is None:
            self.next()
        key = self.key
        while True:
            # noinspection PyBroadException
            try:
                if key not in self.astro_images:
                    image = AstroImage.AstroImage(logger=self.logger)
                    image.load_hdu(self.loaded_hdu)
                    self.astro_images[key] = image

                self.set_image(self.astro_images[key])
                break

            except Exception as ex:
                self.console_box.append_text(str(ex) + '\n')
                self.console_box.append_text("Skipping candidate {} due to load failure\n"
                                             .format(self.candidate[0].provisional_name))
                self.logger.debug(str(ex))
                self.logger.debug("Skipping candidate {} due to load failure\n".format(self.obs_number))
                self.next()

        if self.zoom is not None:
            self.zoom_to(self.zoom)

        self._rotate()

        if self.center is not None:
            self._align()

        # the image cutout is considered the first object on the canvas, this deletes everything over top of it
        self.canvas.delete_objects(self.canvas.get_objects()[1:])

        if key not in self.null_observation:
            self._mark_aperture()

        self.header_box.set_text("Header:\n" + self.info)
        self.console_box.append_text("Loaded: {}\n".format(self.candidate[self.obs_number].comment.frame))

    def _mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object detected.
        """
        ra = self.candidate[self.obs_number].coordinate.ra
        dec = self.candidate[self.obs_number].coordinate.dec
        x, y = WCS(self.header).all_world2pix(ra, dec, 0)
        self.canvas.add(self.circle(x, y, radius=10, color='red'))

    def write_record(self, rejected=False):
        """
        Writing observation lines to a new file.

        :param rejected: Whether or not the candidate set contains a valid celestial object
        :type rejected: bool
        """
        try:
            catalog_dir = os.path.join(storage.CATALOG,
                                       self.header['QRUNID'],
                                       self.candidates.catalog.catalog.dataset_name)

            art = storage.ASTRecord(self.candidate[0].provisional_name,
                                    version='',
                                    catalog_dir=catalog_dir)

            with open(art.filename, 'w+') as fobj:
                for ob in self.candidate:
                    if rejected:
                        ob.null_observation = True
                    fobj.write(ob.to_string() + '\n')

            self.logger.info("Queuing job to write file to VOSpace.")
            self.pool.apply_async(self.downloader.put, (art,))
            msg = "Done Queuing {} for VOSpace write {}".format(self.candidate[0].provisional_name + ".ast",
                                                                art.uri)
            logging.info(msg)
            self.console_box.append_text(msg+"\n")

        except IOError as ex:
            self.console_box.append_text("Unable to write to file.")
            self.console_box.append_text(str(ex) + '\n')
            raise ex

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
        self.logger.debug("Got key: {} from canvas: {} with opn: {} from viewer: {}".format(canvas,
                                                                                            keyname,
                                                                                            opn,
                                                                                            viewer))
        # Only step back if we aren't looking at a comparison images (as determined by the next_image keyword)
        if keyname == 'f':
            if self.next_image is not None:
                self.next_image = None
            else:
                self.obs_number -= 1
                key = self.downloader.image_key(self.candidate[self.obs_number])
                if key in self.comparison_images:
                    self.next_image = self.comparison_images[key]

        # Only step forward if this images doesn't have comparison image in the comparison image list.
        elif keyname == 'g':
            key = self.downloader.image_key(self.candidate[self.obs_number])
            if key in self.comparison_images and self.next_image is None:
                self.next_image = self.comparison_images[key]
            else:
                self.next_image = None
                self.obs_number += 1

        self.zoom = self.get_zoom()
        self.obs_number %= len(self.candidate)
        self._load()

    def clear_candidate_images(self):
        """
        Clear all the images associated with a candidate.
        :return:
        """
        if self.candidate is None:
            return
        for obs_record in self.candidate:
            key = self.downloader.image_key(obs_record)
            del(self.image_list[key])
            if key in self.astro_images:
                del(self.astro_images[key])
            if key in self.comparison_images:
                comp_key = self.comparison_images[key]
                del(self.image_list[comp_key])
                if comp_key in self.astro_images:
                    del(self.astro_images[comp_key])

    @property
    def center(self):
        """
        Returns the center of the image in ra/dec coordinates
        """
        if self._center is not None:
            return self._center

    @property
    def key(self):
        if self.next_image is not None:
            key = self.next_image
        else:
            key = self.downloader.image_key(self.candidate[self.obs_number])
        return key

    @property
    def loaded_hdu(self):
        """
        Return current HDU
        """
        key = self.key
        with self.lock:
            hdu = (isinstance(self.image_list[key], ApplyResult) and self.image_list[key].get()
                   or self.image_list[key])
            self.image_list[key] = hdu

        load_hdu = max_size = None
        for hdu in self.image_list[key]:
            if hdu.header['NAXIS'] == 0:
                continue
            size = hdu.header['NAXIS1'] * hdu.header['NAXIS2']
            if max_size is None or size > max_size:
                max_size = size
                load_hdu = hdu
        return load_hdu

    @property
    def header(self):
        """
        Return current HDU's header
        """
        return self.astro_images[self.key].get_header()

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

    ginga_logger = log.get_logger("ginga", options=params)

    if params.use_opencv:
        from ginga import trcalc
        try:
            trcalc.use('opencv')
        except Exception as ex:
            ginga_logger.warning("Error using OpenCL: {}".format(ex))

    if params.use_opencl:
        from ginga import trcalc
        try:
            trcalc.use('opencl')
        except Exception as ex:
            ginga_logger.warning("Error using OpenCL: {}".format(ex))

    app = Widgets.Application(logger=ginga_logger, host=params.host, port=params.port)

    #  create top level window
    window = app.make_window("Validate", wid='Validate')

    # our own viewer object, customized with methods (see above)
    ValidateGui(logging.getLogger('daomop'), window)

    try:
        app.start()

    except KeyboardInterrupt:
        ginga_logger.info("Terminating viewer...")
        window.close()
