import logging
import storage
from multiprocessing import Lock


class Downloader(object):
    def __init__(self):
        self._downloaded_images = dict()
        self.lock = Lock()
        self.locks = {}

    @staticmethod
    def image_key(obs_record):
        """

        :param obs_record: the observation record to get the key for.
        :type obs_record: mp_ephem.ObsRecord
        :return: the key for the given observation record.
        """
        return obs_record.comment.frame + obs_record.provisional_name

    def get(self, obs_record):

        with self.lock:
            if self.image_key(obs_record) not in self.locks:
                self.locks[self.image_key(obs_record)] = Lock()

        with self.locks[self.image_key(obs_record)]:
            if self.image_key(obs_record) not in self._downloaded_images:
                self._downloaded_images[self.image_key(obs_record)] = self.get_hdu(obs_record)

        return self._downloaded_images[self.image_key(obs_record)]

    def get_hdu(self, obs_record):
        """
        Retrieve a fits image associated with a given obs_record and return the HDU off the assocaited cutout.

        :param obs_record: the ObsRecord for which the image is to be retrieved
        :type obs_record: mp_ephem.ObsRecord
        :return: fits.ImageHDU
        """
        logging.info("Retrieving {}".format(self.image_key(obs_record)))
        try:
            image = storage.FitsImage.from_frame(obs_record.comment.frame)
            hdu = image.ra_dec_cutout(obs_record.coordinate)[-1]
        except Exception as ex:
            logging.info(ex)
            raise ex

        logging.info("Got {}".format(self.image_key(obs_record)))
        return hdu
