from daomop import viewer, candidate
from daomop import storage
import logging
import sys


def run(pixel):
    storage.DBIMAGES = 'vos:cfis/solar_system/dbimages'
    downloader = storage.Downloader()
    with viewer.WebServerFactory() as web_server:
        gui = web_server.get_viewer("validate")
        viewer.ImageViewer(gui, downloader)

        v1 = viewer.ValidateGui(web_server.get_viewer('Validate'), candidates, downloader)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run(int(sys.argv[1]))

