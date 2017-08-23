"""
Small command line application for using plot_viewer.py to plot/visualize data from generated .FITS files.
Also used to run st_deviation.py for generating magnitude standard deviation data files.

Usage:
1. from ~/octarine/src:
    $ python plot.py
     - generate sky coverage scatter plots using data from /fits_data/

Optional:
1. [--std]
    - generates standard deviation files from master catalogs and store them in /mag_std_data/

2. [--hist]
    - generates histograms using data from the directory /mag_std_data/
"""

import argparse
from daomop import storage
from plotting import plot_viewer, st_deviation


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist", dest="hist", default=False, action="store_true",
                        help="Plot histograms")
    parser.add_argument("--std", dest="std", default=False, action="store_true",
                        help="Run st_deviation_plot.py to generate standard deviation data files.")

    args = parser.parse_args()

    if args.std:
        st_deviation.main()
    else:
        plot_viewer.main(args)
