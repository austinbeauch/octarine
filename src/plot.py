"""
Small command line application for using plot_viewer.py to plot/visualize data from generated .FITS files.
Also used to run st_deviation.py for generating magnitude standard deviation data files.

Usage:
1. from ~/octarine/src:
    $ python plot.py
     - generate sky coverage scatter plots using data from /fits_data/

Optional:
1. [--hist]
    - plot histograms using data from the directory /mag_std_data/

2. [--std]
    - generates standard deviation data files from master catalogs and store them in /mag_std_data/

3. [--mag]
    - generate flux limit coverage/overlaps data files and store them in /fits_data/
    - data files available to download from VOSpace under /cfis/solar_system/characterization
"""

import argparse
from daomop import storage
from plotting import plot_viewer, st_deviation, mag_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--hist", dest="hist", default=False, action="store_true",
                        help="Plot histograms")
    parser.add_argument("--std", dest="std", default=False, action="store_true",
                        help="Run st_deviation_plot.py to generate standard deviation data files.")
    parser.add_argument("--mag", dest="mag", default=False, action="store_true",
                        help="Run mag_data.py to generate coverage data files.")

    args = parser.parse_args()

    if args.std:
        st_deviation.main()
    elif args.mag:
        mag_data.main()
    else:
        plot_viewer.main(args)
