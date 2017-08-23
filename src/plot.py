description = """
Small command line application for using plot_viewer.py to plot/visualize data from generated .FITS files.
Also used to run st_deviation.py and mag_data.py for generating magnitude standard deviation data files.
"""

import argparse
from daomop import storage
from plotting import plot_viewer, st_deviation, mag_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--cov", dest="cov", default=False, action="store_true",
                        help="Plot sky coverage per magnitude, overlaps, and stellar densities from data in "
                             "/fits_data/")
    parser.add_argument("--hist", dest="hist", default=False, action="store_true",
                        help="Plot histograms using data in /mag_std_data/")
    parser.add_argument("--std", dest="std", default=False, action="store_true",
                        help="Run st_deviation_plot.py to generate standard deviation data files.")
    parser.add_argument("--mag", dest="mag", default=False, action="store_true",
                        help="Run mag_data.py to generate coverage data files.")

    args = parser.parse_args()

    if args.std:
        st_deviation.main()
    elif args.mag:
        mag_data.main()
    elif args.cov:
        plot_viewer.main(args)
    else:
        parser.error('No action requested, add exactly one argument')

__doc__ = description
