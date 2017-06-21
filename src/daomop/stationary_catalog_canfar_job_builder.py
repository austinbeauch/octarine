"""
Build an input file to be run via calling condor on canfar.
"""
import storage
import sys
import os
import stat
from astropy.time import Time

_COMMAND_FILENAME = "stationary.sh"

CFHT_QRUNS = {
    '17AQ01': (Time('2017-02-01 22:00:00'), Time('2017-02-02 22:00:00')),
    '17AQ03': (Time('2017-02-17 22:00:00'), Time('2017-03-02 22:00:00')),
    '17AQ06': (Time('2017-03-20 22:00:00'), Time('2017-04-03 22:00:00')),
    '17AQ08': (Time('2017-04-18 22:00:00'), Time('2017-05-03 22:00:00')),
    '17AQ10': (Time('2017-05-16 22:00:00'), Time('2017-05-30 22:00:00')),
    '17AQ13': (Time('2017-06-15 22:00:00'), Time('2017-06-27 22:00:00')),
    '17AQ16': (Time('2017-07-14 22:00:00'), Time('2017-07-31 22:00:00'))
}

def _create_shell_script(filename):
    """
    Create the shell script that is the processing step.

    :param filename:  target_name of file to contain the processing shell script.
    :return:
    """
    with open(filename, 'w') as fout:
        fout.write("""#!/bin/bash -i
export HOME=`cd ~ ; pwd`
source activate ossos
getCert

echo "Building Catalog ",$1

stattionary $1 --verbose --catalogs catalogs_20170527

""")

    os.chmod(filename, stat.S_IXGRP | stat.S_IRWXU | stat.S_IXOTH | stat.S_IROTH | stat.S_IRGRP)


def _create_job_file_header(command_filename):
    """
    Create the header of the Condor job.in file.
    :param command_filename: Name of the executable script that will be run.
    :return:
    """
    return """Universe   = vanilla
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
RunAsOwner = True
transfer_output_files = /dev/null

Executable = {}

""".format(command_filename)


def main(qrun):
    """Build the stationary catalog builder job submission script."""

    _create_shell_script(_COMMAND_FILENAME)
    with open('job.in', 'w') as job:
        job.write(_create_job_file_header(_COMMAND_FILENAME))

        for healpix in storage.list_healpix():

            params = {"Arguments": "{} {} {}".format(healpix, CFHT_QRUNS[qrun][0].mjd, CFHT_QRUNS[qrun][1].mjd),
                      "Log": "{}.log".format(healpix),
                      "Output": "{}.out".format(healpix),
                      "Error": "{}.err".format(healpix)
                      }

            for param in params:
                job.write("{} = {}\n".format(param, params[param]))
            job.write("Qeueu\n\n")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
