"""
Build an input file to be run via calling condor on canfar.
"""
import storage
import sys


for healpix in storage.list_healpix():

    params = {"Arguments": healpix,
              "Log": "{}.log".format(healpix),
              "Output": "{}.out".format(healpix),
              "Error": "{}.err".format(healpix)
              }

    for param in params:
        sys.stdout.write("{} = {}\n".format(param, params[param]))
    sys.stdout.write("Qeueu\n\n")
