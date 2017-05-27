"""
Build an input file to be run via calling condor on canfar.
"""
import storage
import sys

command="stationary.sh"

with open(command, 'w') as fout:
   fout.write("""#!/bin/bash -i
export HOME=`cd ~ ; pwd`
source activate ossos
getCert

echo "Building Catalog ",$1

stattionary $1 --verbose --catalogs catalogs_20170527

""")

os.chmod(command, stat.S_IXGRP|stat.S_IRWXU|stat.S_IXOTH|stat.S_IROTH|stat.S_IRGRP)



print("""Universe   = vanilla
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
RunAsOwner = True
transfer_output_files = /dev/null""")



for healpix in storage.list_healpix():

    params = {"Arguments": healpix,
              "Log": "{}.log".format(healpix),
              "Output": "{}.out".format(healpix),
              "Error": "{}.err".format(healpix)
              }

    for param in params:
        sys.stdout.write("{} = {}\n".format(param, params[param]))
    sys.stdout.write("Qeueu\n\n")
