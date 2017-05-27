url="""http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap/sync?LANG=ADQL&REQUEST=doQuery&QUERY=SELECT%20Observation.observationURI%20AS%20%22Preview%22%2C%20Observation.collection%20AS%20%22Collection%22%2C%20Observation.sequenceNumber%20AS%20%22Sequence%20Number%22%2C%20Plane.productID%20AS%20%22Product%20ID%22%2C%20COORD1(CENTROID(Plane.position_bounds))%20AS%20%22RA%20(J2000.0)%22%2C%20COORD2(CENTROID(Plane.position_bounds))%20AS%20%22Dec.%20(J2000.0)%22%2C%20Observation.target_name%20AS%20%22Target%20Name%22%2C%20Plane.time_bounds_lower%20AS%20%22Start%20Date%22%2C%20Plane.time_exposure%20AS%20%22Int.%20Time%22%2C%20Observation.instrument_name%20AS%20%22Instrument%22%2C%20Plane.energy_bandpassName%20AS%20%22Filter%22%2C%20Plane.calibrationLevel%20AS%20%22Cal.%20Lev.%22%2C%20Observation.type%20AS%20%22Obs.%20Type%22%2C%20Observation.proposal_id%20AS%20%22Proposal%20ID%22%2C%20Observation.proposal_pi%20AS%20%22P.I.%20Name%22%2C%20Plane.dataRelease%20AS%20%22Data%20Release%22%2C%20Observation.observationID%20AS%20%22Obs.%20ID%22%2C%20Plane.energy_bounds_lower%20AS%20%22Min.%20Wavelength%22%2C%20Plane.energy_bounds_upper%20AS%20%22Max.%20Wavelength%22%2C%20AREA(Plane.position_bounds)%20AS%20%22Field%20of%20View%22%2C%20Plane.position_bounds%20AS%20%22Polygon%22%2C%20Plane.position_sampleSize%20AS%20%22Pixel%20Scale%22%2C%20Plane.energy_resolvingPower%20AS%20%22Resolving%20Power%22%2C%20Plane.time_bounds_upper%20AS%20%22End%20Date%22%2C%20Plane.dataProductType%20AS%20%22Data%20Type%22%2C%20Observation.target_moving%20AS%20%22Moving%20Target%22%2C%20Plane.provenance_name%20AS%20%22Provenance%20Name%22%2C%20Plane.provenance_keywords%20AS%20%22Provenance%20Keywords%22%2C%20Observation.intent%20AS%20%22Intent%22%2C%20Observation.target_type%20AS%20%22Target%20Type%22%2C%20Observation.target_standard%20AS%20%22Target%20Standard%22%2C%20Plane.metaRelease%20AS%20%22Meta%20Release%22%2C%20Observation.algorithm_name%20AS%20%22Algorithm%20Name%22%2C%20Observation.proposal_title%20AS%20%22Proposal%20Title%22%2C%20Observation.proposal_keywords%20AS%20%22Proposal%20Keywords%22%2C%20Plane.position_resolution%20AS%20%22IQ%22%2C%20Observation.instrument_keywords%20AS%20%22Instrument%20Keywords%22%2C%20Plane.energy_transition_species%20AS%20%22Molecule%22%2C%20Plane.energy_transition_transition%20AS%20%22Transition%22%2C%20Observation.proposal_project%20AS%20%22Proposal%20Project%22%2C%20Plane.energy_emBand%20AS%20%22Band%22%2C%20Plane.provenance_reference%20AS%20%22Prov.%20Reference%22%2C%20Plane.provenance_version%20AS%20%22Prov.%20Version%22%2C%20Plane.provenance_project%20AS%20%22Prov.%20Project%22%2C%20Plane.provenance_producer%20AS%20%22Prov.%20Producer%22%2C%20Plane.provenance_runID%20AS%20%22Prov.%20Run%20ID%22%2C%20Plane.provenance_lastExecuted%20AS%20%22Prov.%20Last%20Executed%22%2C%20Plane.provenance_inputs%20AS%20%22Prov.%20Inputs%22%2C%20Plane.energy_restwav%20AS%20%22Rest-frame%20Energy%22%2C%20Observation.requirements_flag%20AS%20%22Quality%22%2C%20Plane.planeID%20AS%20%22planeID%22%2C%20isDownloadable(Plane.planeURI)%20AS%20%22DOWNLOADABLE%22%2C%20Plane.planeURI%20AS%20%22CAOM%20Plane%20URI%22%20FROM%20caom2.Plane%20AS%20Plane%20JOIN%20caom2.Observation%20AS%20Observation%20ON%20Plane.obsID%20%3D%20Observation.obsID%20WHERE%20%20(%20Plane.calibrationLevel%20%3D%20%271%27%20AND%20Plane.energy_bandpassName%20IN%20(%20%27r.MP9602%27%2C%27r.MP9601%27%20)%20AND%20Observation.instrument_name%20%3D%20%27MegaPrime%27%20AND%20Observation.collection%20%3D%20%27CFHT%27%20AND%20lower(Observation.proposal_title)%20LIKE%20%27%25cfis%25%27%20AND%20%20(%20Plane.quality_flag%20IS%20NULL%20OR%20Plane.quality_flag%20!%3D%20%27junk%27%20)%20)&FORMAT=csv"""

import requests
import os
import stat
from astropy.table import Table
from cStringIO import StringIO
from vos import Client
import sys

Client().copy('vos:sgwyn/tkBAD', 'tkBAD')

bad_expnum_list = []
for lines in open('tkBAD').readlines():
   try:
      bad_expnum_list.append(lines.split()[0]) 
   except:
      sys.stderr.write(lines)


command="myCommand.sh"

with open(command, 'w') as fout:
   fout.write("""#!/bin/bash -i
export HOME=`cd ~ ; pwd`
source activate ossos
getCert

echo "Processing ",$1

populate $1 --verbose
build_cat $1 --verbose

""")

os.chmod(command, stat.S_IXGRP|stat.S_IRWXU|stat.S_IXOTH|stat.S_IROTH|stat.S_IRGRP)


print("""Universe   = vanilla
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
RunAsOwner = True
transfer_output_files = /dev/null

Executable = {}

""".format(command))

for row in Table.read(StringIO(requests.get(url).content), format='csv'):
    if row['Sequence Number'] in bad_expnum_list:
       continue
    print("Arguements = {}".format(row['Sequence Number']))
    print("Output = {}.out".format(row['Sequence Number']))
    print("Log = {}.log".format(row['Sequence Number']))
    print("Error = {}.err".format(row['Sequence Number']))
    print("Queue")
    print("")
