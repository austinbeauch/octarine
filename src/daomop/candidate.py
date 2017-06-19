import sys

import mp_ephem
from mp_ephem import ObsRecord
import storage
from astropy.time import Time
import name


class ObservationSet(object):
    def __init__(self, provisional_name, record):
        self.provisional_name = provisional_name
        self.record = record
        self.current_observation = -1

    def __iter__(self):
        return self

    def next(self):
        self.current_observation += 1
        if not self.current_observation < len(self.record['mag']):
            raise StopIteration
        return ObsRecord(provisional_name=self.provisional_name,
                         discovery=True,
                         note1=None,
                         note2='C',
                         date=Time(self.record['mjd'][self.current_observation], format='mjd'),
                         ra=self.record['ra'][self.current_observation],
                         dec=self.record['dec'][self.current_observation],
                         mag=self.record['mag'][self.current_observation],
                         band=self.record['filterid'][self.current_observation],
                         observatory_code=568,
                         mag_err=self.record['magerr'][self.current_observation],
                         comment='cand',
                         xpos=None,
                         ypos=None,
                         frame=self.record['fitsname'][self.current_observation],
                         plate_uncertainty=None,
                         astrometric_level=4
                         )


class Target(object):
    def __init__(self, hpx, mjdate, record):
        self.hpx = hpx
        self.mjdate = mjdate
        self.record = record
        self.current_observation = -1

    def __iter__(self):
        return self

    @property
    def observation_sets(self):
        return self.record.keys()

    @property
    def provisional_name(self):
        return name.provisional(self.mjdate, self.hpx, self.current_observation)

    def next(self):
        self.current_observation += 1
        if not self.current_observation < len(self.observation_sets):
            raise StopIteration
        return ObservationSet(self.provisional_name,
                              self.record[self.observation_sets[self.current_observation]])


class Catalog(object):
    def __init__(self, pixel):
        self.catalog = storage.JSONCatalog(pixel)
        self.current_target = -1

    def __iter__(self):
        return self

    @property
    def mjdates(self):
        return self.catalog.json.keys()

    def next(self):
        self.current_target += 1
        if not self.current_target < len(self.mjdates):
            raise StopIteration
        mjdate = self.mjdates[self.current_target]
        return Target(self.catalog.pixel, mjdate, self.catalog.json[mjdate])


class CandidateSet(object):

    def __init__(self, pixel):
        self.catalog = Catalog(pixel)
        self.target = self.catalog.next()

    def __iter__(self):
        return self

    def next(self):
        try:
            observation_set = self.target.next()
            obs = []
            for ob in observation_set:
                obs.append(ob)
            return mp_ephem.BKOrbit(obs)
        except StopIteration:
            self.target = self.catalog.next()
            return self.next()


if __name__ == "__main__":
    print sys.argv
    for target in Catalog(int(sys.argv[1])):
        for obs_set in target:
            for observation in obs_set:
                print(observation)
            print("")
