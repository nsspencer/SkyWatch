import unittest

import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude import LVLH
from skywatch.look_angles import LocalTangentENU, XForwardZNadir
from skywatch.skypath import SkyPath

from .tests import get_ephem_as_skypath


class Test1(unittest.TestCase):
    def test_ground_station_to_ground_station(self):
        gs1 = SkyPath.from_geodetic(Time("J2000"), 10 * u.deg, 10 * u.deg, 0 * u.m)
        gs2 = SkyPath.from_geodetic(Time("J2000"), 19 * u.deg, 1 * u.deg, 155000 * u.m)

        gs1_look_angles = LocalTangentENU(gs1)
        gs2_look_angles = LocalTangentENU(gs2)
        gs1_to_gs2 = gs1_look_angles.get_look_angles(gs2, Time("J2000"))
        gs2_to_gs1 = gs2_look_angles.get_look_angles(gs1, Time("J2000"))
        print(gs1_to_gs2)
        print(gs2_to_gs1)

    def test_ground_station_to_satellite(self):
        sat = get_ephem_as_skypath()
        gs1 = SkyPath.from_geodetic(Time("J2000"), 10 * u.deg, 10 * u.deg, 0 * u.m)
        gs2 = SkyPath.from_geodetic(sat.obstime[0], 10 * u.deg, 10 * u.deg, 0 * u.m)

        gs_look_angles = LocalTangentENU(gs1)
        gs2sat = gs_look_angles.get_look_angles(sat, sat.obstime)

        gs2_look_angles = LocalTangentENU(gs2)
        gs2sat2 = gs2_look_angles.get_look_angles(sat, sat.obstime)

        self.assertTrue(np.allclose(gs2sat.azimuth, gs2sat2.azimuth))
        self.assertTrue(np.allclose(gs2sat.elevation, gs2sat2.elevation))
        self.assertTrue(np.allclose(gs2sat.range, gs2sat2.range))

    def test_satellite_to_ground_station(self):
        sat = get_ephem_as_skypath()
        ground_station = SkyPath.from_geodetic(
            sat.obstime, 34.5 * u.deg, -77.0 * u.deg, 0 * u.m
        )

        sat_look_angles_fn = XForwardZNadir(sat, LVLH(sat, "itrs"), "itrs")
        sat_look_angles = sat_look_angles_fn.get_look_angles(
            ground_station, sat.obstime
        )

        gs_look_angles_fn_gcrs = XForwardZNadir(sat, LVLH(sat, "gcrs"), "gcrs")
        sat_look_angles_gcrs = gs_look_angles_fn_gcrs.get_look_angles(
            ground_station, sat.obstime
        )

        _diff_azimuths = np.abs(sat_look_angles.azimuth - sat_look_angles_gcrs.azimuth)
        self.assertTrue(
            np.max(_diff_azimuths[np.where(_diff_azimuths.value < 358.1)[0]])
            < 2.0 * u.deg
        )

        self.assertTrue(
            np.allclose(sat_look_angles.elevation, sat_look_angles_gcrs.elevation)
        )
        self.assertTrue(np.allclose(sat_look_angles.range, sat_look_angles_gcrs.range))

        sat_look_angles_offset_fn = XForwardZNadir(
            sat, LVLH(sat, "itrs", Rotation.from_euler("z", -45, degrees=True))
        )
        sat_2_gs_offset = sat_look_angles_offset_fn.get_look_angles(
            ground_station, sat.obstime
        )
        pass
