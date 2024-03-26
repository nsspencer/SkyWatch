import unittest

import astropy.units as u
from astropy.time import Time

from skywatch.look_angles import LocalTangentENU
from skywatch.platform import Platform
from skywatch.skypath import SkyPath


class Test1(unittest.TestCase):
    def test_ground_station_look_angles(self):
        ground_station = Platform(
            SkyPath.from_geodetic(Time("J2000"), 10 * u.deg, 10 * u.deg, 0 * u.m),
            look_angles_strategy=LocalTangentENU(),
            attitude=None,
        )

        ground_station2 = Platform(
            SkyPath.from_geodetic(Time("J2000"), 19 * u.deg, 1 * u.deg, 155000 * u.m),
            look_angles_strategy=LocalTangentENU(),
            attitude=None,
        )

        print(ground_station.look_angles(Time("J2000"), ground_station2))
        print(ground_station2.look_angles(Time("J2000"), ground_station))
