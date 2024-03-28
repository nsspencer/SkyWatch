import datetime
import time
import unittest

import astropy.units as u
import numpy as np
from astropy.time import Time, TimeDelta
from scipy.spatial.transform import Rotation

from skywatch.access import Access, TimeInterval
from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.attitude import LVLH
from skywatch.look_angles import LocalTangentENU, XForwardZNadir
from skywatch.skypath import SkyPath
from skywatch.utils.coverage import GeoFilter, calculate_coverage

from .tests import get_ephem_as_skypath


class SmokeTests(unittest.TestCase):

    def test_access(self):
        """
        Calculate very specific access from a ground station to a satellite.
        """

        t_start = Time("2024-02-01T00:00:00")
        t_end = Time("2024-02-02T00:00:00")
        low_fidelity_times = np.linspace(t_start, t_end, 1440)
        high_fidelity_times = np.linspace(t_start, t_end, 86400)

        earth_pos = SkyPath.from_body(low_fidelity_times, "earth")
        sun_pos = SkyPath.from_body(low_fidelity_times, "sun")
        moon_pos = SkyPath.from_body(low_fidelity_times, "moon")

        ground_station_pos = SkyPath.from_geodetic(
            low_fidelity_times[0], 34.5 * u.deg, -77.0 * u.deg, 0 * u.m
        )

        sat_position = get_ephem_as_skypath()

        t0 = time.time()
        access_times = (
            Access(
                LineOfSight(ground_station_pos, sat_position, earth_pos),
                AzElRange(
                    ground_station_pos, sun_pos, min_el=0 * u.deg, max_el=5 * u.deg
                ),
                LineOfSight(sat_position, sun_pos, earth_pos),
                Temporal(
                    min_time=Time("2024-02-01T22:15:00"),
                    max_time=Time("2024-02-01T22:17:00"),
                ),
                LineOfSight(moon_pos, sun_pos, earth_pos),
            )
            .use_precise_endpoints(True)
            .set_precision(0.001 * u.s)
            .set_min_duration(60 * u.s)
            .calculate_at(high_fidelity_times)
        )
        t1 = time.time()

        print(f"Access calculation took: {t1-t0} seconds")
        print(access_times.total_duration)
        print(access_times)

    def test_eclipse_access(self):
        """
        calculate approximate lunar and solar eclipse times until 2025
        """

        t_start = Time("2024-01-01T00:00:00")
        t_end = Time("2025-01-01T00:00:00")
        dt = TimeDelta(datetime.timedelta(days=1), format="datetime")
        num_steps = int((t_end - t_start) / dt)
        times = np.linspace(t_start, t_end, num_steps)

        earth_pos = SkyPath.from_body(times, "earth")
        sun_pos = SkyPath.from_body(times, "sun")
        moon_pos = SkyPath.from_body(times, "moon")

        t0 = time.time()
        lunar_eclipse_times = (
            Access(
                LineOfSight(
                    sun_pos, moon_pos, earth_pos, when_obstructed=True, use_frame="icrs"
                )
            )
            .use_precise_endpoints(True)
            .set_precision(1 * u.s)
            .calculate_at(times)
        )
        t1 = time.time()

        print(f"Lunar eclipse access calculation took: {t1-t0} seconds")
        print(lunar_eclipse_times.total_duration)
        print(lunar_eclipse_times)

        t0 = time.time()
        solar_eclipse_times = (
            Access(
                LineOfSight(
                    earth_pos,
                    sun_pos,
                    moon_pos,
                    sma=1079.6 * 6 * u.km,  # use a body size ~= earth
                    smi=1079.6 * 6 * u.km,  # use a body size ~= earth
                    when_obstructed=True,
                    use_frame="icrs",
                )
            )
            .use_precise_endpoints(True)
            .set_precision(1 * u.s)
            .calculate_at(times)
        )
        t1 = time.time()

        print(f"Solar eclipse access calculation took: {t1-t0} seconds")
        print(solar_eclipse_times.total_duration)
        print(solar_eclipse_times)

    def test_coverage(self):
        t_start = Time("2024-02-01T00:00:00")
        t_end = Time("2024-02-02T00:00:00")
        times = np.linspace(t_start, t_end, 8640)

        sat_position = get_ephem_as_skypath()

        worldwide_coverage_result = calculate_coverage([sat_position], times, 300)
        print(
            f"Worldwide min revisit time: {worldwide_coverage_result.min_revisit_time}"
        )
        print(
            f"Worldwide max revisit time: {worldwide_coverage_result.max_revisit_time}"
        )

        at_equator_coverage = worldwide_coverage_result.filter(
            GeoFilter(min_latitude=-1 * u.deg, max_latitude=1 * u.deg)
        )
        print(f"Equatorial min revisit time: {at_equator_coverage.min_revisit_time}")
        print(f"Equatorial max revisit time: {at_equator_coverage.max_revisit_time}")


if __name__ == "__main__":
    unittest.main()
