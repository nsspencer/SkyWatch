import datetime
import time
import unittest

import astropy.units as u
import numpy as np
from astropy.coordinates import CartesianRepresentation, EarthLocation
from astropy.time import Time, TimeDelta

from skywatch.access import Access
from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.skypath import SkyPath
from skywatch.utils.coverage import GeoFilter, calculate_coverage
from skywatch.utils.funcs import lat_lon_to_xyz, FibonacciPointGenerator

from .utils import get_ephem_as_skypath


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

        ground_station_pos = SkyPath.from_geodetic(low_fidelity_times[0], 34.5 * u.deg, -77.0 * u.deg, 0 * u.m)

        sat_position = get_ephem_as_skypath()

        t0 = time.time()
        access_times = (
            Access(
                LineOfSight(ground_station_pos, sat_position, earth_pos),
                AzElRange(ground_station_pos, sun_pos, min_el=0 * u.deg, max_el=5 * u.deg),
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
            .only_check_failed_constraints(True)
            .calculate_at(high_fidelity_times)
        )
        t1 = time.time()

        print(f"Access calculation took: {t1-t0} seconds")
        print(access_times.total_duration)
        print(access_times)

    def test_eclipse_access(self):
        """
        calculate approximate lunar and solar eclipse (total & partial) times until 2025
        """

        t_start = Time("2024-01-01T00:00:00")
        t_end = Time("2025-01-01T00:00:00")
        dt = TimeDelta(datetime.timedelta(hours=1), format="datetime")
        num_steps = int((t_end - t_start) / dt)
        times = np.linspace(t_start, t_end, num_steps)

        # define a sample of position around the surface of the moon
        moon_position_offsets = [lat_lon_to_xyz(i[0], i[1], radius=1079.6) for i in FibonacciPointGenerator(100)()]

        # get the moons position and add the offsets to the position to define
        # coordinates on the lunar surface
        moon_pos = SkyPath.from_body(times, "moon")
        moon_surface_samples = [
            SkyPath(
                CartesianRepresentation(*(moon_pos.cartesian.xyz.T + (moon_pos_offset * u.m)).T),
                frame="gcrs",
                obstime=times,
            )
            for moon_pos_offset in moon_position_offsets
        ]

        earth_pos = SkyPath.from_body(times, "earth")
        sun_pos = SkyPath.from_body(times, "sun")

        # calculate line of sight access from the sun to the moon when
        # obstructed by the earth
        t0 = time.time()
        lunar_access_interval = None
        for moon_sample in moon_surface_samples:
            lunar_eclipse_times = (
                Access(
                    LineOfSight(
                        sun_pos,
                        moon_sample,
                        earth_pos,
                        when_obstructed=True,
                        use_frame="gcrs",
                    )
                )
                .use_precise_endpoints(True)
                .set_precision(1 * u.s)
                .calculate_at(times)
            )

            # union the intervals together to get cumulative time
            if lunar_access_interval is None:
                lunar_access_interval = lunar_eclipse_times
            else:
                lunar_access_interval = lunar_access_interval | lunar_eclipse_times
        t1 = time.time()

        print(f"Lunar eclipse access calculation took: {t1-t0} seconds")
        print(f"There will be ~{lunar_access_interval.total_duration} seconds of lunar eclipse in 2024.")
        print(lunar_access_interval)
        self.assertTrue(len(lunar_access_interval) >= 1)

        # define a sample of points around the earth
        earth_positions = [
            SkyPath(EarthLocation(i[1] * u.deg, i[0] * u.deg, height=0 * u.m).get_itrs(t_start))
            for i in FibonacciPointGenerator(200)()
        ]

        # calculate line of sight access from points around the earth to the sun
        # when obstructed by the moon.
        solar_access_interval = None
        t0 = time.time()
        for earth_point in earth_positions:
            solar_eclipse_times = (
                Access(
                    LineOfSight(
                        earth_point,
                        sun_pos,
                        moon_pos,
                        sma=1079.6 * u.km,  # moon radius
                        smi=1079.6 * u.km,  # moon radius
                        when_obstructed=True,
                        use_frame="gcrs",
                    )
                )
                .use_precise_endpoints(True)
                .set_precision(1 * u.s)
                .calculate_at(times)
            )

            # union the intervals together to get cumulative time
            if solar_access_interval is None:
                solar_access_interval = solar_eclipse_times
            else:
                solar_access_interval = solar_access_interval | solar_eclipse_times
        t1 = time.time()

        print(f"Solar eclipse access calculation took: {t1-t0} seconds")
        print(f"There will be ~{solar_access_interval.total_duration} seconds of solar eclipse in 2024.")
        print(solar_access_interval)
        self.assertTrue(len(solar_access_interval) >= 2)

    def test_coverage(self):
        t_start = Time("2024-02-01T00:00:00")
        t_end = Time("2024-02-02T00:00:00")
        times = np.linspace(t_start, t_end, 8640)

        sat_position = get_ephem_as_skypath()

        worldwide_coverage_result = calculate_coverage([sat_position], times)
        print(f"Worldwide min revisit time: {worldwide_coverage_result.min_revisit_time}")
        print(f"Worldwide max revisit time: {worldwide_coverage_result.max_revisit_time}")

        at_equator_coverage = worldwide_coverage_result.filter(
            GeoFilter(min_latitude=-1 * u.deg, max_latitude=1 * u.deg)
        )
        print(f"Equatorial min revisit time: {at_equator_coverage.min_revisit_time}")
        print(f"Equatorial max revisit time: {at_equator_coverage.max_revisit_time}")
