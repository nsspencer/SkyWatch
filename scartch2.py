import time

import astropy.units as u
import numpy as np
from astropy.time import Time

from skypath.tests.tests import get_ephem_data
from skywatch.access import Access
from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.coordinates import SkyPath


def test1():
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

    leo_csv = get_ephem_data()
    leo_csv_times = Time(leo_csv["utc"])
    leo_csv_position = (
        np.array(
            [leo_csv["x_eci_km"], leo_csv["y_eci_km"], leo_csv["z_eci_km"]]
        ).astype(float)
        * u.km
    )
    leo_csv_velocities = np.array(
        [leo_csv["vx_eci_km_s"], leo_csv["vy_eci_km_s"], leo_csv["vz_eci_km_s"]]
    ).astype(float) * (u.km / u.s)

    sat_position = SkyPath.from_ECI(
        leo_csv_times, *leo_csv_position, leo_csv_velocities
    )

    t0 = time.time()
    access_times = (
        Access(
            # when my ground station can see my satellite without being obstructed by the earth
            LineOfSight(ground_station_pos, sat_position, earth_pos),
            # when the sun angle is low to the horizon
            AzElRange(ground_station_pos, sun_pos, min_el=0 * u.deg, max_el=5 * u.deg),
            # when my satellite can see the sun
            LineOfSight(sat_position, sun_pos, earth_pos),
            # when its after 10:15pm but before 10:17pm
            Temporal(
                min_time=Time("2024-02-01T22:15:00"),
                max_time=Time("2024-02-01T22:17:00"),
            ),
            # when the moon can see the sun without being blocked by the earth
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


def test2():
    """
    calculate approximate lunar and eclipse times
    """

    t_start = Time("2024-01-01T00:00:00")
    t_end = Time("2025-01-01T00:00:00")
    times = np.linspace(t_start, t_end, 8640 * 2)

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
                sma=1079.6 * 6 * u.km,  # use a body size ~ the size of the earth
                smi=1079.6 * 6 * u.km,  # use a body size ~ the size of the earth
                when_obstructed=True,
                use_frame="icrs",
            )
        )
        .use_precise_endpoints(True)
        .calculate_at(times)
    )
    t1 = time.time()

    print(f"Solar eclipse access calculation took: {t1-t0} seconds")
    print(solar_eclipse_times.total_duration)
    print(solar_eclipse_times)


if __name__ == "__main__":
    test1()
    test2()
