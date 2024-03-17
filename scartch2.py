import time

import astropy.units as u
import numpy as np
import tqdm
from astropy.time import Time

from skypath.tests.tests import get_ephem_data
from skywatch.access import Access
from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.coordinates import SkyPath


def test1():
    t_start = Time("2024-02-01T00:00:00")
    t_end = Time("2024-02-02T00:00:00")
    low_fidelity_times = np.linspace(t_start, t_end, 1440)
    high_fidelity_times = np.linspace(t_start, t_end, 8640)

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

    gs_to_sat_access_algorithm = (
        Access(
            LineOfSight(ground_station_pos, sat_position, earth_pos),
            # AzElRange(ground_station_pos, sat_position, min_el=25 * u.deg),
            # Temporal(low_fidelity_times[10], low_fidelity_times[700], True),
            LineOfSight(sat_position, sun_pos, earth_pos),
            # AzElRange(ground_station_pos, sun_pos, max_el=10 * u.deg),
            LineOfSight(sat_position, moon_pos, earth_pos),
            # AzElRange(ground_station_pos, moon_pos, min_el=0 * u.deg),
        )
        .use_precise_endpoints(True)
        .set_precision(0.1 * u.s)
    )

    for _ in tqdm.tqdm(range(100)):
        access_times = gs_to_sat_access_algorithm(high_fidelity_times)
    print(access_times.total_duration)
    print(access_times)
    pass


if __name__ == "__main__":
    test1()
