import math
import time

import astropy.units as u
import numpy as np
import pandas as pd
import pymap3d
import tqdm
from astropy.time import Time

from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.look_angles import ENUAER, NadirWithVelocityConstraint
from skywatch.skypath import SkyPath
from skywatch.tests.tests import get_ephem_data
from skywatch.utils.funcs import fibonacci_latitude_longitude

if __name__ == "__main__":
    scenario_start = Time("2024-02-01T00:00:00.000")
    scenario_end = Time("2024-02-02T00:00:00.000")
    num_elements = int(8640)
    times = np.linspace(scenario_start, scenario_end, num_elements)

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

    earth_positions = []
    num_fib_points = 1000
    point_times = np.linspace(times[0], times[-1], 250)
    for lat, lon in tqdm.tqdm(
        fibonacci_latitude_longitude(num_fib_points), desc="Making Facilities"
    ):
        earth_positions.append(
            SkyPath.from_geodetic(
                time=point_times,
                latitude=lat * u.deg,
                longitude=lon * u.deg,
                altitude=0 * u.m,
            )
        )

    sat_position = SkyPath.from_ECI(
        leo_csv_times, *leo_csv_position, leo_csv_velocities
    )
    sat_position: SkyPath
    interp_sat_position = sat_position.state_at(times[0], "itrs")
    interp_sat_position2 = interp_sat_position.state_at(times[0], "gcrs")

    constraints = []
    # constraints.append(Temporal(times[3000], times[4000], inner=True))
    constraints.append(AzElRange(min_el=40 * u.deg))
    constraints.append(LineOfSight(SkyPath.from_body(times, "earth"), use_frame="itrs"))
    constraints.append(
        LineOfSight(
            SkyPath.from_body(times, "moon"), sma=1737.1 * u.km, smi=1737.1 * u.km
        )
    )

    # calculate current sun az el from DC
    current_time = np.linspace(
        Time("2024-03-14 21:23:45.430137"), Time("2024-03-15 02:23:45.430137"), 100
    )
    sun_pos = SkyPath.from_body(current_time, "sun")
    earth_pos = SkyPath.from_geodetic(
        time=current_time,
        latitude=38.907192 * u.deg,
        longitude=-77.036873 * u.deg,
        altitude=0 * u.m,
    )
    earth_to_sun = earth_pos.look_angles_to(sun_pos, current_time, ENUAER(True))
    eart_sun_access = earth_pos.access_to(
        sun_pos,
        current_time,
        [LineOfSight(SkyPath.from_body(current_time, "earth"), use_frame="itrs")],
    )

    all_access = []
    for point in tqdm.tqdm(earth_positions, desc="Calculating access"):
        point: SkyPath
        # sat_access = sat_position.access_to(point, times, constraints, True, 0.1 * u.s)
        point_access = point.access_to(
            sat_position, times, constraints, True, 0.1 * u.s
        )

        look_angles = []
        for interval in point_access:
            l_times = np.linspace(interval.lower, interval.upper, 400)

            az, el, range = sat_position.look_angles_to(
                point, l_times, NadirWithVelocityConstraint()
            )
            # az1, el1, range1 = sat_position.look_angles_to(point, l_times, NadirWithVelocityConstraint(frame='itrs'))
            # point_az, point_el, point_range = point.look_angles_to(sat_position, l_times, ENUAER(high_fidelity=True))
            look_angles.append((l_times, az, el, range))

        # if len(sat_access) != len(point_access):
        #     raise ValueError("Not commutative")
        # for index in range(len(sat_access)):
        #     sat_interval = (sat_access._intervals[index].upper - sat_access._intervals[index].lower).datetime.total_seconds()
        #     point_interval = (point_access._intervals[index].upper - point_access._intervals[index].lower).datetime.total_seconds()
        #     if sat_interval != point_interval:
        #         raise ValueError("Not commutative")

        all_access.append((point, point_access, look_angles))

    access_seconds = []
    for point, access, look_angles in all_access:
        access_seconds.append(access.total_duration.value)

    print(pd.DataFrame(access_seconds).describe())

    points = np.linspace(0, 6437100, num_elements) * u.m
    velocities = np.zeros_like(points).value * (u.m / u.s)
    lat_points = np.linspace(-90, 90, num_elements)
    lon_points = np.linspace(-180, 180, num_elements)

    # position only

    t0 = time.time()
    earth_point = SkyPath.from_geodetic(scenario_start, 0 * u.deg, 0 * u.deg, 0 * u.m)
    print(f"single geodetic: {time.time() - t0}")

    t0 = time.time()
    multi_earth_points = SkyPath.from_geodetic(
        times, lat_points * u.deg, lon_points * u.deg, points
    )
    print(f"multiple geodetic: {time.time() - t0}")

    t0 = time.time()
    ecef_points = SkyPath.from_ECEF(times, points, points, points)
    print(f"ECEF: {time.time() - t0}")

    t0 = time.time()
    eci_points = SkyPath.from_ECI(times, points, points, points)
    print(f"ECI -> ECEF: {time.time() - t0}")

    # now with velocities

    t0 = time.time()
    ecef_points_vels = SkyPath.from_ECEF(
        times, points, points, points, velocities, velocities, velocities
    )
    print(f"ECEF w/ vel: {time.time() - t0}")

    t0 = time.time()
    eci_points_vels = SkyPath.from_ECI(
        times, points, points, points, velocities, velocities, velocities
    )
    print(f"ECI -> ECEF w/ vel: {time.time() - t0}")

    # now cross check computation time with pymap3d (without velocities)

    t0 = time.time()
    ecef_vals = pymap3d.eci2ecef(points.value, points.value, points.value, times)
    print(f"pymap3d ECI -> ECEF: {time.time() - t0}")

    t0 = time.time()
    eci_vals = pymap3d.ecef2eci(points.value, points.value, points.value, times)
    print(f"pymap3d ECEF -> ECI: {time.time() - t0}")

    # test with the single earth point
    cs0 = SkyPath(earth_point)
    cs0_itrs = cs0.state_at(times, "itrs")
    cs0_gcrs = cs0.state_at(times, "gcrs")
    cs0_icrs = cs0.state_at(times, "icrs")

    # now create a CoordinateInterpolator

    cs1 = SkyPath(ecef_points_vels)
    new_times = np.linspace(scenario_start, scenario_end, num_elements * 1)

    print(f"\n\nNum interp times: {len(new_times):,}")

    t0 = time.time()
    eci_pos = cs1.state_at(new_times, "gcrs")
    print(f"Initial to GCRS: {time.time() - t0}")

    t0 = time.time()
    eci_pos2 = cs1.state_at(new_times, "gcrs")
    print(f"Interp to GCRS: {time.time() - t0}")

    t0 = time.time()
    ecef_pos = cs1.state_at(new_times, "itrs")
    print(f"Initial to ITRS: {time.time() - t0}")

    t0 = time.time()
    ecef_pos2 = cs1.state_at(new_times, "itrs")
    print(f"Interp to ITRS: {time.time() - t0}")

    t0 = time.time()
    icrs_pos = cs1.state_at(new_times, "icrs")
    print(f"initial to ICRS: {time.time() - t0}")

    t0 = time.time()
    icrs_pos2 = cs1.state_at(new_times, "icrs")
    print(f"Interp to ICRS: {time.time() - t0}")

    # now test the geodedic point with no velocity

    print("\n\nGeodetic test no velocity: ")

    cs2 = SkyPath(multi_earth_points)
    print(f"\n\nNum interp times: {len(new_times):,}")

    t0 = time.time()
    eci_pos = cs2.state_at(new_times, "gcrs")
    print(f"Initial to GCRS: {time.time() - t0}")

    t0 = time.time()
    eci_pos2 = cs2.state_at(new_times, "gcrs")
    print(f"Interp to GCRS: {time.time() - t0}")

    t0 = time.time()
    ecef_pos = cs2.state_at(new_times, "itrs")
    print(f"Initial to ITRS: {time.time() - t0}")

    t0 = time.time()
    ecef_pos2 = cs2.state_at(new_times, "itrs")
    print(f"Interp to ITRS: {time.time() - t0}")

    t0 = time.time()
    icrs_pos = cs2.state_at(new_times, "icrs")
    print(f"initial to ICRS: {time.time() - t0}")

    t0 = time.time()
    icrs_pos2 = cs2.state_at(new_times, "icrs")
    print(f"Interp to ICRS: {time.time() - t0}")
    pass
