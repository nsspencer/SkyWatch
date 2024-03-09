from astro_access import FrameInterpolator
from astro_access import AccessAlgorithm
from astro_access.constraints import LineOfSight


import astropy.units as u
from astropy.time import Time
import numpy as np
import time

from astropy.coordinates import SkyCoord, EarthLocation
import pymap3d


@u.quantity_input(latitude=u.deg, longitude=u.deg, altitude=u.m)
def from_geodetic(time: Time,
                  latitude: u.deg, longitude: u.deg, altitude: u.m) -> SkyCoord:
    return SkyCoord(EarthLocation(lat=latitude, lon=longitude, height=altitude).get_itrs(obstime=time))


@u.quantity_input(x=u.m, y=u.m, z=u.m, v_x=u.m/u.s, v_y=u.m/u.s, v_z=u.m/u.s)
def from_ecef(time: Time,
              x: u.m, y: u.m, z: u.m,
              v_x: u.m/u.s = None, v_y: u.m/u.s = None, v_z: u.m/u.s = None) -> SkyCoord:
    return SkyCoord(x=x, y=y, z=z,
                    v_x=v_x, v_y=v_y, v_z=v_z,
                    frame='itrs', representation_type='cartesian',
                    obstime=time)


@u.quantity_input(x=u.m, y=u.m, z=u.m, v_x=u.m/u.s, v_y=u.m/u.s, v_z=u.m/u.s)
def from_eci(time: Time,
             x: u.m, y: u.m, z: u.m,
             v_x: u.m/u.s = None, v_y: u.m/u.s = None, v_z: u.m/u.s = None) -> SkyCoord:
    return SkyCoord(x=x, y=y, z=z,
                    v_x=v_x, v_y=v_y, v_z=v_z,
                    frame='gcrs', representation_type='cartesian',
                    obstime=time)




if __name__ == "__main__":
    scenario_start = Time("2024-03-06T00:00:00")
    scenario_end = Time("2024-03-07T00:00:00")

    num_elements = int(86400/10)
    points = np.linspace(0, 6437100, num_elements) * u.m
    velocities = np.zeros_like(points).value * (u.m/u.s)
    lat_points = np.linspace(-90, 90, num_elements)
    lon_points = np.linspace(-180, 180, num_elements)
    times = np.linspace(scenario_start, scenario_end, num_elements)
    
    # position only
    
    t0 = time.time()
    earth_point = from_geodetic(scenario_start, 0 * u.deg, 0 * u.deg, 0 * u.m)
    print(f"single geodetic: {time.time() - t0}")
    
    t0 = time.time()
    multi_earth_points = from_geodetic(times, lat_points * u.deg, lon_points * u.deg, points)
    print(f"multiple geodetic: {time.time() - t0}")

    t0 = time.time()
    ecef_points = from_ecef(times, points, points, points)
    print(f"ECEF: {time.time() - t0}")
    
    t0 = time.time()
    eci_points = from_eci(times, points, points, points)
    print(f"ECI -> ECEF: {time.time() - t0}")

    # now with velocities
    
    t0 = time.time()
    ecef_points_vels = from_ecef(times, points, points, points, velocities, velocities, velocities)
    print(f"ECEF w/ vel: {time.time() - t0}")
    
    t0 = time.time()
    eci_points_vels = from_eci(times, points, points, points, velocities, velocities, velocities)
    print(f"ECI -> ECEF w/ vel: {time.time() - t0}")


    # now cross check computation time with pymap3d (without velocities)

    t0 = time.time()
    ecef_vals = pymap3d.eci2ecef(points.value, points.value, points.value, times)
    print(f"pymap3d ECI -> ECEF: {time.time() - t0}")
    
    t0 = time.time()
    eci_vals = pymap3d.ecef2eci(points.value, points.value, points.value, times)
    print(f"pymap3d ECEF -> ECI: {time.time() - t0}")
    
    
    # test with the single earth point
    cs0 = FrameInterpolator(earth_point)
    cs0_itrs = cs0.state_at(times, 'itrs')
    cs0_gcrs = cs0.state_at(times, 'gcrs')
    cs0_icrs = cs0.state_at(times, 'icrs')
    
    # now create a FrameInterpolator
    
    cs1 = FrameInterpolator(ecef_points_vels)
    new_times = np.linspace(scenario_start, scenario_end, num_elements*1)
    
    print(f"\n\nNum interp times: {len(new_times):,}")
    
    t0 = time.time()
    eci_pos = cs1.state_at(new_times, 'gcrs')
    print(f"Initial to GCRS: {time.time() - t0}")
    
    t0 = time.time()
    eci_pos2 = cs1.state_at(new_times, 'gcrs')
    print(f"Interp to GCRS: {time.time() - t0}")
    
    t0 = time.time()
    ecef_pos = cs1.state_at(new_times, 'itrs')
    print(f"Initial to ITRS: {time.time() - t0}")
    
    t0 = time.time()
    ecef_pos2 = cs1.state_at(new_times, 'itrs')
    print(f"Interp to ITRS: {time.time() - t0}")
    
    t0 = time.time()
    icrs_pos = cs1.state_at(new_times, 'icrs')
    print(f"initial to ICRS: {time.time() - t0}")
    
    t0 = time.time()
    icrs_pos2 = cs1.state_at(new_times, 'icrs')
    print(f"Interp to ICRS: {time.time() - t0}")
    
    # now test the geodedic point with no velocity
    
    print("\n\nGeodetic test no velocity: ")
    
    cs2 = FrameInterpolator(multi_earth_points)
    print(f"\n\nNum interp times: {len(new_times):,}")
    
    t0 = time.time()
    eci_pos = cs2.state_at(new_times, 'gcrs')
    print(f"Initial to GCRS: {time.time() - t0}")
    
    t0 = time.time()
    eci_pos2 = cs2.state_at(new_times, 'gcrs')
    print(f"Interp to GCRS: {time.time() - t0}")
    
    t0 = time.time()
    ecef_pos = cs2.state_at(new_times, 'itrs')
    print(f"Initial to ITRS: {time.time() - t0}")
    
    t0 = time.time()
    ecef_pos2 = cs2.state_at(new_times, 'itrs')
    print(f"Interp to ITRS: {time.time() - t0}")
    
    t0 = time.time()
    icrs_pos = cs2.state_at(new_times, 'icrs')
    print(f"initial to ICRS: {time.time() - t0}")
    
    t0 = time.time()
    icrs_pos2 = cs2.state_at(new_times, 'icrs')
    print(f"Interp to ICRS: {time.time() - t0}")
    pass