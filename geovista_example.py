import time

import astropy.units as u
import geovista as gv
import numpy as np
import pyvista as pv
from astropy.time import Time

from skywatch.coordinates import SkyPath
from skywatch.tests.tests import get_ephem_as_skypath
from skywatch.utils import coverage

if __name__ == "__main__":
    t_start = Time("2024-02-01T00:00:00")
    t_end = Time("2024-02-02T00:00:00")
    times = np.linspace(t_start, t_end, 8640)
    sat_path = get_ephem_as_skypath()

    reversed_sat_path = SkyPath(
        sat_path.cartesian[::-1],
        representation_type="cartesian",
        frame="gcrs",
        obstime=sat_path.obstime,
    )

    t0 = time.time()
    coverage_result = coverage.calculate_coverage(
        [sat_path, reversed_sat_path] * 1,
        times,
        2500,
        use_precise_endpoints=False,
    )
    t1 = time.time()
    print(f"Coverage calculation took: {t1-t0} seconds")

    t0 = time.time()
    max_num_ball = coverage_result.get_max_simulatenous_visibility()
    print(f"Max Num Ball: {max_num_ball}")
    t1 = time.time()
    print(f"Max num ball took: {t1-t0} seconds")

    t0 = time.time()
    data, points = [], []
    for point in coverage_result:
        xyz = (
            SkyPath.from_geodetic(times[0], point.latitude, point.longitude, 0 * u.m)
            .cartesian.xyz.to(u.km)
            .value
        )
        xyz[2] += 0.01
        points.append(xyz / (np.linalg.norm(xyz)))
        data.append(point.cumulative_access.value)

    # Create a GeoPlotter object
    plotter = gv.GeoPlotter()
    m = pv.PolyData(points).delaunay_3d()
    m["data"] = data
    plotter.add_mesh(m, smooth_shading=True, scalars=data)

    plotter.add_coastlines(zlevel=1, color="black")
    # plotter.add_base_layer(texture=gv.natural_earth_hypsometric())

    t1 = time.time()
    print(f"Geovista plotting took: {t1-t0} seconds")

    # Display the plot
    plotter.show()
