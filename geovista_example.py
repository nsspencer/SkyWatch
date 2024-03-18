import astropy.units as u
import geovista as gv
import numpy as np
from astropy.time import Time

from skywatch.coordinates import SkyPath
from skywatch.tests.tests import get_ephem_as_skypath
from skywatch.utils import coverage

if __name__ == "__main__":
    t_start = Time("2024-02-01T00:00:00")
    t_end = Time("2024-02-02T00:00:00")
    times = np.linspace(t_start, t_end, 8640)
    sat_path = get_ephem_as_skypath()
    coverage_result = coverage.calculate_coverage(
        [sat_path], times, 2500, use_precise_endpoints=False
    )

    # Generate Fibonacci latitude and longitude points
    data, points = [], []
    for point in coverage_result:
        xyz = (
            SkyPath.from_geodetic(times[0], point.latitude, point.longitude, 0 * u.m)
            .cartesian.xyz.to(u.km)
            .value
        )
        points.append(xyz / (np.linalg.norm(xyz)))
        data.append(point.interval.total_duration.value)

    # Create a GeoPlotter object
    plotter = gv.GeoPlotter()

    # Add the points to the plotter
    plotter.add_points(
        points=np.array(points), scalars=data, point_size=25, style="points"
    )
    plotter.add_coastlines(zlevel=1, color="black")
    plotter.add_base_layer(texture=gv.natural_earth_hypsometric())

    # Display the plot
    plotter.show()
