import datetime
import time
import unittest
from typing import Union

import astropy.units as u
import numpy as np
from astropy.coordinates import BaseCoordinateFrame, SkyCoord
from astropy.time import Time, TimeDelta
from scipy.spatial.transform import Rotation

from skywatch.access import Access, TimeInterval
from skywatch.access.constraints import AzElRange, LineOfSight, Temporal
from skywatch.attitude import LVLH
from skywatch.coordinates import SkyPath
from skywatch.look_angles import LookAngles
from skywatch.tests.tests import get_ephem_as_skypath
from skywatch.utils.coverage import GeoFilter, calculate_coverage


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

    def test_look_angles(self):
        t_start = Time("2024-02-01T00:00:00")
        t_end = Time("2024-02-02T00:00:00")
        times = np.linspace(t_start, t_end, 8640)

        sat_position = get_ephem_as_skypath()
        ground_station_pos = SkyPath.from_geodetic(
            times, 34.5 * u.deg, -77.0 * u.deg, 0 * u.m
        )

        gs_look_angles_fn = LookAngles.from_local_enu(ground_station_pos, False)
        gs_look_angles_fn_astropy = LookAngles.from_local_enu(ground_station_pos, True)
        gs_to_sat_look_angles = gs_look_angles_fn(times, sat_position)
        gs_to_sat_look_angles_astropy = gs_look_angles_fn_astropy(times, sat_position)

        sat_look_angles_fn = LookAngles.from_body_frame(
            sat_position, LVLH(sat_position)
        )
        sat_look_angles_fn_offset = LookAngles.from_body_frame(
            sat_position,
            LVLH(sat_position),
            Rotation.from_euler("z", -45, degrees=True),
        )
        sat_to_ground_station_look_angles = sat_look_angles_fn(
            times, ground_station_pos
        )
        sat_to_ground_station_look_angles_offset = sat_look_angles_fn_offset(
            times, ground_station_pos
        )

        gs_body_frame = LookAngles.from_body_frame(
            ground_station_pos, LVLH(ground_station_pos, "itrs")
        )
        plot_ground_track_and_points(
            sat_position,
            [
                ground_station_pos.itrs.earth_location.lat.value[0],
                ground_station_pos.itrs.earth_location.lon.value[0],
            ],
        )

        plot_look_angles(
            sat_to_ground_station_look_angles_offset.azimuth,
            sat_to_ground_station_look_angles_offset.elevation,
            sat_to_ground_station_look_angles_offset.range,
            sat_to_ground_station_look_angles_offset.time,
            "Offset plot",
        )

        plot_look_angles(
            sat_to_ground_station_look_angles.azimuth,
            sat_to_ground_station_look_angles.elevation,
            sat_to_ground_station_look_angles.range,
            sat_to_ground_station_look_angles.time,
            "LVLH plot",
        )

        try:
            gs_body_to_sat = gs_body_frame(times, sat_position)
            self.assertTrue(False)  # fail
        except AttributeError:
            self.assertTrue(True)  # pass because static position has no velocity


# class BodyFrameBase:
#     def __init__(self, reference_frame: str) -> None:
#         self.reference_frame = reference_frame

#     def calculate(self, coordinate: SkyPath) -> Rotation:
#         pass


# class LocalTangent(BodyFrameBase):
#     pass


# class LVLH(BodyFrameBase):
#     pass


# class Platform:
#     def __init__(self, coordinates: SkyPath, body_frame: BodyFrameBase = None) -> None:
#         self.__body_frame = body_frame
#         self.__coordinates = coordinates

#     def set_body_frame(self, body_frame: BodyFrameBase) -> "Platform":
#         if not isinstance(body_frame, BodyFrameBase):
#             raise TypeError("body_frame must subclass BodyFrameBase.")
#         self.__body_frame = body_frame
#         return self

#     def set_coordinates(
#         self,
#         coordinates: Union[SkyPath, SkyCoord, BaseCoordinateFrame],
#     ) -> "Platform":
#         if isinstance(coordinates, SkyPath):
#             self.__coordinates = coordinates
#         else:
#             self.__coordinates = SkyPath(coordinates)
#         return self

#     def look_angles(
#         self,
#         time: Time,
#         target: Union["Platform", SkyPath],
#         attitude_offset: Rotation = None,
#     ) -> AzElRange:
#         if self.__coordinates is None:
#             raise AttributeError(
#                 "You must define a body frame for this platform before you can use the look_angles function."
#             )


def plot_look_angles(
    azimuth: np.ndarray,
    elevation: np.ndarray,
    range: np.ndarray,
    times: np.ndarray,
    title: str = "Look Angles Plot",
) -> None:
    """Simple plotter to show the azimuth, elevation, and range look angles over time.

    Args:
        azimuth (np.ndarray): azimuth values
        elevation (np.ndarray): elevation values
        range (np.ndarray): range values
        times (np.ndarray): associated time values
    """
    import plotly.graph_objects as go

    times = [t.to_datetime() for t in times]
    # Create a trace for each line
    trace1 = go.Scatter(
        x=times, y=azimuth, mode="lines", name="Azimuth", line=dict(color="green")
    )
    trace2 = go.Scatter(
        x=times, y=elevation, mode="lines", name="Elevation", line=dict(color="blue")
    )
    trace3 = go.Scatter(
        x=times, y=range, mode="lines", name="Range", yaxis="y2", line=dict(color="red")
    )

    data = [trace1, trace2, trace3]

    layout = go.Layout(
        title=title,
        yaxis=dict(title="Angle", side="left"),
        yaxis2=dict(title="Distance", overlaying="y", side="right"),
        xaxis=dict(
            title="Date",
            showgrid=True,
            zeroline=True,
            showline=True,
            linecolor="rgb(204, 204, 204)",
            linewidth=2,
            showticklabels=True,
            tickformat="%H:%M",
        ),
    )

    fig = go.Figure(data=data, layout=layout)
    fig.show()


def plot_ground_track_and_points(satellite: SkyPath, observation_llas: list = None):
    import plotly.graph_objects as go

    # Create the line trace
    sat_earth_pos = satellite.transform_to("itrs").earth_location
    lats = sat_earth_pos.lat.value
    lons = sat_earth_pos.lon.value
    ground_track = go.Scattergeo(
        lat=lats,
        lon=lons,
        mode="lines",
        line=dict(width=2, color="blue"),
        name="Line",
        text=[str(t) for t in satellite.obstime],
    )

    data = [ground_track]

    if observation_llas != None:
        if isinstance(observation_llas[0], list):
            observer_lats = [i[0] for i in observation_llas]
            observer_lons = [i[1] for i in observation_llas]
        else:
            observer_lats = [observation_llas[0]]
            observer_lons = [observation_llas[1]]
            print(
                f"plotting single point: [{observation_llas[0]}, {observation_llas[1]}]"
            )

        # Create the points trace
        observers = go.Scattergeo(
            lat=observer_lats,
            lon=observer_lons,
            mode="markers",
            marker=dict(size=20, color="red"),
            name="Points",
        )

        data.append(observers)

    # Create the plot
    fig = go.Figure(data=data)

    # Update the layout
    fig.update_layout(
        title_text="Line and Points on the Earth",
        showlegend=True,
        geo=dict(
            resolution=50,
            showland=True,
            showlakes=True,
            landcolor="rgb(204, 204, 204)",
            countrycolor="rgb(204, 204, 204)",
            lakecolor="rgb(255, 255, 255)",
            projection_type="equirectangular",
            coastlinewidth=2,
            lataxis=dict(range=[-90, 90], showgrid=True, dtick=10),
            lonaxis=dict(range=[-180, 180], showgrid=True, dtick=20),
        ),
    )

    fig.show()


if __name__ == "__main__":
    unittest.main()
