import unittest

import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude import LVLH
from skywatch.look_angles import LocalTangentENU, XForwardZNadir
from skywatch.skypath import SkyPath

from .tests import get_ephem_as_skypath


def angle_differences(angles1, angles2):
    """
    Calculate the differences between corresponding angles in two arrays.
    """
    # Normalize angles to be between 0 and 360 degrees
    angles1 = angles1 % 360
    angles2 = angles2 % 360

    # Calculate the absolute difference between the angles
    diffs = np.abs(angles1 - angles2)

    return diffs


class Test1(unittest.TestCase):
    def test_ground_station_to_ground_station(self):
        gs1 = SkyPath.from_geodetic(Time("J2000"), 10 * u.deg, 10 * u.deg, 0 * u.m)
        gs2 = SkyPath.from_geodetic(Time("J2000"), 19 * u.deg, 1 * u.deg, 155000 * u.m)

        gs1_look_angles = LocalTangentENU(gs1)
        gs2_look_angles = LocalTangentENU(gs2)
        gs1_to_gs2 = gs1_look_angles.get_look_angles(gs2, Time("J2000"))
        gs2_to_gs1 = gs2_look_angles.get_look_angles(gs1, Time("J2000"))
        print(gs1_to_gs2)
        print(gs2_to_gs1)

    def test_ground_station_to_satellite(self):
        sat = get_ephem_as_skypath()
        gs1 = SkyPath.from_geodetic(Time("J2000"), 10 * u.deg, 10 * u.deg, 0 * u.m)
        gs2 = SkyPath.from_geodetic(sat.obstime[0], 10 * u.deg, 10 * u.deg, 0 * u.m)

        gs_look_angles = LocalTangentENU(gs1)
        gs2sat = gs_look_angles.get_look_angles(sat, sat.obstime)

        gs2_look_angles = LocalTangentENU(gs2)
        gs2sat2 = gs2_look_angles.get_look_angles(sat, sat.obstime)

        self.assertTrue(np.allclose(gs2sat.azimuth, gs2sat2.azimuth))
        self.assertTrue(np.allclose(gs2sat.elevation, gs2sat2.elevation))
        self.assertTrue(np.allclose(gs2sat.range, gs2sat2.range))

    def test_satellite_to_ground_station(self):
        sat = get_ephem_as_skypath()
        ground_station = SkyPath.from_geodetic(
            sat.obstime, 34.5 * u.deg, -77.0 * u.deg, 0 * u.m
        )

        sat_look_angles_fn = XForwardZNadir(sat, LVLH(sat, "itrs"), "itrs")
        sat_look_angles = sat_look_angles_fn.get_look_angles(
            ground_station, sat.obstime
        )

        gs_look_angles_fn_gcrs = XForwardZNadir(sat, LVLH(sat, "gcrs"), "gcrs")
        sat_look_angles_gcrs = gs_look_angles_fn_gcrs.get_look_angles(
            ground_station, sat.obstime
        )

        self.assertTrue(
            np.min(
                angle_differences(
                    sat_look_angles.azimuth.value, sat_look_angles_gcrs.azimuth.value
                )
            )
            < 2.0
        )

        self.assertTrue(
            np.allclose(sat_look_angles.elevation, sat_look_angles_gcrs.elevation)
        )
        self.assertTrue(np.allclose(sat_look_angles.range, sat_look_angles_gcrs.range))

        sat_look_angles_offset_fn = XForwardZNadir(
            sat, LVLH(sat, "itrs", Rotation.from_euler("z", -45, degrees=True))
        )
        sat_2_gs_offset = sat_look_angles_offset_fn.get_look_angles(
            ground_station, sat.obstime
        )
        pass

    def test_look_angles(self):
        t_start = Time("2024-02-01T22:10:00")
        t_end = Time("2024-02-01T22:17:00")
        times = np.linspace(t_start, t_end, 8640)

        sat_position = get_ephem_as_skypath()
        ground_station_pos = SkyPath.from_geodetic(
            times, 34.5 * u.deg, -77.0 * u.deg, 0 * u.m
        )

        # gs_look_angles_fn = LocalTangentENU(ground_station_pos, False)
        # gs_look_angles_fn_astropy = LocalTangentENU(ground_station_pos, True)
        # gs_to_sat_look_angles = gs_look_angles_fn.get_look_angles(sat_position, times)
        # gs_to_sat_look_angles_astropy = gs_look_angles_fn_astropy.get_look_angles(
        #     sat_position, times
        # )

        sat_look_angles_fn = XForwardZNadir(sat_position, LVLH(sat_position))
        sat_look_angles_fn_offset = XForwardZNadir(
            sat_position,
            LVLH(
                sat_position,
                offset=Rotation.from_euler("z", -90, degrees=True),
            ),
        )
        sat_to_ground_station_look_angles = sat_look_angles_fn.get_look_angles(
            ground_station_pos, times
        )
        sat_to_ground_station_look_angles_offset = (
            sat_look_angles_fn_offset.get_look_angles(ground_station_pos, times)
        )

        extra_text = [
            f"norm: {az1}\noffset: {az2}"
            for az1, az2 in zip(
                sat_to_ground_station_look_angles.azimuth,
                sat_to_ground_station_look_angles_offset.azimuth,
            )
        ]

        plot_ground_track_and_points(
            times,
            sat_position,
            [
                ground_station_pos.itrs.earth_location.lat.value[0],
                ground_station_pos.itrs.earth_location.lon.value[0],
            ],
            extra_text=extra_text,
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


def plot_ground_track_and_points(
    time: Time,
    satellite: SkyPath,
    observation_llas: list = None,
    extra_text: list = None,
):
    import plotly.graph_objects as go

    # Create the line trace
    sat_earth_pos = satellite.state_at(time, "itrs").earth_location
    lats = sat_earth_pos.lat.value
    lons = sat_earth_pos.lon.value
    if extra_text is not None:
        text_values = [f"{str(t)}\n{str(extra)}" for t, extra in zip(time, extra_text)]
    else:
        text_values = [str(t) for t in time]

    ground_track = go.Scattergeo(
        lat=lats,
        lon=lons,
        mode="lines",
        line=dict(width=2, color="blue"),
        name="Line",
        text=text_values,
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
