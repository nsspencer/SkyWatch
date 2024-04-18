from typing import Tuple

import astropy.units as u
from astropy.time import Time

from skywatch.coverage.funcs import fibonacci_latitude_longitude_generator
from skywatch.skypath import SkyPath


class CoverageDefinitionMixin:

    @classmethod
    def from_skypath(cls, skypath: SkyPath):
        return cls(skypath)

    @classmethod
    def from_geodetic(
        cls,
        time: Time,
        latitude: u.Quantity,
        longitude: u.Quantity,
        altitude: u.Quantity,
    ):
        return cls(SkyPath.from_geodetic(time, latitude, longitude, altitude))

    @classmethod
    def from_fibonacci_earth(
        cls,
        time: Time,
        num_samples: int = 1000,
        altitude: u.Quantity = 0 * u.m,
    ):
        points = (
            SkyPath.from_geodetic(time, lat, lon, altitude)
            for lat, lon in fibonacci_latitude_longitude_generator(num_samples)
        )
        return cls(*points)

    #######
    # TODO: Implement the following:
    #
    # @classmethod
    # def from_country(cls, time: Time, name: str, num_samples: int):
    #     raise NotImplementedError("TODO")

    # @classmethod
    # def from_us_state(cls, time: Time, name: str):
    #     raise NotImplementedError("TODO")

    # @classmethod
    # def from_geojson(cls, time: Time, geojson: Any):
    #     raise NotImplementedError("TODO and use for states, countries, etc.")

    # @classmethod
    # def from_body(
    #     cls,
    #     time: Time,
    #     name: str = "earth",
    #     num_samples: int = 1000,
    #     radius: u.Quantity = 6371 * u.km,
    # ):
    #     raise NotImplementedError(
    #         "This requires more thought, like spin rate, axial tilt (?), clockwise/counter-clockwise rotation."
    #     )

    #     # define a sample of position around the surface of the moon
    #     body_position_offsets = [
    #         lat_lon_to_xyz(lat, lon, radius=radius.value)
    #         for lat, lon in fibonacci_latitude_longitude_generator(num_samples)
    #     ]

    #     # get the body position and add the offsets to the position to define
    #     # coordinates on the body's surface
    #     body_pos = SkyPath.from_body(time, name)
    #     body_surface_points = (
    #         SkyPath(
    #             CartesianRepresentation(
    #                 *(body_pos.cartesian.xyz.T + (body_pos_offset * u.m)).T
    #             ),
    #             frame="itrs",
    #             obstime=time,
    #         )
    #         for body_pos_offset in body_position_offsets
    #     )

    #     return cls(*body_surface_points)


class CoverageDefinition(CoverageDefinitionMixin):
    def __init__(self, *points: SkyPath) -> None:
        self.points = list()
        for point in points:
            if not isinstance(point, SkyPath):
                raise TypeError(
                    "All points in coverage definition must be a SkyPath object."
                )
            self.points.append(point)

    def calculate(
        self,
        time: Time,
        observers: Tuple[SkyPath],
        min_elevation: u.Quantity = 0 * u.deg,
    ) -> None:
        pass
