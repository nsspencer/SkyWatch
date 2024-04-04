import math
from typing import List

import astropy.units as u
import numpy as np


def fibonacci_latitude_longitude(samples=1000) -> List[tuple]:
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        # Convert Cartesian coordinates to latitude and longitude
        lat = math.asin(y) * 180 / math.pi
        lon = math.atan2(z, x) * 180 / math.pi

        points.append((lat, lon))

    return points


def lat_lon_to_xyz(latitude, longitude, radius=6371):
    """
    Convert latitude and longitude to Cartesian coordinates.

    Parameters:
    latitude (float): Latitude in degrees.
    longitude (float): Longitude in degrees.
    radius (float): Radius of the sphere. Defaults to Earth's radius (6371 km).

    Returns:
    tuple: A tuple representing the Cartesian coordinates (X, Y, Z).
    """
    # Convert latitude and longitude from degrees to radians
    lat_rad = math.radians(latitude)
    lon_rad = math.radians(longitude)

    # Calculate Cartesian coordinates
    x = radius * math.cos(lat_rad) * math.cos(lon_rad)
    y = radius * math.cos(lat_rad) * math.sin(lon_rad)
    z = radius * math.sin(lat_rad)

    return x, y, z


@u.quantity_input(altitude=u.km, body_radius=u.km)
def max_elevation_angle(altitude: u.km, body_radius: u.km = 6371 * u.km):
    """Calculates the maximum elevation angle from NADIR
    a satellite can point to and still intersect the earth.
    """
    psi = np.arccos(body_radius / (body_radius + altitude))  # Horizon angle in radians
    e = (np.pi / 2) * u.rad - psi  # Maximum elevation angle in radians
    return np.degrees(e)  # Convert to degrees


# https://stephenhartzell.medium.com/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6
def los_to_earth(position, pointing):
    """Find the intersection of a pointing vector with the Earth
    Finds the intersection of a pointing vector u and starting point s with the WGS-84 geoid
    Args:
        position (np.array): length 3 array defining the starting point location(s) in meters
        pointing (np.array): length 3 array defining the pointing vector(s) (must be a unit vector)
    Returns:
        np.array: length 3 defining the point(s) of intersection with the surface of the Earth in meters

    Usage:
        To get the ground track of a satellite in cartesian coordinates, you would do the following calculation:
        ```
        los_to_earth(sat_pos_ecef.T.to(u.m).value, -sat_pos_ecef.T / np.linalg.norm(sat_pos_ecef.T))
        ```
        This gets the unit vector in the Z (nadir) direction for an ECEF position vector, and returns
        the position on the surface of the earth (in cartesian coords) that intersects the Z axis.
    """

    a = 6378137.0
    b = 6378137.0
    c = 6356752.314245
    x = position[0]
    y = position[1]
    z = position[2]
    u = pointing[0]
    v = pointing[1]
    w = pointing[2]

    value = -(a**2) * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x
    radical = (
        a**2 * b**2 * w**2
        + a**2 * c**2 * v**2
        - a**2 * v**2 * z**2
        + 2 * a**2 * v * w * y * z
        - a**2 * w**2 * y**2
        + b**2 * c**2 * u**2
        - b**2 * u**2 * z**2
        + 2 * b**2 * u * w * x * z
        - b**2 * w**2 * x**2
        - c**2 * u**2 * y**2
        + 2 * c**2 * u * v * x * y
        - c**2 * v**2 * x**2
    )
    magnitude = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2

    if radical < 0:
        raise ValueError("The Line-of-Sight vector does not point toward the Earth")
    d = (value - a * b * c * np.sqrt(radical)) / magnitude

    if d < 0:
        raise ValueError("The Line-of-Sight vector does not point toward the Earth")

    return np.array(
        [
            x + d * u,
            y + d * v,
            z + d * w,
        ]
    )
