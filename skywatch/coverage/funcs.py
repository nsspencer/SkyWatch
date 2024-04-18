from math import asin, atan2, cos, pi, radians, sin, sqrt
from typing import Generator

import numpy as np
from shapely import Point, Polygon


def fibonacci_latitude_longitude_generator(
    samples=1000,
) -> Generator[float, float, float]:
    phi = pi * (3.0 - sqrt(5.0))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = cos(theta) * radius
        z = sin(theta) * radius

        # Convert Cartesian coordinates to latitude and longitude
        lat = asin(y) * 180 / pi
        lon = atan2(z, x) * 180 / pi

        yield lat, lon


def lat_lon_to_xyz(latitude: float, longitude: float, radius: float = 6371):
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
    lat_rad = radians(latitude)
    lon_rad = radians(longitude)

    # Calculate Cartesian coordinates
    x = radius * cos(lat_rad) * cos(lon_rad)
    y = radius * cos(lat_rad) * sin(lon_rad)
    z = radius * sin(lat_rad)

    return x, y, z


def points_in_polygon(polygon: Polygon, N: int):
    raise NotImplementedError("Needs testing")

    # Define the golden ratio
    golden_ratio = (np.sqrt(5) - 1) / 2

    # Get the area of the polygon
    area = polygon.area

    # Calculate the average area per point
    area_per_point = area / N

    # Calculate the average distance between points
    average_distance = np.sqrt(area_per_point)

    # Get the bounds of the polygon
    minx, miny, maxx, maxy = polygon.bounds

    # Calculate the number of points in the x and y directions
    nx = int((maxx - minx) / average_distance)
    ny = int((maxy - miny) / average_distance)

    # Create a list to store the points
    points = []

    # Create the Fibonacci sequence of points
    for i in range(N):
        # Calculate the position of the point in the x direction
        x = minx + (i % nx) * average_distance

        # Calculate the position of the point in the y direction
        y = miny + ((i * golden_ratio) % ny) * average_distance

        # Create a point at the calculated position
        point = Point(x, y)

        # If the point is inside the polygon, add it to the list
        if polygon.contains(point):
            points.append(point)

    # Return the list of points
    return points
