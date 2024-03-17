import math
from typing import List


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
