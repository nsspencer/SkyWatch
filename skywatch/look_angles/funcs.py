import astropy.units as u
import numpy as np


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


@u.quantity_input(latitude=u.deg, longitude=u.deg)
def calculate_north(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    north_direction = np.array(
        [
            -np.sin(latitude) * np.cos(longitude),
            -np.sin(latitude) * np.sin(longitude),
            np.cos(latitude),
        ]
    )
    return north_direction


def _elevation_between(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Gets the elevation angle (in degrees) between the two vectors.

    Args:
        v1 (np.ndarray): vector 1
        v2 (np.ndarray): vector 2

    Returns:
        np.ndarray: elevation angle in degrees
    """
    # Normalize the vectors to unit vectors
    unit_vector1 = v1 / np.linalg.norm(v1, axis=1)[:, np.newaxis]
    unit_vector2 = v2 / np.linalg.norm(v2, axis=1)[:, np.newaxis]

    # Compute the dot product
    dot_product = np.sum((unit_vector1 * unit_vector2), axis=1)

    # Calculate the angle (in radians) & convert to degrees
    elevation = np.degrees(np.arccos(dot_product))
    return elevation


def _clockwise_angle_between(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Calculates the unambiguous clockwise angle between two vectors.

    Args:
        v1 (np.ndarray): vectors 1
        v2 (np.ndarray): vectors 2

    Returns:
        np.ndarray: angles in degrees (between 0 and 360)
    """
    angle_1 = np.arctan2(v1[:, 1], v1[:, 0]) * 180 / np.pi
    angle_2 = np.arctan2(v2[:, 1], v2[:, 0]) * 180 / np.pi
    angle = angle_2 - angle_1
    angle = np.where(angle < 0, angle + 360, angle)
    return angle


def _counterclockwise_angle_between(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Calculates the unambiguous counterclockwise angle between two vectors.

    Args:
        v1 (np.ndarray): vectors 1
        v2 (np.ndarray): vectors 2

    Returns:
        np.ndarray: angles in degrees (between 0 and 360)
    """
    angle_1 = np.arctan2(v1[:, 1], v1[:, 0]) * 180 / np.pi
    angle_2 = np.arctan2(v2[:, 1], v2[:, 0]) * 180 / np.pi
    angle = angle_1 - angle_2  # Swapped the order here to make counterclockwise
    angle = np.where(angle < 0, angle + 360, angle)
    return angle
