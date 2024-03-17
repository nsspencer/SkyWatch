import astropy.units as u
import numpy as np
from scipy.spatial.transform.rotation import Rotation


@u.quantity_input(altitude=u.km, body_radius=u.km)
def max_elevation_angle(altitude: u.km, body_radius: u.km = 6371 * u.km):
    """Calculates the maximum elevation angle from NADIR
    a satellite can point to and still intersect the earth.

    Args:
        h (height): _description_

    Returns:
        _type_: _description_
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


def calculate_reference_frame(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    v_x: np.ndarray,
    v_y: np.ndarray,
    v_z: np.ndarray,
    check_validity: bool = True,
    abs_tolerance: float = 1.0e-9,
) -> np.ndarray:
    """
    Returns a matrix representing the X,Y,Z transform of the reference frame.
    This reference frame is defined with Z being NADIR alignment, X facing
    the direction of the velocity, and Y as the cross product between Z and X
    complete the right handed system.

    The final matrix looks like such:
    [X1,Y1,Z1]
    [X2,Y2,Z2]
    [X3,Y3,Z3]

    Returns an array of matrices representing the reference frames for the given
    positions and velocities.

    Args:
        x (np.ndarray): position X (in any coordinate system)
        y (np.ndarray): position Y (in any coordinate system)
        z (np.ndarray): position Z (in any coordinate system)
        v_x (np.ndarray): velocity X (in any coordinate system)
        v_y (np.ndarray): velocity Y (in any coordinate system)
        v_z (np.ndarray): velocity Z (in any coordinate system)
        check_validity (bool, optional): assert that the dot product of the resulting rotation matrices are identities. Defaults to True.
        abs_tolerance (float, optional): tolerance value for the identity matrix validity check. Defaults to 1.e-2.

    Returns:
        np.ndarray: rotation matrix defining the reference frame.
    """

    # stack the vectors into a uniform format
    pos = np.column_stack([x, y, z])
    vel = np.column_stack([v_x, v_y, v_z])
    assert pos.ndim == vel.ndim

    r_norms = np.linalg.norm(pos, axis=1)
    v_norms = np.linalg.norm(vel, axis=1)
    z_axis = -pos / r_norms[:, np.newaxis]
    x_axis = vel / v_norms[:, np.newaxis]

    # Subtract the projection from the x-axis to make it orthogonal to the Z direction
    x_axis -= np.sum(x_axis * z_axis, axis=1)[:, np.newaxis] * z_axis
    x_axis /= np.linalg.norm(x_axis, axis=1)[:, np.newaxis]

    # define the y axis
    y_axis = np.cross(z_axis, x_axis, axis=1)

    # transpose the rotation matrix to make it 3x3xN
    rot_matrix = np.transpose(np.array([x_axis, y_axis, z_axis]), axes=(1, 2, 0))

    if check_validity:
        # Verify that each axis is orthogonal
        assert np.allclose(np.cross(x_axis, y_axis), z_axis, atol=abs_tolerance)
        assert np.allclose(np.cross(z_axis, x_axis), y_axis, atol=abs_tolerance)
        assert np.allclose(np.cross(y_axis, z_axis), x_axis, atol=abs_tolerance)

    return rot_matrix


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
    return angle * u.deg


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
    return angle * u.deg


def get_look_angles_to(
    observer_pos: np.ndarray, target_pos: np.ndarray, observer_frame: Rotation
) -> tuple:
    """Calculates the look angles from a satellites point of view to a geodedic location.

    Azimuth is defined as degrees clockwise.
        0 degrees faces the direction of velocity.
        90 degrees faces right of the direction of velocity.
        180 degrees faces opposite direction of velocity.
        270 degrees faces left of the direction of velocity.

    Elevation is defined as degrees of rise from the center of the earth.
        -90 degrees faces the center of the earth (directly below the satellite)
        0 degrees is in the direction of the velocity
        90 degrees faces the opposite direction of the center of the earth (directly above the satellite)

    Range is defined as straight line (euclidean) distance from the satellite to the target in kilometers.

    Returns:
        tuple: azimuth, elevation, range
    """
    assert (
        observer_pos.ndim == target_pos.ndim
    ), "Dimensions of observer and target positions must be equivalent."
    if observer_pos.ndim == 1:
        observer_pos = np.reshape(observer_pos, (-1, 3))
        target_pos = np.reshape(target_pos, (-1, 3))

    # get the axis we need from the local frame
    _obs_frame_matrix = observer_frame.as_matrix()
    x_axis = _obs_frame_matrix[:, :, 0]
    z_axis = _obs_frame_matrix[:, :, 2]

    # Transform the target position to the local reference frame
    observer_to_target = target_pos - observer_pos

    # Project the satellite_to_target vector onto the plane perpendicular to the z-axis
    component_along_z = (
        np.sum(observer_to_target * z_axis, axis=1)[:, np.newaxis] * z_axis
    )
    observer_to_target_projected = observer_to_target - component_along_z

    # calculate the azimuth from the x axis
    azimuth = _counterclockwise_angle_between(
        x_axis, observer_to_target_projected[:, :2].value
    )

    # calculate the elevation angle
    elevation = _elevation_between(observer_to_target, component_along_z)

    # calculate the range
    rng = np.sqrt(np.sum(observer_to_target**2, axis=1)).to(u.m)

    return azimuth, elevation, rng
