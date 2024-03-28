import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude.base_attitude import BaseAttitudeStrategy
from skywatch.look_angles.aert import AzElRangeTime
from skywatch.look_angles.base_look_angle import BaseLookAngleStrategy
from skywatch.skypath.skypath import SkyPath


class XForwardZNadir(BaseLookAngleStrategy):
    def __init__(
        self,
        observer: SkyPath,
        attitude_strategy: BaseAttitudeStrategy,
        frame: str = "gcrs",
    ) -> None:
        super().__init__()
        self.observer = observer
        self.attitude_strategy = attitude_strategy
        self.frame = frame

    def get_look_angles(
        self,
        target: SkyPath,
        time: Time,
    ) -> AzElRangeTime:
        observer_position = (
            self.observer.state_at(time, self.frame).cartesian.xyz.to(u.m).value
        )
        target_position = target.state_at(time, self.frame).cartesian.xyz.to(u.m).value
        attitude = self.attitude_strategy.at(time)
        az, el, rng = XForwardZNadir.get_look_angles_to(
            observer_position.T, target_position.T, attitude
        )

        return AzElRangeTime(az * u.deg, el * u.deg, rng * u.m, time)

    @staticmethod
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
            0 degrees faces the center of the earth (directly below the satellite)
            90 degrees is in the direction of the velocity
            180 degrees faces the opposite direction of the center of the earth (directly above the satellite)

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
        azimuth = XForwardZNadir._counterclockwise_angle_between(
            x_axis, observer_to_target_projected[:, :2]
        )

        # calculate the elevation angle
        elevation = XForwardZNadir._elevation_between(
            observer_to_target, component_along_z
        )

        # calculate the range
        rng = np.linalg.norm(observer_to_target, axis=1)

        return azimuth, elevation, rng

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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
