import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation

from skywatch.attitude._base_attitude import BaseAttitudeStrategy
from skywatch.coordinates.skypath import SkyPath
from skywatch.look_angles import funcs
from skywatch.look_angles.aer import AzElRangeTime
from skywatch.look_angles.strategies._base_strategy import BaseLookAngleStrategy


class BodyFrame(BaseLookAngleStrategy):
    def __init__(
        self,
        observer: SkyPath,
        attitude_strategy: BaseAttitudeStrategy,
        attitude_offset: Rotation = None,
    ) -> None:
        super().__init__()
        self.observer = observer
        self.attitude_strategy = attitude_strategy
        self.attitude_offset = attitude_offset

    def __call__(self, time: Time, target: SkyPath, *args, **kwargs) -> AzElRangeTime:
        observer_position = (
            self.observer.state_at(time, self.attitude_strategy.frame)
            .cartesian.xyz.to(u.m)
            .value
        )
        target_position = (
            target.state_at(time, self.attitude_strategy.frame)
            .cartesian.xyz.to(u.m)
            .value
        )

        observer_attitude = self.attitude_strategy(time)
        if self.attitude_offset is not None:
            observer_attitude = observer_attitude * self.attitude_offset

        az, el, rng = BodyFrame.get_look_angles_to(
            observer_position.T, target_position.T, observer_attitude
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
        azimuth = funcs._counterclockwise_angle_between(
            x_axis, observer_to_target_projected[:, :2]
        )

        # calculate the elevation angle
        elevation = funcs._elevation_between(observer_to_target, component_along_z)

        # calculate the range
        rng = np.linalg.norm(observer_to_target, axis=1)

        return azimuth, elevation, rng
