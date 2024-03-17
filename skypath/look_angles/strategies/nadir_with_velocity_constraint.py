import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation

from skypath.attitude._base_attitude_strategy import BaseAttitudeStrategy
from skypath.attitude.identity import Identity
from skypath.coordinates import CoordinateInterpolator
from skypath.look_angles.look_angle_funcs import get_look_angles_to
from skypath.look_angles.strategies._base_look_angle_stragtegy import (
    BaseLookAngleStrategy,
)


class NadirWithVelocityConstraint(BaseLookAngleStrategy):
    def __init__(
        self,
        frame: str = "gcrs",
        attitude: BaseAttitudeStrategy = Identity(),
        check_orthogonal: bool = True,
    ) -> None:
        """
        Define the look angles as Z facing nadir, Y as the cross of X and Y, and X in the direction of velocity.

        Args:
            frame (str, optional): Coordinate frame to compute the look angles in. NOTE: Your attitude must also be defined in this coordinate frame.. Defaults to 'gcrs'.
            attitude (BaseAttitudeStrategy, optional): Strategy for getting attitude of the body at a given time. Defaults to Identity().
            check_orthogonal (bool, optional): Checks orthogonal results are produced from the rotation matrices. Defaults to True.
        """
        super().__init__()
        self.frame = frame
        self.attitude = attitude
        self.check_orthogonal = check_orthogonal

    def __call__(
        self,
        observer: CoordinateInterpolator,
        target: CoordinateInterpolator,
        time: Time,
    ) -> tuple:
        observer_rotation = self.get_rotation(observer, time)
        attitude_rotation = self.attitude.get_rotation(observer, time)
        observer_orientation = attitude_rotation * observer_rotation

        observer_position = observer.state_at(time, self.frame).cartesian.xyz.to(u.m)
        target_position = target.state_at(time, self.frame).cartesian.xyz.to(u.m)
        return get_look_angles_to(
            observer_position.T, target_position.T, observer_orientation
        )

    def get_rotation(self, coord: CoordinateInterpolator, time: Time) -> Rotation:
        # stack the vectors into a uniform format
        state = coord.state_at(time, self.frame)
        pos = state.cartesian.xyz.T
        differentials = state.cartesian.differentials.get("s", None)
        if differentials is None:
            vel = np.zeros_like(pos).value * (u.m / u.s)
        else:
            vel = differentials.d_xyz.T.to(u.m / u.s)
        assert pos.ndim == vel.ndim

        if pos.ndim == 1:
            pos = np.reshape(pos, (-1, 3))
            vel = np.reshape(vel, (-1, 3))

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

        if self.check_orthogonal:
            # Verify that each axis is orthogonal
            assert np.allclose(np.cross(x_axis, y_axis), z_axis, atol=1e-6)
            assert np.allclose(np.cross(z_axis, x_axis), y_axis, atol=1e-6)
            assert np.allclose(np.cross(y_axis, z_axis), x_axis, atol=1e-6)

        return Rotation.from_matrix(rot_matrix)
