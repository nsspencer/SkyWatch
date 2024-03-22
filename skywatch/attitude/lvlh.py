import astropy.units as u
import numpy as np
from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude._base_attitude import BaseAttitudeStrategy
from skywatch.coordinates.skypath import SkyPath


class LVLH(BaseAttitudeStrategy):
    def __init__(self, observer: SkyPath, frame: str = "gcrs") -> None:
        super().__init__(frame)
        self.observer = observer

    def __call__(self, time: Time) -> Rotation:
        state = self.observer.state_at(time, frame=self.frame)
        pos = state.cartesian.xyz.to(u.m).value
        differentials = state.cartesian.differentials.get("s", None)
        if differentials is None:
            raise AttributeError(
                "Observer has no velocity component. Velocity is needed to calculate the LVLH frame."
            )
        vel = differentials.d_xyz.to(u.m / u.s).value

        # Now you have the
        return Rotation.from_matrix(LVLH.calculate_reference_frame(*pos, *vel))

    @staticmethod
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
