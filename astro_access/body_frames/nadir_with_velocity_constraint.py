from astro_access.body_frames._base_body_frame_strategy import BaseBodyFrameStrategy
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time
import numpy as np
from scipy.spatial.transform.rotation import Rotation
from typing import Union
import astropy.units as u

class NadirWithVelocityConstraint(BaseBodyFrameStrategy):
    def __init__(self, frame: str = 'gcrs', check_orthogonal:bool = True) -> None:
        super().__init__()
        self.frame = frame
        self.check_orthogonal = check_orthogonal
        
    def get_rotation(self, coord: CoordinateInterpolator, time: Time) -> Rotation:
        # stack the vectors into a uniform format
        state = coord.state_at(time, self.frame)
        pos = state.cartesian.xyz.T
        differentials = state.cartesian.differentials.get('s', None)
        if differentials is None:
            vel = np.zeros_like(pos).value*(u.m/u.s)
        else:
            vel = differentials.d_xyz.T.to(u.m/u.s)
        assert pos.ndim == vel.ndim
        
        if pos.ndim == 1:
            pos = np.reshape(pos, (-1,3))
            vel = np.reshape(vel, (-1,3))
        
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
        rot_matrix = np.transpose(np.array([x_axis, y_axis, z_axis]), axes=(1,2,0))

        if self.check_orthogonal:
            # Verify that each axis is orthogonal
            assert np.allclose(np.cross(x_axis, y_axis), z_axis, atol=1e-6)
            assert np.allclose(np.cross(z_axis, x_axis), y_axis, atol=1e-6)
            assert np.allclose(np.cross(y_axis, z_axis), x_axis, atol=1e-6)
            
        return Rotation.from_matrix(rot_matrix)