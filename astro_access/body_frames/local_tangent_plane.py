from astro_access.body_frames._base_body_frame_strategy import BaseBodyFrameStrategy
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time
import numpy as np
from scipy.spatial.transform.rotation import Rotation
import astropy.units as u


class LocalTangentPlane(BaseBodyFrameStrategy):
    def __init__(self, frame: str = 'itrs') -> None:
        super().__init__()
        self.frame = frame
        
    def get_rotation(self, coord: CoordinateInterpolator, time: Time) -> Rotation:
        # stack the vectors into a uniform format
        state = coord.state_at(time, self.frame)
        
        lon, lat, alt = state.earth_location.to_geodetic()
        
        # Extract observer's position
        observer_lon = lon
        observer_lat = lat

        # Calculate the local tangent plane's orientation
        north = np.array([-np.sin(observer_lat) * np.cos(observer_lon),
                          -np.sin(observer_lat) * np.sin(observer_lon),
                          np.cos(observer_lat)])

        east = np.array([-np.sin(observer_lon),
                         np.cos(observer_lon),
                         np.zeros_like(observer_lon)])

        up = np.array([np.cos(observer_lat) * np.cos(observer_lon),
                       np.cos(observer_lat) * np.sin(observer_lon),
                       np.sin(observer_lat)])

        # Normalize the basis vectors
        east /= np.linalg.norm(east, axis=0)
        north /= np.linalg.norm(north, axis=0)
        up /= np.linalg.norm(up, axis=0)

        # Create the rotation matrix from ENU to the local tangent plane
        rot_matrix = np.dstack((east.T, north.T, up.T))
        
        return Rotation.from_matrix(rot_matrix)