from skypath.attitude._base_attitude_strategy import BaseAttitudeStrategy
from skypath.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation
import numpy as np



class Identity(BaseAttitudeStrategy):
    def __init__(self) -> None:
        super().__init__()
        
    def get_rotation(self, coordinate: CoordinateInterpolator, time: Time) -> Rotation:
        return Rotation.from_matrix(np.array(np.eye(3)))