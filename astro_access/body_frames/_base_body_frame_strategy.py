
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
from abc import ABC, abstractclassmethod
from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation


class BaseBodyFrameStrategy(ABC):
    
    @abstractclassmethod
    def get_rotation(self, coord: CoordinateInterpolator, time: Time) -> Rotation:
        pass