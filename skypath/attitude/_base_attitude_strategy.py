from skypath.coordinates.coordinate_interpolator import CoordinateInterpolator

from astropy.time import Time
from abc import ABC, abstractclassmethod
from scipy.spatial.transform.rotation import Rotation


class BaseAttitudeStrategy(ABC):
    @abstractclassmethod
    def get_rotation(self, coordinate: CoordinateInterpolator, time: Time) -> Rotation:
        pass