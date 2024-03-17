from abc import ABC, abstractclassmethod

from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation

from skypath.coordinates import CoordinateInterpolator


class BaseAttitudeStrategy(ABC):
    @abstractclassmethod
    def get_rotation(self, coordinate: CoordinateInterpolator, time: Time) -> Rotation:
        pass
