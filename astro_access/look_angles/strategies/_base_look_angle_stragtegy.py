from abc import ABC, abstractclassmethod
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time



class BaseLookAngleStrategy(ABC):
    @abstractclassmethod
    def __call__(self, observer: CoordinateInterpolator, target: CoordinateInterpolator, time: Time) -> tuple:
        pass