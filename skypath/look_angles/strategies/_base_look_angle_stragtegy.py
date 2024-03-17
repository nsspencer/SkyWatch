from abc import ABC, abstractclassmethod

from astropy.time import Time

from skypath.coordinates import CoordinateInterpolator


class BaseLookAngleStrategy(ABC):
    @abstractclassmethod
    def __call__(
        self,
        observer: CoordinateInterpolator,
        target: CoordinateInterpolator,
        time: Time,
    ) -> tuple:
        pass
