from abc import ABC, abstractmethod

from astropy.time import Time

from skywatch.coordinates import SkyPath
from skywatch.look_angles.aer import AzElRangeTime


class BaseLookAngleStrategy(ABC):
    @abstractmethod
    def __call__(self, time: Time, target: SkyPath, *args, **kwargs) -> AzElRangeTime:
        pass
