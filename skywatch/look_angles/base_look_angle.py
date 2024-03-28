from abc import ABC, abstractmethod

from astropy.time import Time

from skywatch.look_angles.aert import AzElRangeTime
from skywatch.skypath import SkyPath


class BaseLookAngleStrategy(ABC):
    @abstractmethod
    def get_look_angles(self, target: SkyPath, time: Time) -> AzElRangeTime:
        pass
