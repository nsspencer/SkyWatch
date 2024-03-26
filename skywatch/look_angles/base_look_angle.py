from abc import ABC, abstractmethod

from astropy.time import Time

from skywatch.attitude import BaseAttitudeStrategy
from skywatch.look_angles.aer import AzElRangeTime
from skywatch.skypath import SkyPath


class BaseLookAngleStrategy(ABC):
    @abstractmethod
    def calculate(
        self,
        time: Time,
        target: SkyPath,
        observer: SkyPath,
        observer_attitude: BaseAttitudeStrategy,
    ) -> AzElRangeTime:
        pass
