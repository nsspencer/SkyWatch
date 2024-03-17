from typing import Any

from astropy.time import Time

from skywatch.coordinates.skypath import SkyPath
from skywatch.look_angles.aer import AzElRangeTime
from skywatch.look_angles.creation import LookAnglesCreationMixin
from skywatch.look_angles.strategies._base_strategy import BaseLookAngleStrategy


class LookAngles(LookAnglesCreationMixin):
    def __init__(self, strategy: BaseLookAngleStrategy) -> None:
        self.strategy = strategy

    def __call__(
        self, time: Time, target: SkyPath, *args: Any, **kwargs: Any
    ) -> AzElRangeTime:
        return self.strategy(time, target, *args, **kwargs)

    def calculate(self, time: Time, target: SkyPath, *args, **kwargs) -> AzElRangeTime:
        return self.__call__(time, target, *args, **kwargs)
