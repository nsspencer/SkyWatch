from typing import Union

from astropy.time import Time

from skywatch.attitude import BaseAttitudeStrategy
from skywatch.look_angles import BaseLookAngleStrategy
from skywatch.skypath import SkyPath


class Platform:
    def __init__(
        self,
        coordinate: SkyPath,
        look_angles_strategy: BaseLookAngleStrategy,
        attitude: BaseAttitudeStrategy = None,
    ) -> None:
        self.coordinate = coordinate
        self.attitude = attitude
        self.look_angle_strategy = look_angles_strategy

    def look_angles(
        self,
        time: Time,
        target: Union["Platform", SkyPath],
    ) -> tuple:
        if isinstance(target, Platform):
            _target_coordinate = target.coordinate
        elif isinstance(target, SkyPath):
            _target_coordinate = target
        else:
            raise TypeError(
                "Target needs to be an instance of a Platform or a SkyPath."
            )
        return self.look_angle_strategy.calculate(
            time, self.coordinate, _target_coordinate, self.attitude
        )
