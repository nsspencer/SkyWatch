from abc import ABC, abstractmethod
from typing import Any

import astropy.units as u
from astropy.time import Time

from skypath import SkyPath


class TerrainSlope:
    def __init__(self) -> None:
        pass


class AttitudeStrategy:
    pass


class Fixed(AttitudeStrategy):
    pass


class Slerped(AttitudeStrategy):
    pass


class VelocityWithNadirConstraint(AttitudeStrategy):
    def __init__(self, observer: SkyPath, frame: str = "gcrs") -> None:
        super().__init__()


class LookAnglesStrategy:
    pass


class LocalTangent(LookAnglesStrategy):
    def __init__(self, terrain_slope: TerrainSlope = None) -> None:
        super().__init__()
        self.terrain_slope = terrain_slope


class BodyFrame(LookAnglesStrategy):
    pass


class LookAnglesMixin:
    @classmethod
    def local_tangential(cls, terrain_slope: TerrainSlope = None):
        return LocalTangent(terrain_slope)

    @classmethod
    def body_frame(cls, *args):
        return BodyFrame(*args)


class LookAngles(LookAnglesMixin, LookAnglesStrategy):
    def __call__(self, time):
        return


class AccessStrategy:
    pass


class Spatial(AccessStrategy):
    def __init__(
        self,
        observer,
        target,
        constraints,
        precise_endpoints: bool = True,
        precision: u.m = 0.1 * u.m,
    ) -> None:
        self.observer = observer
        self.target = target
        self.constraints = constraints
        self.precise_endpoints = precise_endpoints
        self.precision = precision


class AccessMixin:
    @classmethod
    def spatial(
        cls,
        observer,
        target,
        constraints,
        precise_endpoints: bool = True,
        precision: u.m = 0.1 * u.m,
    ):
        return Spatial(observer, target, constraints, precise_endpoints, precision)


class Access(AccessMixin):
    pass


def usage_test():
    sat1 = SkyPath(...)
    ground_station1 = SkyPath(...)

    gs1_look_angles_fn = LocalTangent(
        observer=ground_station1, terrain_slope=None, frame="itrs"
    )
    sat1_look_angles_fn = BodyFrame(
        observer=sat1, frame="gcrs", attitude_strategy=VelocityWithNadirConstraint()
    )

    times = Time.now()
    gs1_to_sat1_access = Access.spatial(
        observer=ground_station1,
        target=sat1,
        time=times,
        constraints=[],
        precise_endpoints=True,
        precision=0.1 * u.s,
    )


class AER:
    pass


class LookAngleStrategy(ABC):
    @abstractmethod
    def __call__(self, target: "SkyPath2", time: Time) -> AER:
        pass


class LocalTangential(LookAnglesStrategy):
    def __init__(self, terrain_slope: TerrainSlope = None) -> None:
        super().__init__()
        self.terrain_slope = terrain_slope

    def __call__(self, target: "SkyPath2", time: Time) -> AER:
        pass


class BodyFrame(LookAnglesStrategy):
    def __init__(
        self, attitude: AttitudeStrategy, coordinate_frame: str = "gcrs"
    ) -> None:
        super().__init__()

    def __call__(self, target: "SkyPath2", time: Time) -> AER:
        pass


class SkyPath2:
    def look_angles_to(
        self, target: "SkyPath2", time: Time, strategy: LookAngleStrategy
    ) -> AER:
        return strategy(target, time)


class AccessStrategy(ABC):
    @abstractmethod
    def __call__(self, target: "SkyPath2", time: Time) -> Access:
        pass


if __name__ == "__main__":
    pass
