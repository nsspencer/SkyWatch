from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation, Slerp

from skywatch.attitude._base_attitude import BaseAttitudeStrategy


class Slerped(BaseAttitudeStrategy):
    def __init__(self, times: Time, rotations: Rotation, frame: str = "gcrs") -> None:
        super().__init__(frame)
        self.slerper = Slerp(times, rotations)

    def __call__(self, time: Time) -> Rotation:
        return self.slerper(time)
