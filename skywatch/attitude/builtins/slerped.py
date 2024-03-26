from astropy.time import Time
from scipy.spatial.transform import Rotation, Slerp

from skywatch.attitude.base_attitude import BaseAttitudeStrategy


class Slerped(BaseAttitudeStrategy):
    def __init__(self, times: Time, rotations: Rotation) -> None:
        """
        Using known attitudes (usually from real observations), do
        spherical linear interpolation between them.

        Args:
            times (Time): Times of known observations
            rotations (Rotation): Known attitudes
        """
        self.slerper = Slerp(times, rotations)

    def at(self, time: Time) -> Rotation:
        return self.slerper(time)
