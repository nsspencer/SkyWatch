from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude.base_attitude import BaseAttitudeStrategy


class Fixed(BaseAttitudeStrategy):
    def __init__(self, rotation: Rotation) -> None:
        """
        Provide a fixed single attitude as the attitude for all time durations.

        Args:
            rotation (Rotation): Rotation to return at all time steps.
        """
        self.rotation = rotation

    def at(self, time: Time) -> Rotation:
        return self.rotation
