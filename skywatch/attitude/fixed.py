from astropy.time import Time
from scipy.spatial.transform import Rotation

from skywatch.attitude._base_attitude import BaseAttitudeStrategy


class Fixed(BaseAttitudeStrategy):
    def __init__(self, rotation: Rotation, frame: str = "gcrs") -> None:
        super().__init__(frame)
        self.rotation = rotation

    def __call__(self, time: Time) -> Rotation:
        return self.rotation
