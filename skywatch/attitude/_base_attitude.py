from abc import ABC, abstractmethod

from astropy.time import Time
from scipy.spatial.transform.rotation import Rotation


class BaseAttitudeStrategy(ABC):
    def __init__(self, frame: str = "gcrs") -> None:
        super().__init__()
        self.frame = frame

    @abstractmethod
    def __call__(self, time: Time) -> Rotation:
        pass
