from abc import ABC, abstractmethod

from astropy.time import Time
from scipy.spatial.transform import Rotation


class BaseAttitudeStrategy(ABC):
    @abstractmethod
    def at(self, time: Time) -> Rotation:
        """
        Return the attitude Rotation instance at the specific time(s).
        By default, the attitude class has no reference frame attached to it,
        so it is up to the user to understand what reference frame their attitudes
        utilize.

        Args:
            time (Time): Times to retrieve attitude data for.

        Returns:
            Rotation: scipy.spatial.transform.Rotation instance defining the attitude.
        """
        pass
