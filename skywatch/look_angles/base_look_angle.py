from abc import ABC, abstractmethod

from astropy.time import Time

from skywatch.look_angles.aert import AzElRangeTime
from skywatch.skypath import SkyPath


class BaseLookAngleStrategy(ABC):
    @abstractmethod
    def get_look_angles(self, target: SkyPath, time: Time) -> AzElRangeTime:
        """
        Get the Azimuth, Elevation, Range, and Time to a target at the given time(s).

        Args:
            target (SkyPath): SkyPath object to calculate look angles to.
            time (Time): Time(s) to calculate look angles to the target during.

        Returns:
            AzElRangeTime: Look angles to the target during the provided times.
        """
        pass
