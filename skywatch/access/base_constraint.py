from abc import ABC, abstractmethod

import numpy as np
from astropy.time import Time


class BaseAccessConstraint(ABC):
    def __init__(self) -> None:
        """
        Represents an access constraint which must be satisfied when provided to
        and access calcualtion.
        """
        super().__init__()

    @abstractmethod
    def __call__(self, time: Time) -> np.ndarray:
        """
        Method called by the Access algorithm which passes in a

        Args:
            time (Time): Astropy Time(s) to check this access constraint during.

        Returns:
            np.ndarray: boolean array matching the shape of the Time argument.
            True means access was successful at that time, False means access
            failed at that time.
        """
        pass
