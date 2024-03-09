from abc import ABC, abstractclassmethod
from astropy.time import Time
import numpy as np

from astro_access.frame_interpolator import FrameInterpolator


class BaseAccessConstraint(ABC):
    @abstractclassmethod
    def __call__(self, observer: FrameInterpolator, target: FrameInterpolator, time: Time) -> np.ndarray:
        """
        Return truth array of times that satisfy the access constraint between the observer and target.

        Args:
            observer (FrameInterpolator): Coordinate to get access from
            target (FrameInterpolator): Coordinate to get access to
            time (Time): Time range of the access calculation

        Returns:
            np.ndarray: Truth array of times that satisfy the constraint.
        """
        pass
    
