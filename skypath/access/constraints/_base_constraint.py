from abc import ABC, abstractclassmethod

import numpy as np
from astropy.time import Time

from skypath.coordinates import CoordinateInterpolator


class BaseAccessConstraint(ABC):
    @abstractclassmethod
    def __call__(
        self,
        observer: CoordinateInterpolator,
        target: CoordinateInterpolator,
        time: Time,
        bounds_check: bool = True,
    ) -> np.ndarray:
        """
        Return truth array of times that satisfy the access constraint between the observer and target.

        Args:
            observer (CoordinateInterpolator): Coordinate to get access from
            target (CoordinateInterpolator): Coordinate to get access to
            time (Time): Time range of the access calculation
            bounds_check (bool, optional): Whether or not to check the provided times are within the bounds
            of the observer and target CoordinateInterpolator.

        Returns:
            np.ndarray: Truth array of times that satisfy the constraint.
        """
        pass
