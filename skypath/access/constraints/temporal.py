from skypath.access.constraints._base_constraint import BaseAccessConstraint
from skypath.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time
import numpy as np


class Temporal(BaseAccessConstraint):
    def __init__(self, from_time: Time, to_time: Time, inner: bool = True) -> None:
        """
        Only allows times in the time array that are within self.from_time and self.to_time (if self.inner is True)
        to pass the constraint.
        
        If self.inner is False, only times that do not fall between the self.from_time and self.to_time are allowed
        to pass this constraint.

        Args:
            from_time (Time): Lower time bound.
            to_time (Time): Upper time bound.
            inner (bool, optional): If True, only times within the lower time bound and upper time bound will pass
            this constraint. If False, only times outside the time bounds will pass. Defaults to True.
        """
        super().__init__()
        self.from_time = from_time
        self.to_time = to_time
        self.inner = inner
    
    def __call__(self, observer: CoordinateInterpolator, target: CoordinateInterpolator, time: Time, bounds_check: bool = True) -> np.ndarray:
        """
        Compute times that pass this constraint.

        Args:
            observer (CoordinateInterpolator): unused, can be None.
            target (CoordinateInterpolator): unused, can be None.
            time (Time): times to check against this constraints allowed times.
            bounds_check (bool, optional): unused. Defaults to True.

        Returns:
            np.ndarray: Boolean array representing times that pass this constraint.
        """
        if self.inner:
            return (time >= self.from_time) & (time <= self.to_time)
        return (time <= self.from_time) | (time >= self.to_time)