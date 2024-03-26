import numpy as np
from astropy.time import Time

from skywatch.access.base_constraint import BaseAccessConstraint


class Temporal(BaseAccessConstraint):
    def __init__(
        self, min_time: Time = None, max_time: Time = None, inner: bool = True
    ) -> None:
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
        self.min_time = min_time
        self.max_time = max_time
        self.inner = inner

        if self.min_time is None and self.max_time is None:
            raise ValueError(
                "Min time and Max time cannot both be None. You must set one, or both."
            )

    def __call__(self, time: Time, *args, **kwargs) -> np.ndarray:
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
        if self.min_time is not None and self.max_time is not None:
            if self.inner:
                return (time >= self.min_time) & (time <= self.max_time)
            else:
                return (time <= self.min_time) | (time >= self.max_time)

        if self.min_time is not None:
            return time >= self.min_time

        if self.max_time is not None:
            return time <= self.max_time
