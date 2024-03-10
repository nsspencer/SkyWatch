from astro_access.constraints._base_constraint import BaseAccessConstraint
from astro_access.frame_interpolator import FrameInterpolator
from astropy.time import Time
import astropy.units as u
import numpy as np


class Temporal(BaseAccessConstraint):
    def __init__(self, from_time: Time, to_time: Time, inner: bool = True) -> None:
        super().__init__()
        self.from_time = from_time
        self.to_time = to_time
        self.inner = inner
    
    def __call__(self, observer: FrameInterpolator, target: FrameInterpolator, time: Time, bounds_check: bool = True) -> np.ndarray:
        if self.inner:
            return (time >= self.from_time) & (time <= self.to_time)
        return (time <= self.from_time) | (time >= self.to_time)