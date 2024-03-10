from astro_access.constraints._base_constraint import BaseAccessConstraint
from astro_access.frame_interpolator import FrameInterpolator
import astropy.units as u
from astropy.time import Time
import numpy as np
import pymap3d


class AzElRange(BaseAccessConstraint):
    @u.quantity_input(min_az = u.deg, max_az = u.deg,
                      min_el = u.deg, max_el = u.deg,
                      min_range = u.m, max_range = u.m)
    def __init__(self,
                 min_az: u.deg = 0 * u.deg, max_az: u.deg = 360 * u.deg,
                 min_el: u.deg = 0 * u.deg, max_el: u.deg = 90 * u.deg,
                 min_range: u.m = 0 * u.m, max_range: u.m = np.inf * u.m) -> None:
        """
        Determine access using ENU azimuth, elevation, and range values from an earth based position.

        Args:
            min_az (u.deg, optional): _description_. Defaults to 0*u.deg.
            max_az (u.deg, optional): _description_. Defaults to 360*u.deg.
            min_el (u.deg, optional): _description_. Defaults to 0*u.deg.
            max_el (u.deg, optional): _description_. Defaults to 90*u.deg.
            min_range (u.m, optional): _description_. Defaults to 0*u.m.
            max_range (u.m, optional): _description_. Defaults to np.inf*u.m.
        """
        super().__init__()
        self.min_az = min_az
        self.max_az = max_az
        self.min_el = min_el
        self.max_el = max_el
        self.min_range = min_range
        self.max_range = max_range
    
    def __call__(self, observer: FrameInterpolator, target: FrameInterpolator, time: Time) -> np.ndarray:
        az, el, rng  = pymap3d.ecef2aer(*observer.state_at(time, 'itrs').cartesian.xyz.to(u.m).value, *pymap3d.ecef2geodetic(*target.state_at(time, 'itrs').cartesian.xyz.to(u.m).value))
        az = az * u.deg
        el = el * u.deg
        rng = rng * u.m
        return (az >= self.min_az) & (az <= self.max_az) & (el >= self.min_el) & (el <= self.max_el) & (rng >= self.min_range) & (rng <= self.max_range)
