from astro_access.access.constraints._base_constraint import BaseAccessConstraint
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
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
            min_az (u.deg, optional): minimum allowable azimuth angle to target. Defaults to 0*u.deg.
            max_az (u.deg, optional): maximum allowable azimuth angle to target. Defaults to 360*u.deg.
            min_el (u.deg, optional): minimum allowable elevation angle to target. Defaults to 0*u.deg.
            max_el (u.deg, optional): maximum allowable elevation angle to target. Defaults to 90*u.deg.
            min_range (u.m, optional): minimum allowable range angle to target. Defaults to 0*u.m.
            max_range (u.m, optional): maximum allowable range to target. Defaults to np.inf*u.m.
        """
        super().__init__()
        self.min_az = min_az
        self.max_az = max_az
        self.min_el = min_el
        self.max_el = max_el
        self.min_range = min_range
        self.max_range = max_range
    
    def __call__(self, observer: CoordinateInterpolator, target: CoordinateInterpolator, time: Time, bounds_check: bool = True) -> np.ndarray:
        """Constraints the times when the azimuth, elevation, and range values are satisfied from the observer to the target.

        Args:
            observer (CoordinateInterpolator): _description_
            target (CoordinateInterpolator): _description_
            time (Time): _description_
            bounds_check (bool, optional): Whether or not to check the provided times are within the bounds 
            of the observer and target CoordinateInterpolator.

        Returns:
            np.ndarray: boolean array representing times when this constraint succeeds vs. fails.
        """
        
        # NOTE: Using pymap3d is tremendously faster than having astropy calculate the AltAz frame from the observer to the target,
        # however there are cases for distant objects that this loses accuracy.
        az, el, rng  = pymap3d.ecef2aer(*target.state_at(time, 'itrs', bounds_check=bounds_check).cartesian.xyz.to(u.m).value,\
            *pymap3d.ecef2geodetic(*observer.state_at(time, 'itrs', bounds_check=bounds_check).cartesian.xyz.to(u.m).value))
        return (az >= self.min_az.value) & (az <= self.max_az.value) & (el >= self.min_el.value) & (el <= self.max_el.value) & (rng >= self.min_range.value) & (rng <= self.max_range.value)
