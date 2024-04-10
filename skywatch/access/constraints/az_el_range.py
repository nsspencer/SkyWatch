import astropy.units as u
import numpy as np
import pymap3d
from astropy.time import Time

from skywatch.access.base_constraint import BaseAccessConstraint
from skywatch.skypath import SkyPath


class AzElRange(BaseAccessConstraint):
    @u.quantity_input(
        min_az=u.deg,
        max_az=u.deg,
        min_el=u.deg,
        max_el=u.deg,
        min_range=u.m,
        max_range=u.m,
    )
    def __init__(
        self,
        observer: SkyPath,
        target: SkyPath,
        min_az: u.deg = 0 * u.deg,
        max_az: u.deg = 360 * u.deg,
        min_el: u.deg = 0 * u.deg,
        max_el: u.deg = 90 * u.deg,
        min_range: u.m = 0 * u.m,
        max_range: u.m = np.inf * u.m,
        bounds_check: bool = False,
    ) -> None:
        """
        Determine access using ENU azimuth, elevation, and range values from an earth based position.

        Args:
            min_az (u.deg, optional): minimum allowable azimuth angle to target. Defaults to 0*u.deg.
            max_az (u.deg, optional): maximum allowable azimuth angle to target. Defaults to 360*u.deg.
            min_el (u.deg, optional): minimum allowable elevation angle to target. Defaults to 0*u.deg.
            max_el (u.deg, optional): maximum allowable elevation angle to target. Defaults to 90*u.deg.
            min_range (u.m, optional): minimum allowable range angle to target. Defaults to 0*u.m.
            max_range (u.m, optional): maximum allowable range to target. Defaults to np.inf*u.m.
            bounds_check (bool, optional): check the time bounds of the observer and target to ensure
            new state times are within the times that the objects were created with.
        """
        super().__init__()
        self.observer = observer
        self.target = target
        self.min_az = min_az
        self.max_az = max_az
        self.min_el = min_el
        self.max_el = max_el
        self.min_range = min_range
        self.max_range = max_range
        self.bounds_check = bounds_check

    def __call__(self, time: Time, *args, **kwargs) -> np.ndarray:
        """Constraints the times when the azimuth, elevation, and range values are satisfied from the observer to the target.

        Args:
            time (Time): time to check for azimuth, elevation and range constraints

        Returns:
            np.ndarray: boolean array representing times when this constraint succeeds vs. fails.
        """

        # NOTE: Using pymap3d is tremendously faster than having astropy calculate the AltAz frame from the observer to the target,
        # however there are cases for distant objects that this loses accuracy.
        az, el, rng = pymap3d.ecef2aer(
            *self.target.state_at(time, "itrs", bounds_check=self.bounds_check)
            .cartesian.xyz.to(u.m)
            .value,
            *pymap3d.ecef2geodetic(
                *self.observer.state_at(time, "itrs", bounds_check=self.bounds_check)
                .cartesian.xyz.to(u.m)
                .value
            )
        )
        return (
            (az >= self.min_az.value)
            & (az <= self.max_az.value)
            & (el >= self.min_el.value)
            & (el <= self.max_el.value)
            & (rng >= self.min_range.value)
            & (rng <= self.max_range.value)
        )
