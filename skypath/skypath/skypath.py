from typing import List, Union

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord
from astropy.time import Time

from skypath.access.access_funcs import AccessInterval, get_access
from skypath.access.constraints import BaseAccessConstraint
from skypath.coordinates import CoordinateInterpolator
from skypath.look_angles import BaseLookAngleStrategy
from skypath.skypath.creation import KinematicCreationMixin


class SkyPath(KinematicCreationMixin):

    def __init__(
        self,
        coordinate_frame: Union[
            SkyCoord, BaseCoordinateFrame, CoordinateInterpolator, "SkyPath"
        ],
    ) -> None:
        if isinstance(coordinate_frame, SkyPath):
            self.coordinate_frame = CoordinateInterpolator(
                coordinate_frame.coordinate_frame
            )
        else:
            self.coordinate_frame = CoordinateInterpolator(coordinate_frame)

    def access_to(
        self,
        target: "SkyPath",
        time: Time,
        constraints: List[BaseAccessConstraint] = [],
        find_precise_times: bool = True,
        precision: u.s = 0.1 * u.s,
    ) -> AccessInterval:
        return get_access(
            self.coordinate_frame,
            target.coordinate_frame,
            time,
            constraints,
            find_precise_times,
            precision,
        )

    def state_at(
        self, time: Time, frame: str = "itrs", bounds_check: bool = True
    ) -> "SkyPath":
        return SkyPath(self.coordinate_frame.state_at(time, frame, True, bounds_check))

    def look_angles_to(
        self, target: "SkyPath", time: Time, strategy: BaseLookAngleStrategy
    ) -> tuple:
        return strategy(self.coordinate_frame, target.coordinate_frame, time)
