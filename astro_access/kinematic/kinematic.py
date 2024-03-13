from astro_access.coordinates import CoordinateInterpolator
from astro_access.access.constraints import BaseAccessConstraint
from astro_access.access.access_funcs import AccessInterval, get_access
from astro_access.attitude import BaseAttitudeStrategy, Identity
from astro_access.kinematic.creation import KinematicCreationMixin
from astro_access.look_angles.look_angle_funcs import get_look_angles_to
from astro_access.look_angles import BaseLookAngleStrategy, NadirWithVelocityConstraint

from astropy.coordinates import SkyCoord, BaseCoordinateFrame
import astropy.units as u
from astropy.time import Time
from typing import Union, List
import numpy as np


class Kinematic(KinematicCreationMixin):
    
    def __init__(self, coordinate_frame: Union[SkyCoord, BaseCoordinateFrame, CoordinateInterpolator]) -> None:
        self.coordinate_frame = CoordinateInterpolator(coordinate_frame)
    
    def access_to(self,
                  target: 'Kinematic',
                  time: Time,
                  constraints: List[BaseAccessConstraint] = [],
                  find_precise_times: bool = True,
                  precision: u.s = 0.1 * u.s) -> AccessInterval:
        return get_access(self.coordinate_frame, target.coordinate_frame, time, constraints, find_precise_times, precision)
    
    def coordinates_at(self, time: Time, frame: str = 'itrs', bounds_check:bool = True) -> 'Kinematic':
        return Kinematic(self.coordinate_frame.state_at(time, frame, True, bounds_check))
    
    def look_angles_to(self,
                       target: 'Kinematic',
                       time: Time,
                       strategy: BaseLookAngleStrategy) -> tuple:        
        return strategy(self.coordinate_frame, target.coordinate_frame, time)
            