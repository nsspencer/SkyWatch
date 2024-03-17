from dataclasses import dataclass

import astropy.units as u
from astropy.time import Time


@dataclass
class AzElRangeTime:
    azimuth: u.deg
    elevation: u.deg
    range: u.m
    time: Time
