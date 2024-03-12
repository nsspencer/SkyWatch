from astro_access.look_angles.strategies._base_look_angle_stragtegy import BaseLookAngleStrategy
from astro_access.coordinates.coordinate_interpolator import CoordinateInterpolator
from astropy.time import Time
import pymap3d
import astropy.units as u
from astropy.coordinates import AltAz, ITRS


class ENUAER(BaseLookAngleStrategy):
    
    def __init__(self, high_fidelity: bool = False) -> None:
        """
        Defines the look angles in the earths ENU frame.
        """
        super().__init__()
        self.high_fidelity = high_fidelity
    
    def __call__(self, observer: CoordinateInterpolator, target: CoordinateInterpolator, time: Time) -> tuple:
        
        if self.high_fidelity:
            target_state = target.state_at(time, 'itrs')
            alt_az = ITRS(target_state.cartesian.without_differentials(), obstime=time).transform_to(AltAz(location=observer.state_at(time, 'itrs').earth_location, obstime=time))
            az = alt_az.az
            el = alt_az.alt
            rng =alt_az.distance
        
        else:
            target_pos = target.state_at(time, 'itrs').cartesian.xyz.to(u.m).value
            itrs_state = observer.state_at(time, 'itrs')
            lats = itrs_state.earth_location.lat.value
            lons = itrs_state.earth_location.lon.value
            heights = itrs_state.earth_location.height.to(u.m).value
            
            az, el, rng  = pymap3d.ecef2aer(*target_pos, lats, lons, heights)
            az = az*u.deg
            el = el*u.deg
            rng = rng*u.m
            
        return az, el, rng