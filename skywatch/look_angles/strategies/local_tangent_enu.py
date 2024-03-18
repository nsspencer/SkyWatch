import astropy.units as u
import pymap3d
from astropy.coordinates import ITRS, AltAz
from astropy.time import Time

from skywatch.coordinates.skypath import SkyPath
from skywatch.look_angles.aer import AzElRangeTime
from skywatch.look_angles.strategies._base_strategy import BaseLookAngleStrategy


class LocalTangentENU(BaseLookAngleStrategy):
    def __init__(self, observer: SkyPath, use_astropy: bool = False) -> None:
        """
        Using East, North, Up in reference to Earth. Up is zenith from the plane tangent to the earths
        surface at the observers location.

        Args:
            observer (SkyPath): Observer to calculate the look angles in reference to.
            use_astropy (bool, optional): Do calculations using astropy's AltAz frame. NOTE: This is slower, but potentially
            more accurate. Defaults to False.
        """
        super().__init__()
        self.observer = observer
        self.use_astropy = use_astropy

    def __call__(self, time: Time, target: SkyPath, *args, **kwargs) -> AzElRangeTime:
        if self.use_astropy:
            target_state = target.state_at(time, "itrs")
            alt_az = ITRS(
                target_state.cartesian.without_differentials(), obstime=time
            ).transform_to(
                AltAz(
                    location=self.observer.state_at(time, "itrs").earth_location,
                    obstime=time,
                )
            )
            az = alt_az.az
            el = alt_az.alt
            rng = alt_az.distance

        else:
            target_pos = target.state_at(time, "itrs").cartesian.xyz.to(u.m).value
            itrs_state = self.observer.state_at(time, "itrs")
            lats = itrs_state.earth_location.lat.value
            lons = itrs_state.earth_location.lon.value
            heights = itrs_state.earth_location.height.to(u.m).value

            az, el, rng = pymap3d.ecef2aer(*target_pos, lats, lons, heights)
            az = az * u.deg
            el = el * u.deg
            rng = rng * u.m

        return AzElRangeTime(az, el, rng, time)