import astropy.units as u
import pymap3d
from astropy.coordinates import ITRS, AltAz
from astropy.time import Time

from skywatch.attitude import BaseAttitudeStrategy
from skywatch.look_angles.aert import AzElRangeTime
from skywatch.look_angles.base_look_angle import BaseLookAngleStrategy
from skywatch.skypath.skypath import SkyPath


class LocalTangentENU(BaseLookAngleStrategy):
    def __init__(
        self,
        observer: SkyPath,
        attitude_strategy: BaseAttitudeStrategy = None,
        use_astropy: bool = False,
    ) -> None:
        """
        Using East, North, Up in reference to Earth. Up is zenith from the plane tangent to the earths
        surface at the observers location.

        Args:
            use_astropy (bool, optional): Do calculations using astropy's AltAz frame. NOTE: This is slower, but potentially
            more accurate. Defaults to False.
        """
        super().__init__()
        self.observer = observer
        self.use_astropy = use_astropy
        self.attitude_strategy = attitude_strategy

    def get_look_angles(
        self,
        target: SkyPath,
        time: Time,
    ) -> AzElRangeTime:
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

            enu = pymap3d.ecef2enu(*target_pos, lats, lons, heights)
            if self.attitude_strategy is not None:
                attitude = self.attitude_strategy.at(time)
                enu = attitude.inv().apply(enu)

            az, el, rng = pymap3d.enu2aer(*enu)
            az = az * u.deg
            el = el * u.deg
            rng = rng * u.m

        return AzElRangeTime(az, el, rng, time)
