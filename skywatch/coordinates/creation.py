import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, EarthLocation, SkyCoord, get_body
from astropy.time import Time


class SkyPathCreationMixin:
    def __init__(self, *_, **__):  # HACK stub to make mypy happy
        ...

    @classmethod
    @u.quantity_input(latitude=u.deg, longitude=u.deg, altitude=u.m)
    def from_geodetic(
        cls,
        time: Time,
        latitude: u.Quantity,
        longitude: u.Quantity,
        altitude: u.Quantity,
    ):
        return cls(
            SkyCoord(
                EarthLocation(lat=latitude, lon=longitude, height=altitude).get_itrs(
                    obstime=time
                )
            )
        )

    @classmethod
    @u.quantity_input(x=u.m, y=u.m, z=u.m, v_x=u.m / u.s, v_y=u.m / u.s, v_z=u.m / u.s)
    def from_ECEF(
        cls,
        time: Time,
        x: u.Quantity,
        y: u.Quantity,
        z: u.Quantity,
        v_x: u.Quantity = None,
        v_y: u.Quantity = None,
        v_z: u.Quantity = None,
    ):
        return cls(
            SkyCoord(
                x=x,
                y=y,
                z=z,
                v_x=v_x,
                v_y=v_y,
                v_z=v_z,
                frame="itrs",
                representation_type="cartesian",
                obstime=time,
            )
        )

    @classmethod
    @u.quantity_input(x=u.m, y=u.m, z=u.m, v_x=u.m / u.s, v_y=u.m / u.s, v_z=u.m / u.s)
    def from_ECI(
        cls,
        time: Time,
        x: u.Quantity,
        y: u.Quantity,
        z: u.Quantity,
        v_x: u.Quantity = None,
        v_y: u.Quantity = None,
        v_z: u.Quantity = None,
    ):
        return cls(
            SkyCoord(
                x=x,
                y=y,
                z=z,
                v_x=v_x,
                v_y=v_y,
                v_z=v_z,
                frame="gcrs",
                representation_type="cartesian",
                obstime=time,
            )
        )

    @classmethod
    def from_body(cls, time: Time, body: str = "earth", *args, **kwargs):
        return cls(get_body(body, time, *args, **kwargs))

    @classmethod
    def from_base_coordinate(cls, base: BaseCoordinateFrame):
        return cls(base)
