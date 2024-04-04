import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, EarthLocation, SkyCoord, get_body
from astropy.time import Time


class SkyPathCreationMixin:

    @classmethod
    @u.quantity_input(latitude=u.deg, longitude=u.deg, altitude=u.m)
    def from_geodetic(
        cls,
        time: Time,
        latitude: u.Quantity,
        longitude: u.Quantity,
        altitude: u.Quantity,
    ):
        """
        Construct a SkyPath from a geodetic location.

        Args:
            time (Time): Time of the observation. For stationary coordinates, just use one time.
            latitude (u.Quantity): latitude in degrees.
            longitude (u.Quantity): longitude in degrees.
            altitude (u.Quantity): altitude in meters.

        Returns:
            SkyPath: SkyPath instance in the ITRS frame.
        """
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
        """
        Construct a SkyPath from ECEF X,Y,Z position and/or
        V_X,V_Y,V_Z velocity at a given time.

        NOTE: ECEF is NOT equivalent to the ITRS frame; however, it is
        often considered close enough and functionally comparable.

        Args:
            time (Time): Time of position and/or velocity measurements.
            x (u.Quantity): position X
            y (u.Quantity): position Y
            z (u.Quantity): posizion Z
            v_x (u.Quantity, optional): velocity X. Defaults to None.
            v_y (u.Quantity, optional): velocity Y. Defaults to None.
            v_z (u.Quantity, optional): velocity Z. Defaults to None.

        Returns:
            SkyPath: SkyPath instance in the ITRS frame.
        """
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
        """
        Construct a SkyPath from ECI X,Y,Z position and/or
        V_X,V_Y,V_Z velocity at a given time.

        Args:
            time (Time): Time of position and/or velocity measurements.
            x (u.Quantity): position X
            y (u.Quantity): position Y
            z (u.Quantity): posizion Z
            v_x (u.Quantity, optional): velocity X. Defaults to None.
            v_y (u.Quantity, optional): velocity Y. Defaults to None.
            v_z (u.Quantity, optional): velocity Z. Defaults to None.

        Returns:
            SkyPath: SkyPath instance in the GCRS frame.
        """
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
    def from_body(
        cls,
        time: Time,
        body: str,
        location: EarthLocation = None,
        ephemeris: str = None,
    ):
        """
        Construct a SkyPath for a solar system body as observed
        from a location on Earth in the `~astropy.coordinates.GCRS` reference
        system.

        Parameters
        ----------
        body : str or list of tuple
            The solar system body for which to calculate positions.  Can also be a
            kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
            kernel.
        time : `~astropy.time.Time`
            Time of observation.
        location : `~astropy.coordinates.EarthLocation`, optional
            Location of observer on the Earth.  If not given, will be taken from
            ``time`` (if not present, a geocentric observer will be assumed).
        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            ``astropy.coordinates.solar_system_ephemeris.set`` (which is
            set to 'builtin' by default).

        Returns:
            SkyPath: SkyPath instance representing the coordinates of a planetary body.
        """
        return cls(get_body(body, time, location, ephemeris))

    @classmethod
    def from_base_coordinate(cls, base: BaseCoordinateFrame):
        """
        Construct a SkyPath from a base coordinate frame like ITRS, GCRS, ICRS, etc.

        Args:
            base (BaseCoordinateFrame): base coordinate frame.

        Returns:
            SkyPath: SkyPath instance wrapping the base coordinate frame.
        """
        return cls(base)
