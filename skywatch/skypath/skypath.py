import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from .creation import SkyPathCreationMixin


class SkyPath(SkyCoord, SkyPathCreationMixin):

    class SavedFrame:
        """
        Represents the a coordinate frames interpolation splines for quick position and velocity calculations.
        """

        def __init__(
            self,
            name: str,
            position_spline: CubicHermiteSpline,
            velocity_spline: CubicHermiteSpline,
        ) -> None:
            self.name = name
            self.position_spline = position_spline
            self.velocity_spline = velocity_spline

        def __repr__(self) -> str:
            return self.name

    def _copy_from(self, other: "SkyPath"):
        self._interpolation_allowed = other._interpolation_allowed
        self._min_original_time = other._min_original_time
        self._max_original_time = other._max_original_time
        self._saved_frames = other._saved_frames
        self._original_times = other._original_times
        self._original_frame = other._original_frame

    def __init__(self, *args, **kwargs) -> None:
        """
        SkyPath has equivalent functionality to the Astropy SkyCoord class, except
        it includes interpolation to get the coordinate state(s) at any time
        between the minimum and maximum obstime of the coordinates that this instance
        was created with.

        This class introduces a new *state_at* method that acts as the
        SkyCoord.transform_to(*frame*) method, but saves the state of this SkyCoord
        in the new frame for future calls to the state_at method. Future calls return
        the position of the SkyPath in the given frame using the saved and interpolated
        state. This allows position and velocity to be determined at any point in time
        for a SkyPath, even if the time was not present in the obstime array upon creation
        of the object.

        Therefore, this class is useful for simulations where you may want the coordinates
        of an object in multiple frames at multiple points in time, but dont want to pay the
        cost to transform the object into the new frame every time you need the coordinates.
        """
        super().__init__(*args, **kwargs)

        # call copy constructor
        if len(args) != 0 and isinstance(args[0], SkyPath):
            self._copy_from(args[0])

        # construct new instance
        else:
            # check that the obstime has more than 1 value in it for interpolation
            if self.data.size > 1:
                self._interpolation_allowed = True
            else:
                self._interpolation_allowed = False

            self._min_original_time = None
            self._max_original_time = None
            self._saved_frames = []
            self._original_times = self.obstime
            self._original_frame = self

    @property
    def min_time(self) -> Time:
        if self._min_original_time is not None:
            return self._min_original_time
        self._min_original_time = np.min(self._original_times)
        return self._min_original_time

    @property
    def max_time(self) -> Time:
        if self._max_original_time is not None:
            return self._max_original_time
        self._max_original_time = np.max(self._original_times)
        return self._max_original_time

    def state_at(self, time: Time, frame: str, copy: bool = True, bounds_check: bool = True) -> "SkyPath":
        """
        Interpolates this SkyPath at the given time(s) and saves the
        interpolation spline for later use. Therefore, you only pay the computation penalty
        for frame transforms once, and every subsequent call this method uses the saved frames
        interpolators to provide quick coordinate results for the provided time(s).

        Args:
            time (Time): time(s) to get the coordinate frame values for.
            frame (str): BaseCoordinateFrame name to get the results in.
            Example: "itrs", "gcrs", "icrs", "teme", etc...

        Returns:
            SkyPath: SkyPath representing an Astropy SkyCoord in the coordinate system you requested.
        """
        if not self._interpolation_allowed:
            return SkyPath(SkyCoord(self.frame.cartesian, obstime=time, frame=self.frame.name).transform_to(frame))

        if bounds_check:
            if np.min(time) < self.min_time or np.max(time) > self.max_time:
                raise ValueError(
                    f"Cannot interpolate times that are outside the bounds of the original coordinate frame.\nTime bounds are: [{self._min_original_time}, {self._max_original_time}]"
                )

        # using a list for speed, check if the conversion has already been calculated
        saved_frame = None
        for f in self._saved_frames:
            if f.name == frame:
                saved_frame = f
                break

        # if the frame exists, use it to interpolate the position and velocity to the given time
        unix_time = time.unix
        if saved_frame is not None:
            _position = saved_frame.position_spline(unix_time) * u.m
            if saved_frame.velocity_spline is not None:
                _velocity = saved_frame.velocity_spline(unix_time) * (u.m / u.s)
            else:
                _velocity = [None, None, None]

            # return a new copy of this object in the requested frame
            new_coord = SkyPath(
                x=_position[0],
                y=_position[1],
                z=_position[2],
                v_x=_velocity[0],
                v_y=_velocity[1],
                v_z=_velocity[2],
                frame=saved_frame.name,
                representation_type="cartesian",
                obstime=time,
                copy=False,
            )

            if copy:
                new_coord._copy_from(self)

            return new_coord

        # we have not converted to this frame before, so we need to do the conversion and store it as a SavedFrame
        try:
            converted_frame = self._original_frame.transform_to(frame)
        except AttributeError as err:
            # handle errors when coming from non terrestrial coordinate system
            if err.name == "to_geodetic":
                converted_frame = self.state_at(self._original_times, "itrs", copy=False).transform_to(frame)
            else:
                raise err

        # get the position and velocity in cartesian XYZ coordinates of the converted frame
        position = converted_frame.cartesian.xyz.to(u.m)
        differentials = converted_frame.cartesian.differentials.get("s", None)
        if differentials is not None:
            velocity = differentials.d_xyz.to(u.m / u.s)
        else:
            velocity = None

        # now create the splines for the data that we have
        if velocity is not None:
            position_spline = CubicHermiteSpline(
                self._original_times.unix, position, velocity, axis=1, extrapolate=False
            )
            velocity_spline = position_spline.derivative()
            interpolated_position = position_spline(unix_time) * u.m
            interpolated_velocity = velocity_spline(unix_time) * (u.m / u.s)
        else:
            position_spline = CubicSpline(self._original_times.unix, position, axis=1, extrapolate=False)
            velocity_spline = None
            interpolated_position = position_spline(unix_time) * u.m
            interpolated_velocity = [None, None, None]

        # save the frame and its splines
        new_frame = self.SavedFrame(frame, position_spline, velocity_spline)
        self._saved_frames.append(new_frame)

        # construct the new frame from the interpolated coordinates as a copy of this object
        new_coord = SkyPath(
            x=interpolated_position[0],
            y=interpolated_position[1],
            z=interpolated_position[2],
            v_x=interpolated_velocity[0],
            v_y=interpolated_velocity[1],
            v_z=interpolated_velocity[2],
            frame=frame,
            representation_type="cartesian",
            obstime=time,
            copy=False,
        )

        if copy:
            new_coord._copy_from(self)

        return new_coord
