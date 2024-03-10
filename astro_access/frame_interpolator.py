import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import sky_coordinate_parsers
import numpy as np
from scipy.interpolate import CubicSpline, CubicHermiteSpline


class SavedFrame:
    """
    Represents the a coordinate frames interpolation splines for quick position and velocity calculations.
    """
    def __init__(self, name: str, position_spline: CubicHermiteSpline, velocity_spline: CubicHermiteSpline) -> None:
        self.name = name
        self.position_spline = position_spline
        self.velocity_spline = velocity_spline
    
    def __repr__(self) -> str:
        return self.name
    
    
class FrameInterpolator(SkyCoord):
    
    def _copy_from(self, other: 'FrameInterpolator'):
        self._interpolation_allowed = other._interpolation_allowed
        self._min_original_time = other._min_original_time
        self._max_original_time = other._max_original_time
        self._saved_frames = other._saved_frames
        self._original_times = other._original_times
    
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        
        # call copy constructor
        if len(args) != 0 and isinstance(args[0], FrameInterpolator):
            self._copy_from(args[0])
            
        # construct new instance
        else:
            # check that the obstime has more than 1 value in it for interpolation
            if self.data.size > 1:
                self._interpolation_allowed = True
            else:
                self._interpolation_allowed = False
                
            self._min_original_time = np.min(self.obstime)
            self._max_original_time = np.max(self.obstime)
            self._saved_frames = []
            self._original_times = self.obstime
        
        
    def state_at(self, time: Time, frame: str, copy:bool = True) -> 'FrameInterpolator':
        """
        Interpolates an Astropy.BaseCoordinateFrame object at the given time(s) and saves the
        interpolation spline for later use. Therefore, you only pay the computation penalty
        for frame transforms once, and every subsequent time is multiple orders of magnitude
        quicker.
        
        NOTE: This has the potential to use a decent amount of memory.

        Args:
            time (Time): time(s) to get the coordinate frame values for
            frame (str): BaseCoordinateFrame name to get the results in

        Returns:
            FrameInterpolator: FrameInterpolator representing an Astropy SkyCoord in the coordinate system you requested. 
        """
        if not self._interpolation_allowed:
            new_frame = sky_coordinate_parsers._get_frame_class(frame)
            try:
                new_frame.obstime = time
            except AttributeError:
                pass
            
            return FrameInterpolator(self.transform_to(new_frame, merge_attributes=False))
        
        if np.min(time) < self._min_original_time or np.max(time) > self._max_original_time:
            raise ValueError(f"Cannot interpolate times that are outside the bounds of the original coordinate frame.\nTime bounds are: [{self._min_original_time}, {self._max_original_time}]")
            
        # using a list for speed, check if the conversion has already been calculated
        saved_frame = None
        for f in self._saved_frames:
            if f.name == frame:
                saved_frame = f
                break
        
        # if the frame exists, use it to interpolate the position and velocity to the given time
        if saved_frame is not None:
            _position = saved_frame.position_spline(time.unix) * u.m
            if saved_frame.velocity_spline is not None:
                _velocity = saved_frame.velocity_spline(time.unix) * (u.m / u.s)
            else:
                _velocity = [None, None, None]
            
            # return a new copy of this object in the requested frame
            new_coord = FrameInterpolator(x=_position[0], y=_position[1], z=_position[2],
                            v_x=_velocity[0], v_y=_velocity[1], v_z=_velocity[2],
                            frame=saved_frame.name, representation_type='cartesian', obstime=time, copy=False)
            
            if copy:
                new_coord._copy_from(self)
                
            return new_coord

        # we have not converted to this frame before, so we need to do the conversion and store it as a SavedFrame
        converted_frame = self.transform_to(frame)
        
        # get the position and velocity in cartesian XYZ coordinates of the converted frame
        position = converted_frame.cartesian.xyz.to(u.m)
        differentials = converted_frame.cartesian.differentials.get('s', None)
        if differentials is not None:
            velocity = differentials.d_xyz.to(u.m/u.s)
        else:
            velocity = None
            
        # now create the splines for the data that we have
        if velocity is not None:
            position_spline = CubicHermiteSpline(self._original_times.unix, position, velocity, axis=1, extrapolate=False)
            velocity_spline = position_spline.derivative()
            interpolated_position = position_spline(time.unix) * u.m
            interpolated_velocity = velocity_spline(time.unix) * (u.m / u.s)
        else:
            position_spline = CubicSpline(self._original_times.unix, position, axis=1, extrapolate=False)
            velocity_spline = None
            interpolated_position = position_spline(time.unix) * u.m
            interpolated_velocity = [None, None, None]
        
        # save the frame and its splines
        new_frame = SavedFrame(frame, position_spline, velocity_spline)
        self._saved_frames.append(new_frame)
        
        # construct the new frame from the interpolated coordinates as a copy of this object
        new_coord = FrameInterpolator(x=interpolated_position[0], y=interpolated_position[1], z=interpolated_position[2],
                        v_x=interpolated_velocity[0], v_y=interpolated_velocity[1], v_z=interpolated_velocity[2],
                        frame=frame, representation_type='cartesian', obstime=time, copy=False)
        
        if copy:
            new_coord._copy_from(self)
            
        return new_coord
