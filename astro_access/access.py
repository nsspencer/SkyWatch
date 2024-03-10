import astropy.units as u
import numpy as np
import portion as P

from typing import List
from astropy.time import Time
from astro_access.frame_interpolator import FrameInterpolator
from astro_access.constraints._base_constraint import BaseAccessConstraint



class Access:    
    
    @staticmethod
    def _get_true_ranges(bool_array: np.ndarray) -> List[tuple]:
        # Find the indices where the array changes from False to True or True to False
        change_indices = np.where(np.diff(bool_array))[0] + 1

        # If the first value is True, we need to add a 0 at the beginning
        if bool_array[0]:
            change_indices = np.insert(change_indices, 0, 0)

        # If the last value is True, we need to add the length of the array at the end
        if bool_array[-1]:
            change_indices = np.append(change_indices, len(bool_array))

        # Reshape the array into a 2D array where each row is a range of True values
        ranges = change_indices.reshape(-1, 2)

        # Convert the ranges to a list of tuples
        ranges = [tuple(r) for r in ranges]

        return ranges
    
    @staticmethod
    def _access_times(constrained_times: np.ndarray, time: Time) -> tuple:
        if not np.any(constrained_times):
            return [], [] # no access
        
        window_ranges = Access._get_true_ranges(constrained_times)
        valid_ranges = []
        access_times = []
        for range in window_ranges:
            new_time = time[range[0]:range[1]]
            if len(new_time) == 1:
                continue # need at least 2 times for access to be considered
            valid_ranges.append(range)
            access_times.append(new_time)
            
        return valid_ranges, access_times
    
    
    @u.quantity_input(precision = u.s) # seconds
    @staticmethod
    def get_access(observer: FrameInterpolator,
                   target: FrameInterpolator,
                   time: Time,
                   constraints: List[BaseAccessConstraint] = [],
                   find_precise_times: bool = True,
                   precision: u.s = 1 * u.s) -> P.Interval:
        
        # scale the precision to reflect the it in terms of seconds
        precision = (1 * u.s) / precision
        
        # first pass check of access
        constrined_times = [np.array([True] * len(time))]
        for constraint in constraints:
            constrined_times.append(constraint(observer, target, time, bounds_check=True))
        
        # get all the windows from the first pass
        valid_time = np.all(constrined_times, axis=0)
        valid_ranges, access_times = Access._access_times(valid_time, time)
        
        # if no access or if not using precise timing, set the first pass as the final access
        if not find_precise_times or len(constraints) == 0:
           final_access_times = access_times
           
        else:
            # calculate access at the precise time scale around the start and stop times of each window
            final_access_times = []
            for window_range, access_time in zip(valid_ranges, access_times):
                start_index = window_range[0]
                before_start_index = max(window_range[0] - 1, 0)
                
                if start_index == before_start_index:
                    exact_start_time = time[start_index] # nothing before the start time to interpolate to
                else:
                    # calculate number of steps between times to get desired precision
                    t1 = time[start_index]
                    t0 = time[before_start_index]
                    num_steps = max(int(((t1 - t0).datetime.total_seconds()) * precision.value), 2)                    
                    new_start_times = np.linspace(t0, t1, num_steps)
                    
                    # find the exact start time for this window
                    constrained_times = [np.array([True] * len(new_start_times))]
                    for constraint in constraints:
                        constrained_times.append(constraint(observer, target, new_start_times, bounds_check=False))
                    exact_start_time = new_start_times[np.all(constrained_times, axis=0)][0]

                # calculate the exact end time
                end_index = window_range[1]-1
                after_end_index = min(window_range[1], len(time)-1)
                
                if end_index == after_end_index:
                    exact_end_time = time[end_index] # nothing after the end time to interpolate to
                else:
                    # calculate number of steps between times to get desired precision
                    t0 = time[end_index]
                    t1 = time[after_end_index]
                    num_steps = max(int(((t1 - t0).datetime.total_seconds()) * precision.value), 2)
                    new_end_times = np.linspace(t0, t1, num_steps)
                    
                    # find the exact end time for this window
                    constrained_times = [np.array([True] * len(new_end_times))]
                    for constraint in constraints:
                        constrained_times.append(constraint(observer, target, new_end_times, bounds_check=False))
                    exact_end_time = new_end_times[np.all(constrained_times, axis=0)][-1]
                
                final_access_times.append((exact_start_time, exact_end_time))
       
       # compute the access windows using portion's intervals
        access = P.Interval()
        for access_time in final_access_times:
            t_start, t_end = access_time[0], access_time[-1]
            access = access.union(P.closed(t_start, t_end))
        
        return access
