from typing import List

import astropy.units as u
import numpy as np
import portion as P
from astropy.time import Time

from .constraints._base_constraint import BaseAccessConstraint
from .time_interval import TimeInterval


class Access:
    def __init__(self, *constraints) -> None:
        """
        Access is used to calculate times when all provided constraints
        are satisfied.

        This class follows the builer pattern.

        Raises:
            TypeError: constraints must inherit from BaseAccessConstraint class.

        Args:
            consraints: BaseAccessConstraint objects representing the constraints
            that must be satisfied for access to return a non-empty TimeInterval.
        """
        self.constraints = list()
        for constraint in constraints:
            if not isinstance(constraint, BaseAccessConstraint):
                raise TypeError(
                    "All constraints must inherit from BaseAccessConstraint class"
                )
            self.constraints.append(constraint)

        self._precision = 0.001 * u.s
        self._precise_endpoints = False
        self._min_duration = None
        self._max_duration = None

    def use_precise_endpoints(self, value: bool) -> "Access":
        self._precise_endpoints = value
        return self

    def set_precision(self, precision: u.s) -> "Access":
        self._precision = precision
        return self

    def set_min_duration(self, duration: u.s) -> "Access":
        self._min_duration = duration
        return self

    def set_max_duration(self, duration: u.s) -> "Access":
        self._max_duration = duration
        return self

    def add_constraint(self, constraint: BaseAccessConstraint) -> "Access":
        """
        Add a constraint to this Access instance.

        Args:
            constraint (BaseAccessConstraint): Constraint that must be satisfied.

        Raises:
            TypeError: Must pass a subclass of BaseAccessConstraint.

        Returns:
            self: follows the builder pattern.
        """
        if not isinstance(constraint, BaseAccessConstraint):
            raise TypeError("Constraint must subclass BaseAccessConstraint.")

        self.constraints.append(constraint)
        return self

    def calculate_at(self, time: Time, *args, **kwargs) -> TimeInterval:
        return self.__call__(time, *args, **kwargs)

    def __call__(self, time: Time, *args, **kwargs) -> TimeInterval:
        if not isinstance(time, Time):
            time = Time(time)

        # first pass check of access
        constrined_times = [np.array([True] * time.size)]
        for constraint in self.constraints:
            constrined_times.append(constraint(time, *args, **kwargs))

        # get all the windows from the first pass
        valid_time = np.all(constrined_times, axis=0)
        valid_ranges, access_times = Access._access_times(valid_time, time)

        # if no access or if not using precise timing, set the first pass as the final access
        if not self._precise_endpoints or len(self.constraints) == 0:
            final_access_times = access_times

        else:
            # scale the precision to reflect the it in terms of seconds
            precision = (1 * u.s) / self._precision

            # calculate access at the precise time scale around the start and stop times of each window
            final_access_times = []
            for window_range, access_time in zip(valid_ranges, access_times):
                start_index = window_range[0]
                before_start_index = max(window_range[0] - 1, 0)

                if start_index == before_start_index:
                    exact_start_time = time[
                        start_index
                    ]  # nothing before the start time to interpolate to
                else:
                    # calculate number of steps between times to get desired precision
                    t1 = time[start_index]
                    t0 = time[before_start_index]
                    num_steps = max(
                        int(((t1 - t0).datetime.total_seconds()) * precision.value),
                        2,
                    )
                    new_start_times = np.linspace(t0, t1, num_steps)

                    # find the exact start time for this window
                    constrained_times = [np.array([True] * len(new_start_times))]
                    for constraint in self.constraints:
                        constrained_times.append(
                            constraint(new_start_times, *args, **kwargs)
                        )
                    exact_start_time = new_start_times[
                        np.all(constrained_times, axis=0)
                    ][0]

                # calculate the exact end time
                end_index = window_range[1] - 1
                after_end_index = min(window_range[1], time.size - 1)

                if end_index == after_end_index:
                    exact_end_time = time[
                        end_index
                    ]  # nothing after the end time to interpolate to
                else:
                    # calculate number of steps between times to get desired precision
                    t0 = time[end_index]
                    t1 = time[after_end_index]
                    num_steps = max(
                        int(((t1 - t0).datetime.total_seconds()) * precision.value),
                        2,
                    )
                    new_end_times = np.linspace(t0, t1, num_steps)

                    # find the exact end time for this window
                    constrained_times = [np.array([True] * len(new_end_times))]
                    for constraint in self.constraints:
                        constrained_times.append(
                            constraint(new_end_times, *args, **kwargs)
                        )
                    exact_end_time = new_end_times[np.all(constrained_times, axis=0)][
                        -1
                    ]

                final_access_times.append((exact_start_time, exact_end_time))

        # compute the access windows using portion's intervals
        access = TimeInterval()
        for access_time in final_access_times:
            t_start, t_end = access_time[0], access_time[-1]
            if self._min_duration is not None:
                if (
                    t_end - t_start
                ).datetime.total_seconds() * u.s < self._min_duration:
                    continue
            if self._max_duration is not None:
                if (
                    t_end - t_start
                ).datetime.total_seconds() * u.s > self._max_duration:
                    continue
            access = access.union(P.closed(t_start, t_end))

        return access

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
            return [], []  # no access

        window_ranges = Access._get_true_ranges(constrained_times)
        valid_ranges = []
        access_times = []
        for range in window_ranges:
            new_time = time[range[0] : range[1]]
            if len(new_time) == 1:
                new_time = np.array([time[range[0]]] * 2)
            valid_ranges.append(range)
            access_times.append(new_time)

        return valid_ranges, access_times
