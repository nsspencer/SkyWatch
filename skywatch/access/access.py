from typing import List

import astropy.units as u
import numpy as np
import portion as P
from astropy.time import Time

from .base_constraint import BaseAccessConstraint
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
        self._only_check_failed_constraints = False

    def use_precise_endpoints(self, value: bool) -> "Access":
        """
        Toggle using high precision endpoints to get accurate
        start and stop times for access intervals.

        Args:
            value (bool): True or false

        Returns:
            Access: self
        """
        self._precise_endpoints = value
        return self

    def set_precision(self, precision: u.s = 0.001 * u.s) -> "Access":
        """
        Represents the delta time of the high precision constraint
        checks when *use_precise_endpoints* is True.

        Args:
            precision (u.s): delta time value in seconds.
            Default: 0.001 seconds.

        Returns:
            Access: self
        """
        self._precision = precision
        return self

    def set_min_duration(self, duration: u.s) -> "Access":
        """
        Sets the minimum duration of the resulting TimeInterval
        when running access against all constraints. Successful
        access intervals must be greater than this duration to be
        included in the TimeInterval output.

        Args:
            duration (u.s): minimum duration in seconds.

        Returns:
            Access: self
        """
        self._min_duration = duration
        return self

    def set_max_duration(self, duration: u.s) -> "Access":
        """
        Sets the maximum duration of the resulting TimeInterval
        when running access against all constraints. Successful
        access intervals must be less than this duration to be
        included in the TimeInterval output.

        Args:
            duration (u.s): maximum duration in seconds.

        Returns:
            Access: self
        """
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

    def only_check_failed_constraints(self, value: bool) -> "Access":
        """
        When using precise endpoints, setting this to True will use an optimization
        to only re-check the constraints that have failed before the start and stop
        time of the access interval(s). If a constraint did not fail at the time period
        directly before or after the access interval, then it will not be recalculated.

        NOTE: If you are using a sufficiently large enough gap in your time intervals
        to check access during, then this has the potential to jump over a precise time
        that would have failed for a given constraint. However, if you know that your
        access constraints do not have a precision of less than the delta time between
        your time values, then this will not be an issue.

        NOTE: The defalt is set to False.

        Args:
            value (bool): True or False for whether or not to recheck only the failed
            constraints.

        Returns:
            Access: self
        """
        self._only_check_failed_constraints = bool(value)
        return self

    def calculate_at(self, time: Time, *args, **kwargs) -> TimeInterval:
        """
        Calculates access using all access constraints in this Access instance
        at the specified time(s).

        Args:
            time (Time): Time(s) to check for calculate access intervals.

        Returns:
            TimeInterval: Resulting TimeIntervals that represent time(s) when all the
            access constraints are successful.
        """
        return self.__call__(time, *args, **kwargs)

    def __call__(self, time: Time, *args, **kwargs) -> TimeInterval:
        if not isinstance(time, Time):
            time = Time(time)

        # first pass check of access
        original_constrained_times = [np.array([True] * time.size)]
        for constraint in self.constraints:
            original_constrained_times.append(constraint(time, *args, **kwargs))

        # get all the windows from the first pass
        valid_time = np.all(original_constrained_times, axis=0)
        valid_ranges, access_times = Access._access_times(valid_time, time)

        # if no access or if not using precise timing, set the first pass as the final access
        if not self._precise_endpoints or len(self.constraints) == 0:
            final_access_times = access_times
            return self._compute_final_access_interval(final_access_times)

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
                exact_start_time = self.__get_exact_times(
                    time,
                    start_index,
                    before_start_index,
                    precision,
                    original_constrained_times,
                    *args,
                    **kwargs
                )[0]

            # calculate the exact end time
            end_index = window_range[1] - 1
            after_end_index = min(window_range[1], time.size - 1)

            if end_index == after_end_index:
                exact_end_time = time[
                    end_index
                ]  # nothing after the end time to interpolate to
            else:
                exact_end_time = self.__get_exact_times(
                    time,
                    end_index,
                    after_end_index,
                    precision,
                    original_constrained_times,
                    *args,
                    **kwargs
                )[-1]

            final_access_times.append((exact_start_time, exact_end_time))

        return self._compute_final_access_interval(final_access_times)

    def __get_exact_times(
        self,
        time: Time,
        index_0: int,
        index_1: int,
        precision: u.Quantity,
        original_constrained_times: list,
        *args,
        **kwargs
    ) -> Time:
        """
        Helper to calculate the exact times between one index and another
        from the original times.

        Args:
            time (Time): original time values.
            index_0 (int): index of value 0
            index_1 (int): index of value 1
            precision (u.Quantity): precision in seconds
            original_constrained_times (list): boolean array of original contraint pass/fails

        Returns:
            Time: array of Time values.
        """
        # calculate number of steps between times to get desired precision
        t0 = time[index_0]
        t1 = time[index_1]
        num_steps = max(
            int(((t1 - t0).datetime.total_seconds()) * precision.value),
            2,
        )
        new_times = np.linspace(t0, t1, num_steps)

        # find the exact end time for this window
        constrained_times = [np.array([True] * len(new_times))]
        for index, constraint in enumerate(self.constraints):
            # check if the constraint needs to be computed again like above
            if (
                self._only_check_failed_constraints
                and original_constrained_times[index + 1][index_1]
            ):
                continue

            # need to get the precise time when this constraint turned False
            constrained_times.append(constraint(new_times, *args, **kwargs))

        return new_times[np.all(constrained_times, axis=0)]

    def _compute_final_access_interval(self, final_access_times: list):
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
