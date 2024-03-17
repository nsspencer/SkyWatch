from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List

import astropy.units as u
import numpy as np
import tqdm
from astropy.time import Time

from skywatch.access import Access, AccessInterval
from skywatch.access.constraints import AzElRange
from skywatch.coordinates import SkyPath
from skywatch.utils.funcs import fibonacci_latitude_longitude


@dataclass
class CoveragePoint:
    latitude: u.deg
    longitude: u.deg
    interval: AccessInterval
    time_bound_lower: Time
    time_bound_upper: Time

    def _default_revisit_stats(self) -> u.s:
        if (
            self.interval.lower == self.time_bound_lower
            and self.interval.upper == self.time_bound_upper
        ):
            return 0 * u.s
        return np.inf * u.s

    @property
    def average_revisit_time(self) -> u.s:
        """
        Average revisit time is the time between the end of one access interval and the start of the next, all summed together.

        Returns:
            u.s: Returns np.inf * u.s if this coordinate has no access, returns 0 seconds if access
            was never lost, else returns seconds representing the average of all revisit times for this lat/lon coordinate.
        """
        if len(self.interval) < 2:
            return self._default_revisit_stats()

        lost_access_times = 0 * u.s
        for index, window in enumerate(self.interval[:-1]):
            lost_access_times += (
                self.interval[index + 1].lower - window.upper
            ).datetime.total_seconds() * u.s

        return lost_access_times / (len(self.interval) - 1)

    @property
    def max_revisit_time(self) -> u.s:
        """
        Maximum revisit time is the time between the end of one access interval and the start of the next.

        Returns:
            u.s: Returns None if this coordinate has no access, returns 0 seconds if access
            was never lost, else returns seconds representing the average of all revisit times for this lat/lon coordinate.
        """
        if len(self.interval) < 2:
            return self._default_revisit_stats()

        max_revisit_time = None
        for index, window in enumerate(self.interval[:-1]):
            loss_seconds = (
                self.interval[index + 1].lower - window.upper
            ).datetime.total_seconds() * u.s
            if max_revisit_time == None:
                max_revisit_time = loss_seconds
            elif max_revisit_time < loss_seconds:
                max_revisit_time = loss_seconds

        return max_revisit_time

    @property
    def min_revisit_time(self) -> u.s:
        """
        Maximum revisit time is the time between the end of one access interval and the start of the next.

        Returns:
            u.s: Returns None if this coordinate has no access, returns 0 seconds if access
            was never lost, else returns seconds representing the average of all revisit times for this lat/lon coordinate.
        """
        if len(self.interval) < 2:
            return self._default_revisit_stats()

        min_revisit_time = None
        for index, window in enumerate(self.interval[:-1]):
            loss_seconds = (
                self.interval[index + 1].lower - window.upper
            ).datetime.total_seconds() * u.s
            if min_revisit_time == None:
                min_revisit_time = loss_seconds
            elif min_revisit_time > loss_seconds:
                min_revisit_time = loss_seconds

        return min_revisit_time


class CoverageFilter(ABC):
    @abstractmethod
    def __call__(self, coverage_point: CoveragePoint) -> bool:
        """Return True when filter should remove the CoveragePoint, false when it should remain.

        Args:
            coverage_point (CoveragePoint): Point to check if it satisfies the filter.

        Returns:
            bool: True or False representing filter result.
        """
        pass


class GeoFilter(CoverageFilter):
    def __init__(
        self,
        min_latitude: u.deg = -90 * u.deg,
        max_latitude: u.deg = 90 * u.deg,
        min_longitude: u.deg = -180 * u.deg,
        max_longitude: u.deg = 180 * u.deg,
    ) -> None:
        super().__init__()
        self.min_latitude = min_latitude
        self.max_latitude = max_latitude
        self.min_longitude = min_longitude
        self.max_longitude = max_longitude

    def __call__(self, coverage_point: CoveragePoint) -> bool:
        if coverage_point.latitude < self.min_latitude:
            return True
        if coverage_point.latitude > self.max_latitude:
            return True
        if coverage_point.longitude < self.min_longitude:
            return True
        if coverage_point.longitude > self.max_longitude:
            return True
        return False


class CoverageResult:
    def __init__(self, *coverage_points) -> None:
        # statistics parameters
        self.percent_coverage = 0.0
        self.max_duration = 0.0 * u.s
        self.min_duration = None
        self.average_duration = 0.0 * u.s
        self.total_coverage_duration = 0.0 * u.s
        self.num_points_with_coverage = 0
        self.num_points_total = 0

        self.coverage_points = []
        for result in coverage_points:
            self.add_result(result)

    def add_result(self, coverage: CoveragePoint):
        if not isinstance(coverage, CoveragePoint):
            raise TypeError("Coverage must be of type CoveragePoint")

        self.coverage_points.append(coverage)
        self.num_points_total += 1

        # calculate statistics
        coverage_duration = coverage.interval.total_duration
        self.total_coverage_duration += coverage_duration

        if coverage_duration > self.max_duration:
            self.max_duration = coverage_duration

        if self.min_duration == None:
            self.min_duration = coverage_duration
        elif self.min_duration > coverage_duration:
            self.min_duration = coverage_duration

        if coverage_duration > 0.0 * u.s:
            self.num_points_with_coverage += 1
        self.percent_coverage = self.num_points_with_coverage / self.num_points_total

        self.average_duration = self.total_coverage_duration / self.num_points_total

    @property
    def max_revisit_time(self) -> u.s:
        _inf = np.inf * u.s

        max_revisit = None
        for result in self.coverage_points:
            result: CoveragePoint
            result_max = result.max_revisit_time
            if result_max == _inf:
                continue

            if max_revisit == None:
                max_revisit = result_max
                continue

            if result_max > max_revisit:
                max_revisit = result_max

        return max_revisit

    @property
    def min_revisit_time(self) -> u.s:
        _inf = np.inf * u.s

        min_revisit = None
        for result in self.coverage_points:
            result: CoveragePoint
            result_min = result.min_revisit_time
            if result_min == _inf:
                continue

            if min_revisit == None:
                min_revisit = result_min
                continue

            if result_min < min_revisit:
                min_revisit = result_min

        return min_revisit

    def filter(self, *filters: CoverageFilter):
        # check filter types
        for filter in filters:
            if not isinstance(filter, CoverageFilter):
                raise TypeError("Filters must be of type CoverageFilter")

        # run the filters over the results
        new_coverage = CoverageResult()
        for result in self.coverage_points:
            keep = True
            for filter in filters:
                if filter(result):
                    keep = False
                    break
            if keep:
                new_coverage.add_result(result)

        # return remaining coverage results
        return new_coverage


def calculate_coverage(
    satellites: List[SkyPath],
    time: Time,
    num_earth_points: int = 1000,
    min_elevation: u.deg = 0 * u.deg,
    precision: u.s = 0.1 * u.s,
    min_duration: u.s = 0.0 * u.s,
) -> CoverageResult:

    earth_point_coverages = []
    earth_points = []
    for lat, lon in fibonacci_latitude_longitude(num_earth_points):
        lat, lon = lat * u.deg, lon * u.deg
        earth_points.append(SkyPath.from_geodetic(time[0], lat, lon, 0 * u.m))
        earth_point_coverages.append(
            CoveragePoint(lat, lon, AccessInterval(), time[0], time[-1])
        )

    for point, coverage in tqdm.tqdm(
        zip(earth_points, earth_point_coverages),
        desc="Calculating Coverage",
        total=len(earth_points),
    ):
        for satellite in satellites:
            _access = (
                Access(
                    AzElRange(point, satellite, min_el=min_elevation),
                )
                .use_precise_endpoints(True)
                .set_precision(precision)
                .set_min_duration(min_duration)
                .calculate_at(time)
            )
            coverage.interval = coverage.interval.union(_access)

    return CoverageResult(*earth_point_coverages)
