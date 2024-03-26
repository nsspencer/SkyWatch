import astropy.units as u
import portion as P
from astropy.time import Time


class TimeInterval(P.Interval):
    def __init__(self, *intervals):
        """
        Represents an interval of time.
        """
        super().__init__(*intervals)

    @property
    def total_duration(self) -> u.s:
        duration = 0
        for interval in self:
            duration += (
                Time(interval.upper) - Time(interval.lower)
            ).datetime.total_seconds()
        return duration * u.s
