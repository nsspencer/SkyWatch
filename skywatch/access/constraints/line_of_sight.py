from typing import Union

import astropy.units as u
import numpy as np
from astropy.coordinates import BaseCoordinateFrame
from astropy.time import Time

from skywatch.access.base_constraint import BaseAccessConstraint
from skywatch.skypath import SkyPath


class LineOfSight(BaseAccessConstraint):
    def __init__(
        self,
        observer: SkyPath,
        target: SkyPath,
        obstructor: Union[SkyPath, BaseCoordinateFrame],
        sma=6378.137 * u.km,
        smi=6356.7523 * u.km,
        use_frame: str = "itrs",
        when_obstructed: bool = False,
    ) -> None:
        """
        Line of sight constraint which checks that the line of sight between two positions
        is not blocked by some body with a given semi major and semi minor access.

        Args:
            obstructor (Union[SkyPath, BaseCoordinateFrame]): coordinate of the body which potentially blocks line of sight.
            sma (Quantity, optional): semi major axis of the body. Defaults to 6378.137*u.km.
            smi (Quantity, optional): semi minor axis of the body. Defaults to 6356.7523*u.km.
            use_frame (str, optional): coordinate frame name to calculate access in. Defaults to 'itrs'.
        """
        super().__init__()
        assert (
            observer is not target
        ), "Observer cannot be the same object as the target"
        self.observer = observer
        self.target = target
        self.when_obstructed = when_obstructed

        self.sma = sma
        self.smi = smi
        self.use_frame = use_frame
        if isinstance(obstructor, SkyPath):
            self.obstructor = obstructor
        elif isinstance(obstructor, BaseCoordinateFrame):
            self.obstructor = SkyPath(obstructor)
        else:
            raise TypeError(
                "Body must be a SkyPath, CoordinateInterpolator, or BaseCoordinateFrame."
            )

    def __call__(self, time: Time, *args, **kwargs) -> np.ndarray:
        """
        Returns the subset of times when the earth does not obstruct the two Kinematics.

        Args:
            observer (CoordinateInterpolator): Coordinate to get line of sight access from
            target (CoordinateInterpolator): Coordinate to get line of sight access to
            time (Time): Times to check for line of sight access
            bounds_check (bool, optional): Whether or not to check the provided times are within the bounds
            of the observer and target CoordinateInterpolator.

        Returns:
            Time: boolean array of times when this SkyPath has line of sight access to the target.
        """
        pos1 = self.observer.state_at(
            time, self.use_frame, bounds_check=True
        ).cartesian.xyz
        pos2 = self.target.state_at(
            time, self.use_frame, bounds_check=True
        ).cartesian.xyz
        obstructor_pos = self.obstructor.state_at(
            time, self.use_frame, bounds_check=True
        ).cartesian.xyz

        los_times = self._line_of_sight_body(
            pos1.T, pos2.T, obstructor_pos.T, self.sma, self.smi
        )
        if self.when_obstructed:
            return los_times
        return ~los_times

    @staticmethod
    def _line_of_sight_body(
        pos1: np.ndarray,
        pos2: np.ndarray,
        pos_body: np.ndarray,
        sma: u.km,
        smi: u.km,
    ) -> np.ndarray:
        """
        Determines whether a planetary body intersects the line of sight between two vectors.

        Parameters:
        pos1 (np.ndarray): The position vector of the observer in 3D space.
        pos2 (np.ndarray): The position vector of the target in 3D space.
        pos_body (np.ndarray): The position vector of the body (e.g., Earth) in 3D space.
        sma (u.km): The semi-major axis of the planetary body, default is 6378.137 km (Earth).
        smi (u.km): The semi-minor axis of the planetary body, default is 6356.752 km (Earth).

        Returns:
        np.ndarray: A boolean array indicating True when line of sight is not impeded by the body, and False when it is for each pair of position vectors.
        """
        # Make sure the arrays have the same shape
        assert pos1.shape == pos2.shape == pos_body.shape
        if pos1.ndim == 1:
            pos1 = np.reshape(pos1, (-1, 3))
            pos2 = np.reshape(pos2, (-1, 3))
            pos_body = np.reshape(pos_body, (-1, 3))

        # Calculate the line of sight vectors
        los_vectors = pos2 - pos1

        # Calculate the unit vectors of the line of sight
        los_unit_vectors = (
            los_vectors / np.linalg.norm(los_vectors, axis=1)[:, np.newaxis]
        )

        # Calculate the vectors from the observer to the body's center
        observer_to_body = pos_body - pos1

        # Calculate the scalar projections of observer_to_body onto the line of sight
        scalar_projections = np.einsum("ij,ij->i", observer_to_body, los_unit_vectors)

        # Calculate the closest points to the body's center on the line of sight
        closest_points = pos1 + scalar_projections[:, np.newaxis] * los_unit_vectors

        # Calculate the geocentric latitude of the closest points
        latitudes = np.arctan2(
            closest_points[:, 2],
            np.sqrt(closest_points[:, 0] ** 2 + closest_points[:, 1] ** 2),
        )

        # Calculate the body's radius at the geocentric latitude
        body_radius = np.sqrt(
            ((sma**2 * np.cos(latitudes)) ** 2 + (smi**2 * np.sin(latitudes)) ** 2)
            / ((sma * np.cos(latitudes)) ** 2 + (smi * np.sin(latitudes)) ** 2)
        )

        # Calculate the distances from the closest points to the body's center
        distances_to_body_center = np.linalg.norm(closest_points - pos_body, axis=1)

        # Check if the body intersects the line of sight
        intersections = np.logical_and(
            distances_to_body_center <= body_radius,
            np.logical_and(
                0 <= scalar_projections,
                scalar_projections <= np.linalg.norm(los_vectors, axis=1),
            ),
        )

        # Return the times that the body intersects the line of sight vector
        return intersections
