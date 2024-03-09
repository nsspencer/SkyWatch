from astro_access.constraints._base_constraint import BaseAccessConstraint
from astro_access.frame_interpolator import FrameInterpolator
import astropy.units as u
from astropy.time import Time
import numpy as np


class LineOfSight(BaseAccessConstraint):
    def __init__(self, body_position: FrameInterpolator, sma = 6378.137 * u.km, smi = 6356.7523 * u.km) -> None:
        """
        Line of sight constraint which checks that the line of sight between two positions
        is not blocked by some body with a given semi major and semi minor access.

        Args:
            body_position (FrameInterpolator): coordinate of the body which potentially blocks line of sight.
            sma (_type_, optional): semi major axis of the body. Defaults to 6378.137*u.km.
            smi (_type_, optional): semi minor axis of the body. Defaults to 6356.7523*u.km.
        """
        super().__init__()
        self.sma = sma
        self.smi = smi
        self.body_position = body_position
        
    def __call__(self, observer: FrameInterpolator, target: FrameInterpolator, time: Time) -> np.ndarray:
        """
        Returns the subset of times when the earth does not obstruct the two platforms.

        Args:
            observer (FrameInterpolator): Coordinate to get line of sight access from
            target (FrameInterpolator): Coordinate to get line of sight access to
            time (Time): Times to check for line of sight access

        Returns:
            Time: boolean array of times when this platform has line of sight access to the target.
        """
        assert self != target, "Cannot get line of sight to self"
        pos1 = observer.position.get_position(time) 
        pos2 = target.position.get_position(time)
        body_pos = self.body_position.get_position(time)
        return self._line_of_sight_body(pos1.T, pos2.T, body_pos.T, self.sma, self.smi)
        
    
    @staticmethod
    def _line_of_sight_body(pos1: np.ndarray, pos2: np.ndarray, pos_body: np.ndarray, sma = 6378.137 * u.km, smi = 6356.7523 * u.km) -> np.ndarray:
        """
        Determines whether a planetary body intersects the line of sight between two vectors.

        Parameters:
        pos1 (np.ndarray): The position vector of the observer in 3D space in ECEF coordinates.
        pos2 (np.ndarray): The position vector of the target in 3D space in ECEF coordinates.
        pos_body (np.ndarray): The position vector of the body (e.g., Earth) in 3D space in ECEF coordinates.
        sma (u.km): The semi-major axis of the planetary body, default is 6378.137 km (Earth).
        smi (u.km): The semi-minor axis of the planetary body, default is 6356.752 km (Earth).

        Returns:
        np.ndarray: A boolean array indicating True when line of sight is not impeded by the body, and False when it is for each pair of position vectors.
        """
        # Make sure the arrays have the same shape
        assert pos1.shape == pos2.shape == pos_body.shape

        # Calculate the line of sight vectors
        los_vectors = pos2 - pos1

        # Calculate the unit vectors of the line of sight
        los_unit_vectors = los_vectors / np.linalg.norm(los_vectors, axis=1)[:, np.newaxis]

        # Calculate the vectors from the observer to the body's center
        observer_to_body = pos_body - pos1

        # Calculate the scalar projections of observer_to_body onto the line of sight
        scalar_projections = np.einsum('ij,ij->i', observer_to_body, los_unit_vectors)

        # Calculate the closest points to the body's center on the line of sight
        closest_points = pos1 + scalar_projections[:, np.newaxis] * los_unit_vectors

        # Calculate the geocentric latitude of the closest points
        latitudes = np.arctan2(closest_points[:, 2], np.sqrt(closest_points[:, 0]**2 + closest_points[:, 1]**2))

        # Calculate the body's radius at the geocentric latitude
        body_radius = np.sqrt(((sma**2 * np.cos(latitudes))**2 + (smi**2 * np.sin(latitudes))**2) / ((sma * np.cos(latitudes))**2 + (smi * np.sin(latitudes))**2))

        # Calculate the distances from the closest points to the body's center
        distances_to_body_center = np.linalg.norm(closest_points - pos_body, axis=1)

        # Check if the body intersects the line of sight
        intersections = np.logical_and(distances_to_body_center <= body_radius, np.logical_and(0 <= scalar_projections, scalar_projections <= np.linalg.norm(los_vectors, axis=1)))

        # Return the times that the body does not intersect the line of sight vector
        return ~intersections