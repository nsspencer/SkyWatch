from scipy.spatial.transform.rotation import Rotation

from skywatch.attitude._base_attitude import BaseAttitudeStrategy
from skywatch.coordinates.skypath import SkyPath
from skywatch.look_angles.strategies.body_frame import BodyFrame
from skywatch.look_angles.strategies.local_tangent_enu import LocalTangentENU


class LookAnglesCreationMixin:
    @classmethod
    def from_local_enu(cls, observer: SkyPath, use_astropy: bool = False):
        strategy = LocalTangentENU(observer, use_astropy)
        return cls(strategy)

    @classmethod
    def from_body_frame(
        cls,
        observer: SkyPath,
        attitude_strategy: BaseAttitudeStrategy,
        attitude_offset: Rotation = None,
    ):
        strategy = BodyFrame(observer, attitude_strategy, attitude_offset)
        return cls(strategy)
