from abc import ABC, abstractmethod

import numpy as np
from astropy.time import Time


class BaseAccessConstraint(ABC):
    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def __call__(self, time: Time, *args, **kwargs) -> np.ndarray:
        pass
