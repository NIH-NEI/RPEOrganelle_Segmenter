from pathlib import Path
from typing import Union

import numpy as np

PathLike = Union[str, Path]
SegmentationLike = Union[str, bool]
Stacklike = np.ndarray  # TODO: ndim = 3?, update to numpy typing
ArrayLike = Union[list, tuple, np.ndarray, set]
PathList = Union[list, PathLike]
