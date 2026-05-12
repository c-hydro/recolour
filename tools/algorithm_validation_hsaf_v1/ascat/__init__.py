from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ascat_custom_2018")
except PackageNotFoundError:
    __version__ = "unknown"

from ascat.timeseries import *
from ascat.h_saf import *
