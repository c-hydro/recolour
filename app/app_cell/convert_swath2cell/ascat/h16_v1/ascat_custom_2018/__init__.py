from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ascat_custom_2018")
except PackageNotFoundError:
    __version__ = "unknown"

from ascat_custom_2018.timeseries import *
from ascat_custom_2018.h_saf import *
