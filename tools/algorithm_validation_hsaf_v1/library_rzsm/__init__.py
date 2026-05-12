from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("rzsm")
except PackageNotFoundError:
    __version__ = "unknown"
