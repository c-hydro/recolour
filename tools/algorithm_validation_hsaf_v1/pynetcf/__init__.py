from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pynetcf_local")
except PackageNotFoundError:
    __version__ = "unknown"
