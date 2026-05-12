from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("gldas_local")
except PackageNotFoundError:
    __version__ = "unknown"
