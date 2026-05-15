from importlib.metadata import version, PackageNotFoundError

__all__ = ['metrics', 'scaling', 'temporal_matching',
           'timedate', 'time_series',
           'grid', 'io', 'colormaps']


try:
    __version__ = version("pytesmo_local")
except PackageNotFoundError:
    __version__ = "unknown"
