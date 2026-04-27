"""
Library Features:

Name:           lib_utils_analysis_grid
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np

try:
    from pyresample.geometry import SwathDefinition
    from pyresample.kd_tree import resample_nearest
except ImportError:
    SwathDefinition = None
    resample_nearest = None

from config_utils import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to smooth map
def build_smooth_map(data_map, domain_mask, fill_value,
                       method="gaussian", sigma=1.0, size=3):
    """
    Smooth a map only over valid pixels, preserving no-data.

    Parameters
    ----------
    data_map : np.ndarray
        Input map (e.g. soil_moisture_map)
    domain_mask : np.ndarray (bool)
        Valid domain mask
    fill_value : float or np.nan
        No-data value
    method : str
        "gaussian" or "mean"
    sigma : float
        Sigma for gaussian
    size : int
        Window size for mean filter

    Returns
    -------
    smoothed_map : np.ndarray
    """
    from scipy.ndimage import gaussian_filter, uniform_filter

    data_map = np.asarray(data_map, dtype=np.float32)

    # valid data mask
    if np.isnan(fill_value):
        valid = np.isfinite(data_map) & domain_mask
    else:
        valid = (
            np.isfinite(data_map) &
            (~np.isclose(data_map, fill_value, atol=1e-6)) &
            domain_mask
        )

    # prepare values and weights
    values = np.zeros_like(data_map, dtype=np.float32)
    weights = np.zeros_like(data_map, dtype=np.float32)

    values[valid] = data_map[valid]
    weights[valid] = 1.0

    # apply filter
    if method == "gaussian":
        smooth_values = gaussian_filter(values, sigma=float(sigma))
        smooth_weights = gaussian_filter(weights, sigma=float(sigma))
    elif method == "mean":
        smooth_values = uniform_filter(values, size=int(size))
        smooth_weights = uniform_filter(weights, size=int(size))
    else:
        raise ValueError("method must be 'gaussian' or 'mean'")

    # normalize
    smoothed = np.full_like(data_map, np.float32(fill_value), dtype=np.float32)

    valid_smooth = smooth_weights > 0
    smoothed[valid_smooth] = (
        smooth_values[valid_smooth] / smooth_weights[valid_smooth]
    )

    # enforce domain
    smoothed[~domain_mask] = np.float32(fill_value)

    return smoothed
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to create mask based on pixel extention
def build_mask_by_pixel_extension(rows, cols, domain_mask, radius_pixels):
    mask = np.zeros(domain_mask.shape, dtype=bool)

    if radius_pixels < 0:
        raise ValueError("radius_pixels must be >= 0")

    for r, c in zip(rows, cols):
        r0 = max(0, r - radius_pixels)
        r1 = min(domain_mask.shape[0], r + radius_pixels + 1)
        c0 = max(0, c - radius_pixels)
        c1 = min(domain_mask.shape[1], c + radius_pixels + 1)

        block_domain = domain_mask[r0:r1, c0:c1]
        mask[r0:r1, c0:c1][block_domain] = True

    mask &= domain_mask
    return mask
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to interpolate points to grid
def interpolate_points_to_grid(
        src_lons, src_lats, src_vals,
        grid_lons, grid_lats, domain_mask,
        roi_km, fill_value):

    out = np.full(domain_mask.shape, np.float32(fill_value), dtype=np.float32)

    if roi_km <= 0:
        raise RuntimeError("roi_km must be > 0")

    if SwathDefinition is None or resample_nearest is None:
        raise ImportError(
            "pyresample is required for interpolate_points_to_grid. "
            "Install it with: pip install pyresample"
        )

    tgt_lons = grid_lons[domain_mask]
    tgt_lats = grid_lats[domain_mask]

    if src_lons.size == 0:
        out[~domain_mask] = np.float32(fill_value)
        return out

    src_def = SwathDefinition(lons=src_lons, lats=src_lats)
    tgt_def = SwathDefinition(lons=tgt_lons, lats=tgt_lats)

    mapped = resample_nearest(
        source_geo_def=src_def,
        data=src_vals,
        target_geo_def=tgt_def,
        radius_of_influence=roi_km * 1000.0,
        fill_value=fill_value,
        epsilon=0.0,
    )

    out[domain_mask] = np.asarray(mapped, dtype=np.float32)
    out[~domain_mask] = np.float32(fill_value)

    return out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to apply masks
def apply_mask_filter(data, filter_mask, domain_mask, fill_value):
    out = data.copy()
    out[~filter_mask] = np.float32(fill_value)
    out[~domain_mask] = np.float32(fill_value)
    return out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to build mask boundary
def build_mask_boundary(data_map, domain_mask, fill_value, radius=1):
    """
    Build a boundary map of valid pixels that touch at least one no-data pixel
    inside the domain mask.

    Parameters
    ----------
    data_map : np.ndarray
        Input 2D map (for example soil_moisture_map).
    domain_mask : np.ndarray of bool
        True where the domain is valid, False outside the domain.
    fill_value : float
        No-data value used in the map. Can be np.nan.
    radius : int, optional
        Neighbourhood radius in pixels.
        radius=1 -> 3x3 window
        radius=2 -> 5x5 window

    Returns
    -------
    boundary_mask : np.ndarray of bool
        True where a valid pixel touches at least one no-data pixel
        inside the domain.
    boundary_map : np.ndarray
        Map containing original values only on boundary pixels,
        fill_value elsewhere.
    valid_data : np.ndarray of bool
        Valid-data mask inside the domain.
    no_data : np.ndarray of bool
        No-data mask inside the domain.
    """
    from scipy.ndimage import maximum_filter

    if radius < 1:
        raise ValueError("radius must be >= 1")

    data_map = np.asarray(data_map)
    domain_mask = np.asarray(domain_mask, dtype=bool)

    if data_map.shape != domain_mask.shape:
        raise ValueError("data_map and domain_mask must have the same shape")

    if np.isnan(fill_value):
        valid_data = np.isfinite(data_map) & domain_mask
    else:
        valid_data = (
            np.isfinite(data_map) &
            (~np.isclose(data_map, fill_value, atol=1e-6)) &
            domain_mask
        )

    no_data = (~valid_data) & domain_mask

    window_size = 2 * radius + 1
    no_data_neigh = maximum_filter(no_data.astype(np.uint8), size=window_size) > 0
    no_data_neigh = no_data_neigh & domain_mask

    boundary_mask = valid_data & no_data_neigh

    if np.isnan(fill_value):
        boundary_map = np.full(data_map.shape, np.nan, dtype=np.float32)
    else:
        boundary_map = np.full(data_map.shape, np.float32(fill_value), dtype=np.float32)

    boundary_map[boundary_mask] = data_map[boundary_mask].astype(np.float32)

    return boundary_mask, boundary_map, valid_data, no_data
# ----------------------------------------------------------------------------------------------------------------------

