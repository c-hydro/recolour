# libraries
import warnings
import geopandas as gpd
import numpy as np
import pandas as pd

from repurpose.resample import resample_to_grid
import matplotlib.pylab as plt

import warnings
warnings.filterwarnings(
    "ignore",
    message="Possible more than 7 neighbours within 25000 m"
)

# method to return distance in metric using Haversine formula
def haversine_distance_matrix(lon1, lat1, lon2, lat2):
    """
    Compute Haversine distance matrix between two sets of lon/lat points (in degrees).
    Returns distance in meters.
    """
    R = 6371000  # Earth radius in meters
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    dlon = lon2[:, None] - lon1
    dlat = lat2[:, None] - lat1

    a = np.sin(dlat / 2)**2 + np.cos(lat2[:, None]) * np.cos(lat1) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return R * c

# method to fill missing values in a grid using nearest neighbors
def fill_missing_with_nearest_mean_radius(grid_values, mask_lons, mask_lats, k=8, max_radius_m=25000):
    """
    Fill NaN values using mean of up to `k` nearest neighbors within `max_radius_m` (meters),
    using Haversine distance.
    """
    # Flatten the grid
    grid_flat = grid_values.ravel()
    lat_flat = mask_lats.ravel()
    lon_flat = mask_lons.ravel()

    valid_idx = ~np.isnan(grid_flat)
    invalid_idx = np.isnan(grid_flat)

    if not np.any(invalid_idx):
        return grid_values  # Nothing to fill

    # Coordinates
    valid_lats = lat_flat[valid_idx]
    valid_lons = lon_flat[valid_idx]
    valid_vals = grid_flat[valid_idx]

    invalid_lats = lat_flat[invalid_idx]
    invalid_lons = lon_flat[invalid_idx]

    # Compute distance matrix (invalid_points x valid_points)
    dists = haversine_distance_matrix(valid_lons, valid_lats, invalid_lons, invalid_lats)  # meters

    filled_values = grid_flat.copy()

    for i in range(dists.shape[0]):
        neighbor_distances = dists[i]
        within_radius = neighbor_distances <= max_radius_m

        if np.any(within_radius):
            nearest_idxs = np.argsort(neighbor_distances[within_radius])[:k]
            neighbor_values = valid_vals[within_radius][nearest_idxs]
            mean_value = np.mean(neighbor_values)
            filled_values[np.where(invalid_idx)[0][i]] = mean_value
        # else: leave NaN

    return filled_values.reshape(grid_values.shape)

# method to resample points to grid
def resample_points2grid(
        points_lon, points_lat, points_value,
        mask_lons, mask_lats, mask_values,
        radius_of_influence=25, neighbours_min=1, neighbours_n=7,
        fill_value=np.nan, max_value=100, min_value=0,
        mask=False, plot=False):

    # check types
    if fill_value is None:
        fill_value = np.nan
    # convert km to meters
    radius_of_influence = radius_of_influence * 1000

    # filter points based on limits
    if max_value is not None:
        points_value[points_value > max_value] = np.nan
    if min_value is not None:
        points_value[points_value <= min_value] = np.nan

    # create DataFrame from points
    points_df = pd.DataFrame({'value': points_value, 'lat': points_lat, 'lon': points_lon})
    # Drop rows with NaN values in 'value'
    points_df = points_df.dropna(subset=['value'])
    # return cleaned points
    points_value = points_df['value'].to_numpy()
    points_lat, points_lon = points_df['lat'].to_numpy(), points_df['lon'].to_numpy()

    # target grid geometry
    min_lon, max_lon = np.min(mask_lons), np.max(mask_lons)
    min_lat, max_lat = np.min(mask_lats), np.max(mask_lats)
    height, width = mask_lats.shape
    area_extent = (min_lon, min_lat, max_lon, max_lat)

    # method to resample data to grid
    grid_obj = resample_to_grid(
        {'data': points_value},
        points_lon, points_lat, mask_lons, mask_lats, search_rad=radius_of_influence,
        min_neighbours=neighbours_min, neighbours=neighbours_n, fill_values=np.nan)

    # grid data organization
    grid_resampled = grid_obj['data']
    if max_value is not None:
        grid_resampled[grid_resampled > max_value] = fill_value
    if min_value is not None:
        grid_resampled[grid_resampled < min_value] = fill_value

    if plot:
        plt.figure()
        plt.imshow(grid_resampled, extent=area_extent, origin='lower', cmap='viridis')
        plt.figure()
        plt.imshow(mask_values, extent=area_extent, origin='lower', cmap='viridis')

    # Apply binary mask
    if mask:
        if mask_values is not None:
            mask_values = mask_values.astype(bool)  # Convert 0/1 to boolean
            grid_resampled = np.where(mask_values, grid_resampled, fill_value)

    grid_filled = fill_missing_with_nearest_mean_radius(
        grid_resampled, mask_lons, mask_lats, k=8, max_radius_m=radius_of_influence)

    if plot:
        plt.figure()
        plt.imshow(grid_resampled, extent=area_extent, origin='lower', cmap='viridis')
        plt.show(block=True)
        plt.figure()
        plt.imshow(grid_filled, extent=area_extent, origin='lower', cmap='viridis')
        plt.show(block=True)

    return grid_resampled, grid_filled

# method to filter grid values based on percentage of land coverage
def filter_grid(grid_values, percentage_values, threshold=0.5):
    """
    Filters coarse-grid values based on percentage of fine-grid land coverage.

    Parameters:
        grid_values (np.ndarray): 2D array of coarse-grid values.
        percentage_values (np.ndarray): 2D array, same shape, land percentage per coarse cell (0-1).
        threshold (float): Minimum percentage to retain value (e.g., 0.5 for 50%).

    Returns:
        np.ndarray: Filtered grid with np.nan where threshold not met.
    """
    mask = percentage_values >= threshold
    return np.where(mask, grid_values, np.nan)
