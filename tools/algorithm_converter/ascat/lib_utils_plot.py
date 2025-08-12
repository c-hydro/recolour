# libraries
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import numpy as np

import matplotlib
matplotlib.use('TkAgg')

# method to plot a 2D grid with coastlines using Cartopy
def plot_grid(grid_lon, grid_lat, grid_values, title="Grid Plot with Coastlines"):

    grid_values[grid_values < 0] = np.nan
    grid_values[grid_values > 100] = np.nan# Set negative values to NaN for better visualization

    """
    Plots a 2D grid on a map using EPSG:4326 with coastlines.

    Parameters:
    - grid_lon: 2D array of longitudes
    - grid_lat: 2D array of latitudes
    - grid_values: 2D array of values to plot
    - title: Optional title for the plot
    """
    fig, ax = plt.subplots(figsize=(10, 6),
                           subplot_kw={'projection': ccrs.PlateCarree()})  # EPSG:4326

    # Plot the grid using pcolormesh
    mesh = ax.pcolormesh(grid_lon, grid_lat, grid_values, cmap='viridis', vmin=0, vmax=100,
                         shading='auto', transform=ccrs.PlateCarree())

    # Add coastlines and other features
    # Add map features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # Colorbar
    plt.colorbar(mesh, ax=ax, orientation='vertical', label='Values')

    # Set extent (optional, based on data bounds)
    ax.set_extent([
        np.min(grid_lon), np.max(grid_lon),
        np.min(grid_lat), np.max(grid_lat)
    ], crs=ccrs.PlateCarree())

    # Title and show
    plt.title(title)
    plt.tight_layout()
    plt.show(block=True)

    print()

# method to plot geospatial points with values on a map with coastlines
def plot_points(lons, lats, values, title="Geospatial Plot", cmap="viridis", size=100):
    """
    Plot geospatial points (EPSG:4326) with associated values on a map with coastlines.

    Args:
        lons (list or array): Longitudes (EPSG:4326)
        lats (list or array): Latitudes (EPSG:4326)
        values (list or array): Values corresponding to each point
        title (str): Plot title
        cmap (str): Colormap for values
        size (int): Marker size
    """
    if not (len(lons) == len(lats) == len(values)):
        raise ValueError("Input arrays lons, lats, and values must be the same length.")

    values[values < 0] = np.nan

    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Define bounds with margin
    margin = 2.0
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)
    extent = [lon_min - margin, lon_max + margin, lat_min - margin, lat_max + margin]

    # Set visible map extent and axis limits
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_xlim(extent[0], extent[1])  # Longitudes
    ax.set_ylim(extent[2], extent[3])  # Latitudes

    # Add map features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Plot data
    scatter = ax.scatter(
        lons, lats, c=values, cmap=cmap, s=size, vmin=0, vmax=100,
        edgecolor='k', transform=ccrs.PlateCarree()
    )

    plt.colorbar(scatter, ax=ax, label="Value")
    plt.title(title)
    plt.show(block=True)

# method to plot a 2D grid and geospatial points on the same map with coastlines
def plot_grid_and_points(grid_lon, grid_lat, grid_values,
                         point_lons, point_lats, point_values,
                         cmap='viridis',
                         vmin=0, vmax=100,
                         point_size=100, margin=2.0,
                         nan_point_color='white',
                         nan_point_edge='black',
                         nan_point_size=100,
                         title="Grid and Points Plot with Coastlines"):
    """
    Plots a 2D grid and geospatial points on the same map with a shared colormap and colorbar.
    NaN-valued points are shown explicitly with custom styling.
    """

    # --- Clean grid values ---
    grid_values = np.array(grid_values)
    grid_values[(grid_values < vmin) | (grid_values > vmax)] = np.nan

    # --- Convert point data to arrays ---
    point_lons = np.array(point_lons)
    point_lats = np.array(point_lats)
    point_values = np.array(point_values)

    # --- Split valid and NaN points ---
    valid_mask = (
        (point_values >= vmin) & (point_values <= vmax) &
        ~np.isnan(point_values) & ~np.isnan(point_lons) & ~np.isnan(point_lats)
    )
    nan_mask = ~valid_mask & ~np.isnan(point_lons) & ~np.isnan(point_lats)

    # --- Apply mask ---
    point_lons_valid = point_lons[valid_mask]
    point_lats_valid = point_lats[valid_mask]
    point_values_valid = point_values[valid_mask]

    point_lons_nan = point_lons[nan_mask]
    point_lats_nan = point_lats[nan_mask]

    # --- Shared color normalization ---
    norm = Normalize(vmin=vmin, vmax=vmax)

    # --- Set up map ---
    fig, ax = plt.subplots(figsize=(12, 7),
                           subplot_kw={'projection': ccrs.PlateCarree()})

    # --- Plot grid ---
    grid_mesh = ax.pcolormesh(grid_lon, grid_lat, grid_values,
                              cmap=cmap, shading='auto',
                              norm=norm, transform=ccrs.PlateCarree())

    # --- Plot valid points ---
    ax.scatter(point_lons_valid, point_lats_valid, c=point_values_valid,
               cmap=cmap, norm=norm, s=point_size,
               edgecolor='black', transform=ccrs.PlateCarree(), zorder=3)

    # --- Plot NaN points with special styling ---
    ax.scatter(point_lons_nan, point_lats_nan, c=nan_point_color,
               s=nan_point_size, edgecolor=nan_point_edge,
               transform=ccrs.PlateCarree(), zorder=3)

    # --- Map features ---
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # --- Gridlines ---
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # --- Set extent with margin ---
    lon_min = min(np.nanmin(grid_lon), np.nanmin(point_lons))
    lon_max = max(np.nanmax(grid_lon), np.nanmax(point_lons))
    lat_min = min(np.nanmin(grid_lat), np.nanmin(point_lats))
    lat_max = max(np.nanmax(grid_lat), np.nanmax(point_lats))

    ax.set_extent([lon_min - margin, lon_max + margin,
                   lat_min - margin, lat_max + margin], crs=ccrs.PlateCarree())

    # --- Shared colorbar ---
    cbar = plt.colorbar(grid_mesh, ax=ax, orientation='vertical', shrink=0.8, label="Values")

    # --- Title and show ---
    ax.set_title(title)
    plt.tight_layout()
    plt.show(block=True)
