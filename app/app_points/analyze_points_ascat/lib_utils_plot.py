# lib_utils_plot.py
"""
Plotting routines for point and grid data using GeoPandas and Matplotlib (no Basemap).
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import matplotlib.pyplot as plt
import geopandas as gpd

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# Load once at import
_WORLD = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to plot data points (over a map)
def plot_data_points(
    lons: np.ndarray,
    lats: np.ndarray,
    sm: np.ndarray,
    ax=None,
    world_gdf: gpd.GeoDataFrame = _WORLD,
    vmin: float = 0.0,
    vmax: float = 100.0,
    cmap: str = "viridis",
    marker_size: int = 30,
    title: str = "Soil Moisture (%)",
    zoom_to_data: bool = True,
    output_path: str = None
) -> None:
    """
    Scatter-plot points colored by SM value on a GeoPandas world map.
    """
    # Mask and clip SM
    sm = np.ma.masked_invalid(sm)
    sm = np.clip(sm, vmin, vmax)

    # Prepare axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    # Plot base map
    world_gdf.plot(ax=ax, color='lightgray', edgecolor='white')

    # Scatter points
    sc = ax.scatter(
        lons, lats,
        c=sm,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        s=marker_size,
        marker='o',       # ensure filled markers
        linewidth=0
    )

    # Zoom to data extent
    if zoom_to_data:
        buffer_deg = 0.5  # half-degree buffer
        minx, maxx = lons.min(), lons.max()
        miny, maxy = lats.min(), lats.max()
        ax.set_xlim(minx - buffer_deg, maxx + buffer_deg)
        ax.set_ylim(miny - buffer_deg, maxy + buffer_deg)

    # Colorbar and labels
    cb = plt.colorbar(sc, ax=ax, label=title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(title)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300)
    else:
        plt.show()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to plot data grid (over a map)
def plot_data_grid(
        xi: np.ndarray,
        yi: np.ndarray,
        grid_data: np.ndarray,
        ax: plt.Axes = None,
        world_gdf: gpd.GeoDataFrame = _WORLD,
        output_path: str = None,
        cmap: str = "viridis",
        vmin: float = None,
        vmax: float = None,
        buffer_deg: float = 0.5,
        equal_aspect: bool = True,
) -> None:
    """
    Plot interpolated grid using pcolormesh on a GeoPandas map.

    Parameters
    ----------
    xi, yi : 1D arrays of longitudes and latitudes defining the grid
    grid_data : 2D array of values (shape len(yi)×len(xi))
    ax : existing Matplotlib Axes (optional)
    world_gdf : GeoDataFrame for base map (optional)
    cmap : colormap name
    vmin, vmax : color limits (auto-inferred if None)
    buffer_deg : degrees to pad around grid for plot limits
    equal_aspect : if True, sets aspect='equal'
    output_path : path to save PNG (if provided)
    """
    # Mask invalid data
    data = np.ma.masked_invalid(grid_data)

    # Auto infer vmin/vmax
    if vmin is None:
        vmin = np.ma.min(data)
    if vmax is None:
        vmax = np.ma.max(data)

    # Prepare axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Plot base map
    world_gdf.plot(ax=ax, color='lightgray', edgecolor='white')

    # Build mesh
    lon_grid, lat_grid = np.meshgrid(xi, yi)

    # Plot pcolormesh
    mesh = ax.pcolormesh(
        lon_grid, lat_grid, data,
        cmap=cmap, vmin=vmin, vmax=vmax,
        shading='auto'
    )

    # Zoom to grid extent
    minx, maxx = xi.min(), xi.max()
    miny, maxy = yi.min(), yi.max()
    ax.set_xlim(minx - buffer_deg, maxx + buffer_deg)
    ax.set_ylim(miny - buffer_deg, maxy + buffer_deg)

    # Equal aspect if desired
    if equal_aspect:
        ax.set_aspect('equal', adjustable='box')

    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, label='Interpolated SM')

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Interpolated Soil Moisture Grid')
    plt.tight_layout()

    # Save or show
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved grid plot to {output_path}")
    else:
        plt.show()
# ----------------------------------------------------------------------------------------------------------------------
