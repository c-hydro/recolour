# libraries
import logging
import geopandas as gpd

from shapely.geometry import Point, box
import numpy as np
import geopandas as gpd

# method to convert meters to degrees for bounding box buffer
def convert_meters_to_degrees(mask_lons, mask_lats, buffer_km=50):

    # Base bounding box
    min_lon, max_lon = np.min(mask_lons), np.max(mask_lons)
    min_lat, max_lat = np.min(mask_lats), np.max(mask_lats)

    # Convert km buffer to degrees
    deg_per_km_lat = 1 / 111.0  # 1° latitude ≈ 111 km
    mean_lat = (min_lat + max_lat) / 2.0
    deg_per_km_lon = 1 / (111.0 * np.cos(np.radians(mean_lat)))  # longitude varies with latitude

    buffer_deg_lat = buffer_km * deg_per_km_lat
    buffer_deg_lon = buffer_km * deg_per_km_lon

    return buffer_deg_lon, buffer_deg_lat

# method to compute coverage from coarse and fine grid arrays
def compute_coverage_from_arrays(coarse_lons, coarse_lats, fine_lons, fine_lats, fine_mask, coarse_km=12.5):
    """
    For each coarse grid point (2D), count how many land fine grid points fall within its box.
    """
    half_deg = coarse_km / 111.32 / 2  # approx conversion km to degrees
    logging.info(f" ::: Using ±{half_deg:.5f} degrees box per coarse grid center.")

    # Flatten inputs
    coarse_lats_flat = coarse_lats.ravel()
    coarse_lons_flat = coarse_lons.ravel()
    fine_lats_flat = fine_lats.ravel()
    fine_lons_flat = fine_lons.ravel()
    fine_mask_flat = fine_mask.ravel()

    # Stack fine grid points (lon, lat)
    fine_coords = np.column_stack((fine_lons_flat, fine_lats_flat))

    land_counts, total_counts = [], []
    for clat, clon in zip(coarse_lats_flat, coarse_lons_flat):
        # Bounding box in degrees
        lat_min, lat_max = clat - half_deg, clat + half_deg
        lon_min, lon_max = clon - half_deg, clon + half_deg

        # Find fine grid points within the box
        in_box = (
            (fine_coords[:, 0] >= lon_min) & (fine_coords[:, 0] <= lon_max) &
            (fine_coords[:, 1] >= lat_min) & (fine_coords[:, 1] <= lat_max)
        )

        total = np.sum(in_box)
        land = np.sum(fine_mask_flat[in_box]) if total > 0 else 0

        land_counts.append(land)
        total_counts.append(total)

    percent_land = [lc / tc * 100 if tc > 0 else 0 for lc, tc in zip(land_counts, total_counts)]

    return (
        np.array(land_counts).reshape(coarse_lats.shape),
        np.array(total_counts).reshape(coarse_lats.shape),
        np.array(percent_land).reshape(coarse_lats.shape)
    )

# method to generate a binary land mask for points using Natural Earth polygons
def get_land_mask(lons, lats):
    """
    Generate a binary land mask for points using Natural Earth polygons.
    """
    logging.info(" ::: Loading Natural Earth land boundaries ... ")
    land = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')).to_crs(epsg=4326)
    land_union = land.unary_union
    logging.info(" ::: Loading Natural Earth land boundaries ... DONE")

    logging.info(" ::: Generating land-sea mask ... ")
    points = gpd.GeoSeries([Point(lon, lat) for lon, lat in zip(lons, lats)], crs='EPSG:4326')
    mask = points.within(land_union).astype(int)
    logging.info(" ::: Generating land-sea mask ... DONE")

    return mask.values  # 1 for land, 0 for sea
