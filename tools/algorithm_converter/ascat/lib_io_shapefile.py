# libraries
import os
import numpy as np
import geopandas as gpd
import geodatasets

from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling, calculate_default_transform

import matplotlib.pyplot as plt

try:
    from shapely.validation import make_valid
except ImportError:
    def make_valid(geom): return geom.buffer(0)

from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

from lib_utils_io import compute_coverage_from_arrays

# private method to load a global land polygon
def _load_land_gdf(target_crs):
    """
    Load a global 'land' polygon using the same Natural Earth data behind
    Matplotlib/Cartopy coastlines. Prefer Cartopy 10m for detail.
    """
    try:
        # Preferred: Cartopy natural earth "physical:land"
        from cartopy.io import shapereader as shpreader
        land_path = shpreader.natural_earth(resolution='10m', category='physical', name='land')
        land = gpd.read_file(land_path)
    except Exception:
        # Fallback: dissolve countries into a land mass (less detailed coastline)
        land = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        # Drop Antarctica to avoid huge buffers around the ice shelf if undesired
        land = land[land['name'] != 'Antarctica']
        land = land[['geometry']]

    # Clean & dissolve to one polygon
    land = land[land.geometry.notnull() & (~land.geometry.is_empty)]

    # If CRS is missing, assign it (no reprojection yet)
    if land.crs is None:
        land = land.set_crs("EPSG:4326", allow_override=True)

    land = land.to_crs(target_crs)
    land_union = unary_union(land.geometry.values)
    # Fix any minor topology glitches
    land_union = land_union.buffer(0)
    return gpd.GeoDataFrame(geometry=[land_union], crs=target_crs)

# method to create terrain and sea buffers
def terrain_and_sea_buffers(
    gdf_in: gpd.GeoDataFrame,
    land_buffer_km: float = 50,
    sea_buffer_km: float = 20,
    projected_crs: str | None = None,   # choose a meters CRS; default below
    keep_original: bool = True,
    land_gdf: gpd.GeoDataFrame | None = None,
):
    """
    Create two outward buffers:
      - Land buffer: clipped to land only
      - Sea buffer: 0..sea_buffer_km ring outside the domain, clipped to water (not land)

    Returns a single-row GeoDataFrame with the union of:
      [optional original] ∪ land_buffer ∪ sea_buffer
    in the SAME CRS as input.
    """
    if gdf_in.empty:
        raise ValueError("Input GeoDataFrame is empty.")
    if gdf_in.crs is None:
        raise ValueError("Input GeoDataFrame has no CRS. Set gdf.crs first.")
    if land_buffer_km < 0 or sea_buffer_km < 0:
        raise ValueError("Buffers must be >= 0 km.")

    in_crs = gdf_in.crs
    proj_crs = projected_crs or "EPSG:6933"  # equal-area, meters (good global default)

    # Clean and project input
    gdf = gdf_in[gdf_in.geometry.notnull() & ~gdf_in.geometry.is_empty].copy()
    gdf["geometry"] = gdf.geometry.apply(make_valid)
    gdf_m = gdf.to_crs(proj_crs)

    # Load/prepare land (dissolved) in same CRS
    if land_gdf is None:
        land = _load_land_gdf(proj_crs)  # should return polygons covering land (dissolve not assumed)
    else:
        land = land_gdf.to_crs(proj_crs) if land_gdf.crs != proj_crs else land_gdf
    land = land[land.geometry.notnull() & ~land.geometry.is_empty].copy()
    land_union = unary_union(land.geometry.values)
    if land_union.is_empty:
        raise ValueError("Land mask is empty after dissolve/union.")

    # Unite domain
    dom_union = unary_union(gdf_m.geometry.values)
    dom_union = make_valid(dom_union)
    if dom_union.is_empty:
        raise ValueError("Input geometry became empty after union/clean.")

    m_per_km = 1000.0
    L = float(land_buffer_km) * m_per_km
    S = float(sea_buffer_km) * m_per_km

    # 1) Land buffer (clip buffered domain to land)
    land_buf = dom_union.buffer(L) if L > 0 else dom_union
    land_buf = land_buf.intersection(land_union)
    land_buf = make_valid(land_buf)

    # 2) Sea buffer: the *ring* outside the domain up to S, minus land
    #    - ring = buffer(S) minus the original domain (so "0..S km outward")
    #    - then remove land to keep only water portions
    if S > 0:
        sea_ring = dom_union.buffer(S).difference(dom_union)
        sea_buf = make_valid(sea_ring.difference(land_union))
    else:
        sea_buf = Polygon()  # empty

    # Optionally include the original domain
    pieces = []
    if keep_original:
        pieces.append(dom_union)
    if not land_buf.is_empty:
        pieces.append(land_buf)
    if not sea_buf.is_empty:
        pieces.append(sea_buf)

    if not pieces:
        raise ValueError("Resulting geometry is empty (check buffers and inputs).")

    out_geom = unary_union(pieces)
    out_geom = make_valid(out_geom)
    if out_geom.is_empty:
        raise ValueError("Unioned output is empty after validity fix.")

    out = gpd.GeoDataFrame({"geometry": [out_geom]}, crs=proj_crs).to_crs(in_crs)
    return out

# method to create terrain buffer
def terrain_buffers(gdf_in, buffer_km=50, projected_crs="EPSG:3857", keep_original=True):
    """
    Extend domain outward by buffer_km but only over terrain (no sea).
    Steps:
      1) Reproject to meters CRS
      2) Buffer outward
      3) Intersect with global land polygon
      4) (optional) union back with original domain
    Returns GeoDataFrame in the *same CRS as input*.
    """
    if gdf_in.empty:
        raise ValueError("Input GeoDataFrame is empty.")

    # Remember input CRS (assume set); reproject to meters for buffering
    in_crs = gdf_in.crs
    if in_crs is None:
        raise ValueError("Input GeoDataFrame has no CRS. Set gdf.crs first.")

    gdf_m = gdf_in.to_crs(projected_crs)
    gdf_m = gdf_m[~gdf_m.geometry.is_empty & gdf_m.geometry.notnull()]

    # Load/dissolve global land mask in same projected CRS
    land_gdf = _load_land_gdf(projected_crs)

    # Buffer outward (meters)
    buf_m = float(buffer_km) * 1000.0
    # buffer(0) to clean, then buffer out; dissolve first to avoid self overlaps
    dom_union = unary_union(gdf_m.geometry.values).buffer(0)
    dom_buffered = dom_union.buffer(buf_m)

    # Clip to land to avoid extending into the sea
    extended_on_land = dom_buffered.intersection(land_gdf.geometry.iloc[0]).buffer(0)

    # Optionally add the original domain back in
    if keep_original:
        extended_on_land = unary_union([extended_on_land, dom_union]).buffer(0)

    # Pack back into a GeoDataFrame and reproject to original CRS
    # Ensure geometry is MultiPolygon/Polygon as needed
    if extended_on_land.is_empty:
        raise ValueError("Resulting extended area is empty after land clipping.")

    out = gpd.GeoDataFrame(geometry=[extended_on_land], crs=projected_crs).to_crs(in_crs)
    return out


# method to convert shapefile to mask
def shapefile_to_mask(shapefile_path_in, shapefile_path_out=None,
                      resolution_km_fine=1, resolution_km_coarse=12.5, buffer_km_fine=20,
                      extend_area_flag=True,
                      extend_area_buffer_km_land=50, extend_area_buffer_km_sea=0,
                      update=False, plot=False):

    # Convert shape resolution from km to meters
    resolution_meters_fine = resolution_km_fine * 1000
    # Convert buffer to meters
    buffer_meters_fine = buffer_km_fine * 1000

    # Read and reproject shapefile to EPSG:3857 for accurate rasterization
    tmp = gpd.read_file(shapefile_path_in).to_crs(epsg=3857)
    if tmp.empty:
        raise ValueError("The shapefile contains no valid geometries.")

    # update and remove existing shapefile if update flag is set
    if update:
        if os.path.exists(shapefile_path_out):
            os.remove(shapefile_path_out)

    # If extend_area_flag is True, compute the extended area
    if extend_area_flag:
        if not os.path.exists(shapefile_path_out):
            # compute the extended area
            gdf = terrain_and_sea_buffers(
                tmp,
                land_buffer_km=extend_area_buffer_km_land,
                sea_buffer_km=extend_area_buffer_km_sea)
            #gdf = terrain_buffers(tmp, buffer_km=extend_area_buffer_km)

            # save the extended area to a shapefile
            folder_path, _ = os.path.split(shapefile_path_out)
            os.makedirs(folder_path, exist_ok=True)
            gdf.to_file(shapefile_path_out)
        else:
            # load the existing extended area shapefile
            gdf = gpd.read_file(shapefile_path_out).to_crs(epsg=3857)
            if gdf.empty:
                raise ValueError("The extended area shapefile contains no valid geometries.")
    else:
        # Use the original shapefile without extension
        gdf = tmp

    # compute bbox and dimensions
    minx_fine, miny_fine, maxx_fine, maxy_fine = map(float, gdf.total_bounds)

    # Extend bounding box with buffer
    minx_fine = minx_fine - buffer_meters_fine
    maxx_fine = maxx_fine + buffer_meters_fine
    miny_fine = miny_fine - buffer_meters_fine
    maxy_fine = maxy_fine + buffer_meters_fine

    width_meters_fine = int(np.ceil((maxx_fine - minx_fine) / resolution_meters_fine))
    height_meters_fine = int(np.ceil((maxy_fine - miny_fine) / resolution_meters_fine))

    # compute transform meters
    transform_meters_fine = from_origin(minx_fine, maxy_fine, resolution_meters_fine, resolution_meters_fine)

    # rasterize the shapefile geometries
    mask_meters_fine = rasterize(
        [(geom, 1) for geom in gdf.geometry],
        out_shape=(height_meters_fine, width_meters_fine),
        transform=transform_meters_fine,
        fill=0,
        all_touched=False,
        dtype='uint8'
    )

    # reproject raster to EPSG:4326
    mask_meters_epsg, mask_deg_epsg = 'EPSG:3857', 'EPSG:4326'
    transform_deg_fine, width_deg_fine, height_deg_fine = calculate_default_transform(
      mask_meters_epsg, mask_deg_epsg,
        width_meters_fine, height_meters_fine, minx_fine, miny_fine, maxx_fine, maxy_fine
    )
    mask_deg_fine = np.zeros((height_deg_fine, width_deg_fine), dtype='uint8')

    # object to reproject the mask
    reproject(
        source=mask_meters_fine,
        destination=mask_deg_fine,
        src_transform=transform_meters_fine,
        src_crs=mask_meters_epsg,
        dst_transform=transform_deg_fine,
        dst_crs=mask_deg_epsg,
        resampling=Resampling.nearest
    )

    # optional plot
    if plot:
        plt.figure(figsize=(8, 6))
        plt.imshow(mask_deg_fine)
        plt.grid(True)
        plt.show(block=True)

    # create lon/lat grids for fine-resolution mask
    xs_fine = np.arange(width_deg_fine) * transform_deg_fine.a + transform_deg_fine.c + transform_deg_fine.a / 2
    ys_fine = np.arange(height_deg_fine) * transform_deg_fine.e + transform_deg_fine.f + transform_deg_fine.e / 2
    lons_fine, lats_fine = np.meshgrid(xs_fine, ys_fine)

    # regrid to resolution_km_coarse
    resolution_deg_coarse = resolution_km_coarse / 111.32  # Approximate conversion

    # get bounding box from current transform
    min_lon_fine = transform_deg_fine.c
    max_lat_fine = transform_deg_fine.f
    max_lon_fine = min_lon_fine + width_deg_fine * transform_deg_fine.a
    min_lat_fine = max_lat_fine + height_deg_fine * transform_deg_fine.e  # e is negative

    # compute new grid shape
    width_deg_coarse = int(np.ceil((max_lon_fine - min_lon_fine) / resolution_deg_coarse))
    height_deg_coarse = int(np.ceil((max_lat_fine - min_lat_fine) / resolution_deg_coarse))
    transform_deg_coarse = from_origin(min_lon_fine, max_lat_fine, resolution_deg_coarse, resolution_deg_coarse)

    # create new empty array and reproject
    mask_deg_coarse = np.zeros((height_deg_coarse, width_deg_coarse), dtype='uint8')
    reproject(
        source=mask_deg_fine,
        destination=mask_deg_coarse,
        src_transform=transform_deg_fine,
        src_crs=mask_deg_epsg,
        dst_transform=transform_deg_coarse,
        dst_crs=mask_deg_epsg,
        resampling=Resampling.nearest
    )

    # bbox in degrees for coarse grid
    bbox_deg_coarse = (
        transform_deg_coarse.c,
        transform_deg_coarse.f + height_deg_coarse * transform_deg_coarse.e,
        transform_deg_coarse.c + width_deg_coarse * transform_deg_coarse.a,
        transform_deg_coarse.f
    )

    # create lon/lat grids for coarse-resolution mask
    xs_coarse = np.arange(width_deg_coarse) * transform_deg_coarse.a + transform_deg_coarse.c + transform_deg_coarse.a / 2
    ys_coarse = np.arange(height_deg_coarse) * transform_deg_coarse.e + transform_deg_coarse.f + transform_deg_coarse.e / 2
    lons_coarse, lats_coarse = np.meshgrid(xs_coarse, ys_coarse)

    # Optional plot
    if plot:
        plt.figure(figsize=(8, 6))
        plt.imshow(mask_deg_coarse)
        plt.grid(True)
        plt.show(block=True)

    # define coarse grid over land
    land_coarse, count_coarse, percentage_coarse = compute_coverage_from_arrays(
        lons_coarse, lats_coarse, lons_fine, lats_fine,
        mask_deg_fine, coarse_km=resolution_km_coarse)

    # Optional plot
    if plot:
        plt.figure(figsize=(8, 6))
        plt.imshow(land_coarse)
        plt.grid(True)
        plt.figure(figsize=(8, 6))
        plt.imshow(count_coarse)
        plt.grid(True)
        plt.figure(figsize=(8, 6))
        plt.imshow(percentage_coarse)
        plt.grid(True)
        plt.show(block=True)

    return (mask_deg_coarse, transform_deg_coarse, bbox_deg_coarse,
            lons_coarse, lats_coarse, land_coarse, count_coarse, percentage_coarse)