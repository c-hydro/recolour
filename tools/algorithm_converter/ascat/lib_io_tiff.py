import logging
from rasterio.transform import from_bounds

import geopandas as gpd
import numpy as np

import rasterio

# method to save a GeoTIFF file from grid data
def write_tiff(
    path,
    lons,              # 2D lon grid (same shape as data)
    lats,              # 2D lat grid (same shape as data)
    grid_bands: dict,  # {"band_name": 2D numpy array}
    crs="EPSG:4326",
    nodata=np.nan,
    band_order=None    # optional: custom list like ["precip","temp"]
):
    if not grid_bands:
        raise ValueError("grid_bands is empty")

    # Decide band order (stable & explicit)
    band_names = band_order if band_order is not None else sorted(grid_bands.keys())

    # Validate shapes and get height/width
    first = grid_bands[band_names[0]]
    height, width = first.shape
    for k in band_names:
        if grid_bands[k].shape != (height, width):
            raise ValueError(f"Band '{k}' shape {grid_bands[k].shape} != {(height, width)}")

    # Build geotransform from the lon/lat bounds (assumes rectilinear grid)
    minx, maxx = float(np.nanmin(lons)), float(np.nanmax(lons))
    miny, maxy = float(np.nanmin(lats)), float(np.nanmax(lats))
    transform = from_bounds(minx, miny, maxx, maxy, width, height)

    # Rasterio wants row 0 at the top (north). If your lats increase downward, flip.
    flip_required = lats[0, 0] < lats[-1, 0]  # top row has smaller lat than bottom row
    dtype = np.result_type(*[grid_bands[k].dtype for k in band_names])

    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": len(band_names),
        "dtype": dtype,
        "crs": crs,
        "transform": transform,
        "nodata": nodata,
        "compress": "deflate",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(path, "w", **profile) as dst:
        for i, name in enumerate(band_names, start=1):
            data = grid_bands[name]
            if flip_required:
                data = np.flipud(data)
            dst.write(data.astype(dtype, copy=False), i)
            dst.set_band_description(i, name)      # label the band with dict key
            # Optional per-band metadata:
            # dst.update_tags(i, long_name=name)

        # Optional global tags
        dst.update_tags(source="ascat_to_tif", bands=",".join(band_names))
