#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RECOLOUR TOOLS - SOILGRIDS IMG2CELL - REprocess paCkage for static soil datasets

__date__ = '20260521'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_img2cell_soilgrids.py -settings_file configuration.json

Version(s):
20260521 (1.0.0) --> First development for static SoilGrids datasets
"""
#-----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import re
import time
import json
import glob
import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime, timezone
import h5py
import numpy as np
import rasterio

from lib_utils_grid import GridRegistry

# logger
alg_logger = logging.getLogger("app_soilgrids_tiff2cell")
#-----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# parameters
FILL_I32 = np.int32(-2147483648)
# algorithm information
PROJECT_NAME = 'recolour'
ALG_NAME = 'img2cell_soil_grids'
ALG_TYPE = 'Application'
ALG_VERSION = '1.0.0'
ALG_RELEASE = '2026-05-21'
# -------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
## UTILS SETUP
# helper to get arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="Convert static SoilGrids GeoTIFFs to WARP5 cell NetCDF files"
    )
    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file"
    )
    return parser.parse_args()

# helper to setup logger
def setup_logger(settings):
    log_settings = settings.get("logging", {})
    level_name = str(log_settings.get("level", "INFO")).upper()
    level = getattr(logging, level_name, logging.INFO)

    alg_logger.handlers = []
    alg_logger.setLevel(level)
    alg_logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    alg_logger.addHandler(console_handler)

    folder = log_settings.get("folder")
    filename = log_settings.get("filename")
    if folder and filename:
        folder = fill_template(folder)
        filename = fill_template(filename)
        make_folder(folder)
        file_handler = logging.FileHandler(Path(folder) / filename)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        alg_logger.addHandler(file_handler)
#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
## UTILS I/O
# helper to read file json
def read_file_json(file_name):
    if not os.path.exists(file_name):
        raise FileNotFoundError(f'File "{file_name}" not found')
    with open(file_name, "r", encoding="utf-8") as file_handle:
        return json.load(file_handle)

# method to read grid
def load_cell_grid(grid_name="fibgrid_6.25"):
    registry = GridRegistry()
    grid = registry.get(grid_name)
    return grid

# helper to write cell
def write_cell_file(output_file, cell_id, gpi, lon, lat, soil_data, layer_meta, attrs=None, coordinate_scale=1e-6):
    make_folder(os.path.dirname(output_file))
    n_points = len(gpi)

    lon_i32 = np.round(lon / coordinate_scale).astype(np.int32)
    lat_i32 = np.round(lat / coordinate_scale).astype(np.int32)

    with h5py.File(output_file, "w") as dst:
        base_attrs = {
            "featureType": "point",
            "grid_mapping_name": "WARP5",
            "product_type": "static",
            "cell_id": int(cell_id),
            "n_locations": int(n_points),
        }
        if attrs:
            base_attrs.update(attrs)
        for key, val in base_attrs.items():
            dst.attrs[key] = _attr_value(val)

        #dst.create_dataset("obs", data=np.arange(n_points, dtype=np.int32))

        # define real dimension: locations
        locations = dst.create_dataset(
            "locations",
            data=np.arange(n_points, dtype=np.int32)
        )
        locations.make_scale("locations")
        locations.attrs["long_name"] = _attr_value("locations dimension")

        add_1d_dataset(
            dst, "location_id", gpi, np.int32,
            attrs={
                "long_name": "Location identifier (Grid Point ID)",
                "cf_role": "timeseries_id",
                "coordinates": "latitude longitude",
            }
        )
        lon_dset = add_1d_dataset(
            dst, "longitude", lon_i32, np.int32,
            attrs={
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "axis": "X",
                "missing_value": int(FILL_I32),
            },
            fillvalue=FILL_I32
        )
        lon_dset.attrs["scale_factor"] = np.float64(coordinate_scale)
        lon_dset.attrs["add_offset"] = np.float64(0.0)

        lat_dset = add_1d_dataset(
            dst, "latitude", lat_i32, np.int32,
            attrs={
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
                "axis": "Y",
                "missing_value": int(FILL_I32),
            },
            fillvalue=FILL_I32
        )
        lat_dset.attrs["scale_factor"] = np.float64(coordinate_scale)
        lat_dset.attrs["add_offset"] = np.float64(0.0)

        for var_name, values in soil_data.items():
            meta = layer_meta.get(var_name, {})
            dset = add_1d_dataset(
                dst, var_name, values, np.float32,
                attrs={
                    "long_name": meta.get("long_name", var_name),
                    "units": meta.get("units", ""),
                    "coordinates": "latitude longitude",
                    "source": meta.get("source", "SoilGrids GeoTIFF sampled at WARP5 GPI coordinates"),
                    "missing_value": "nan",
                },
                fillvalue=None
            )
            dset.attrs["_FillValue"] = np.float32(np.nan)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# UTILS BASE
# helper to make folder
def make_folder(folder_name):
    if folder_name:
        os.makedirs(folder_name, exist_ok=True)

# helper to fill template
def fill_template(path_template, tags=None):
    if tags is None:
        tags = {}
    out = str(path_template)
    now = datetime.now(timezone.utc)
    defaults = {
        "yyyy": now.strftime("%Y"),
        "mm": now.strftime("%m"),
        "dd": now.strftime("%d"),
        "date": now.strftime("%Y%m%d"),
    }
    defaults.update(tags)
    for key, value in defaults.items():
        out = out.replace("{" + key + "}", str(value))
    return out

# helper to clean variable name
def clean_var_name(file_path):
    name = Path(file_path).stem
    name = name.replace("(1)", "")
    name = re.sub(r"[^0-9a-zA-Z_]+", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    return name

# method to collect tiff
def collect_tiff(settings):
    source = settings.get("source", {})

    layers = source.get("layers")
    if layers:
        out = []
        for layer in layers:
            file_name = layer.get("file")
            if file_name is None:
                raise RuntimeError("Each source.layers item needs a file field")
            var_name = layer.get("var_name", clean_var_name(file_name))
            out.append({
                "file": fill_template(file_name),
                "var_name": var_name,
                "long_name": layer.get("long_name", var_name),
                "units": layer.get("units", ""),
                "scale_factor": layer.get("scale_factor"),
                "add_offset": layer.get("add_offset", 0.0),
                "valid_min": layer.get("valid_min"),
                "valid_max": layer.get("valid_max"),
            })
        return out

    folder = source.get("folder")
    pattern = source.get("file_pattern", "*.tif")
    if folder is None:
        raise RuntimeError("source.folder is required when source.layers is not provided")

    folder = fill_template(folder)
    files = sorted(glob.glob(str(Path(folder) / pattern)))
    if not files:
        raise RuntimeError(f"No TIFF files found using {Path(folder) / pattern}")

    return [{"file": f, "var_name": clean_var_name(f), "long_name": clean_var_name(f), "units": ""} for f in files]

# method to sample tiff on points
def sample_tiff_on_points(layer, lon, lat, require_epsg4326=True, invalid_below=-30000):
    tif_file = layer["file"]
    with rasterio.open(tif_file) as src:
        if require_epsg4326 and src.crs is not None and src.crs.to_epsg() != 4326:
            raise RuntimeError(f"{tif_file} is not EPSG:4326. Reproject first or disable require_epsg4326.")

        bounds = src.bounds
        inside = (
            (lon >= bounds.left) & (lon <= bounds.right) &
            (lat >= bounds.bottom) & (lat <= bounds.top)
        )

        out = np.full(lon.shape, np.nan, dtype=np.float32)
        if not np.any(inside):
            return out, inside

        coords = list(zip(lon[inside], lat[inside]))
        values = np.asarray([v[0] for v in src.sample(coords)], dtype=np.float32)

        if src.nodata is not None:
            values[values == np.float32(src.nodata)] = np.nan
        if invalid_below is not None:
            values[values <= np.float32(invalid_below)] = np.nan

        scale_factor = layer.get("scale_factor")
        add_offset = layer.get("add_offset", 0.0)
        if scale_factor is not None:
            values = values * np.float32(scale_factor) + np.float32(add_offset)

        valid_min = layer.get("valid_min")
        valid_max = layer.get("valid_max")
        if valid_min is not None:
            values[values < np.float32(valid_min)] = np.nan
        if valid_max is not None:
            values[values > np.float32(valid_max)] = np.nan

        out[inside] = values

    return out, inside

# define attribute value
def _attr_value(value):
    if isinstance(value, bytes):
        return np.bytes_(value)
    if isinstance(value, str):
        return np.bytes_(value.encode("utf-8"))
    if isinstance(value, (int, np.integer)):
        return np.int32(value)
    if isinstance(value, (float, np.floating)):
        return np.float32(value)
    return np.bytes_(str(value).encode("utf-8"))

# helper to add 1d dataset
def add_1d_dataset(dst, name, values, dtype, attrs=None, fillvalue=None, compression=True):
    data = values.astype(dtype)
    kwargs = {}
    if compression:
        kwargs.update({"compression": "gzip", "compression_opts": 4, "shuffle": True})
    if fillvalue is not None:
        kwargs["fillvalue"] = fillvalue

    dset = dst.create_dataset(name, data=data, **kwargs)

    if "locations" in dst:
        dset.dims[0].attach_scale(dst["locations"])
        dset.dims[0].label = "locations"

    if attrs:
        for key, val in attrs.items():
            if val is not None:
                dset.attrs[key] = _attr_value(val)
    return dset
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
## UTILS VARIABLES
# helper to convert bdod units
def convert_bdod_to_g_cm3(bdod, units):
    units = str(units).lower().replace(" ", "")
    if units in ["cgcm-3", "cg/cm3", "cgcm^-3"]:
        return bdod * np.float32(0.01)
    if units in ["kgdm-3", "kg/dm3", "gcm-3", "g/cm3"]:
        return bdod.astype(np.float32)
    raise RuntimeError(f'Unsupported bdod units "{units}"')

# helper to compute porosity
def compute_porosity_from_bdod(bdod, bdod_units="cg cm-3", particle_density=2.65):
    bdod_g_cm3 = convert_bdod_to_g_cm3(bdod.astype(np.float32), bdod_units)
    porosity = 1.0 - bdod_g_cm3 / np.float32(particle_density)
    porosity = np.where((porosity > 0.05) & (porosity < 0.90), porosity, np.nan)
    return porosity.astype(np.float32)

# helprt to compute porosity weighted values
def compute_weighted_average(data_map, var_names, weights):
    weights = np.asarray(weights, dtype=np.float32)
    values = [data_map[name].astype(np.float32) for name in var_names]

    numerator = np.zeros(values[0].shape, dtype=np.float32)
    denominator = np.zeros(values[0].shape, dtype=np.float32)

    for value, weight in zip(values, weights):
        valid = np.isfinite(value)
        numerator[valid] += value[valid] * weight
        denominator[valid] += weight

    out = numerator / denominator
    out[denominator <= 0] = np.nan
    return out.astype(np.float32)

# method to compute extra variables
def compute_extra_variables(settings, soil_all, layer_meta):
    extra = settings.get("extra_variables", {})
    if not bool(extra.get("enabled", False)):
        return soil_all, layer_meta

    # iterate over variables
    for item in extra.get("variables", []):
        var_name = item["var_name"]
        method = item["method"]

        # info start
        alg_logger.info(f' -----> Layer {var_name} ... ')

        # select computing method
        if method == "porosity_from_bdod":
            bdod_name = item["bdod"]
            if bdod_name not in soil_all:
                raise RuntimeError(f'Variable "{bdod_name}" needed by "{var_name}" not found')

            soil_all[var_name] = compute_porosity_from_bdod(
                bdod=soil_all[bdod_name],
                bdod_units=item.get("bdod_units", "cg cm-3"),
                particle_density=float(item.get("particle_density", 2.65))
            )

        elif method == "weighted_average":
            var_names = item["variables"]
            weights = item["weights"]

            for name in var_names:
                if name not in soil_all:
                    raise RuntimeError(f'Variable "{name}" needed by "{var_name}" not found')

            soil_all[var_name] = compute_weighted_average(
                data_map=soil_all,
                var_names=var_names,
                weights=weights
            )

        else:
            # method not supported
            raise RuntimeError(f'Unsupported extra variable method "{method}"')

        # store layer(s)
        layer_meta[var_name] = {
            "long_name": item.get("long_name", var_name),
            "units": item.get("units", ""),
            "source": item.get("source", f'derived using method "{method}"')
        }

        # info end
        alg_logger.info(f' -----> Layer {var_name} ... DONE')

    return soil_all, layer_meta
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# helper to convert soil grids layers to cells
def convert_grid2cell(settings):

    # info start
    alg_logger.info(' ---> Convert grids to cell ... ')

    # get settings
    grid_settings = settings.get("grid", {})
    source_settings = settings.get("source", {})
    destination_settings = settings.get("destination", {})
    parameters = settings.get("parameters", {})
    output_settings = settings.get("output", {})

    grid_settings = settings.get("grid", {})
    grid_name = grid_settings.get("name", "fibgrid_6.25")
    #if grid_file is None:
    #    raise RuntimeError("grid.file or grid.reference_grid is required")
    #grid_file = fill_template(grid_file)

    #variables = grid_settings.get("variables", {})
    #use_land_flag = bool(grid_settings.get("use_land_flag", True))
    require_epsg4326 = bool(source_settings.get("require_epsg4326", True))
    invalid_below = source_settings.get("invalid_below", -30000)

    output_folder = destination_settings.get("folder")
    if output_folder is None:
        raise RuntimeError("destination.folder is required")
    output_folder = fill_template(output_folder)
    filename_template = destination_settings.get("filename", "soilgrids_cell_{cell:04d}.nc")

    keep_only_points_with_data = bool(parameters.get("keep_only_points_with_data", True))
    allowed_cells = parameters.get("cells", [])
    coordinate_scale = float(output_settings.get("coordinate_scale_factor", 1e-6))

    # get grid
    alg_logger.info(f' ----> Get reference grid "{grid_name}" ... ')

    grid = load_cell_grid(grid_name)

    lon = grid.activearrlon.astype(np.float64)
    lat = grid.activearrlat.astype(np.float64)
    gpi = grid.activegpis.astype(np.int32)
    cell = grid.activearrcell.astype(np.int32)

    alg_logger.info(f' ----> Get reference grid "{grid_name}" ... DONE')

    # collect layers
    alg_logger.info(' ----> Collect datasets ... ')
    layers = collect_tiff(settings)
    alg_logger.info(' ----> Collect datasets ... DONE')

    # iterate over layers
    alg_logger.info(' ----> Read datasets ... ')
    soil_all, layer_meta = {}, {}
    inside_any = np.zeros(gpi.shape, dtype=bool)
    for layer in layers:

        # layer name start
        alg_logger.info(f' -----> Layer {layer} ... ')

        tif_file = layer["file"]
        if not os.path.exists(tif_file):
            raise FileNotFoundError(f'TIFF file "{tif_file}" not found')
        var_name = layer["var_name"]

        alg_logger.info(f" ::: Read -- variable: {var_name} - file: {tif_file}")

        values, inside = sample_tiff_on_points(
            layer=layer,
            lon=lon,
            lat=lat,
            require_epsg4326=require_epsg4326,
            invalid_below=invalid_below,
        )
        soil_all[var_name] = values
        layer_meta[var_name] = {
            "long_name": layer.get("long_name", var_name),
            "units": layer.get("units", ""),
            "source": Path(tif_file).name,
        }
        inside_any |= inside & np.isfinite(values)

        # layer name end
        alg_logger.info(f' -----> Layer {layer} ... DONE')

    # layer group end
    alg_logger.info(' ----> Read datasets ... DONE')

    # layer group end
    alg_logger.info(' ----> Compute datasets ... ')
    soil_all, layer_meta = compute_extra_variables(
        settings=settings,
        soil_all=soil_all,
        layer_meta=layer_meta
    )
    alg_logger.info(' ----> Compute datasets ... DONE')

    # merge layers start
    alg_logger.info(' ----> Merge datasets ... ')

    # define valid points
    valid = np.ones(gpi.shape, dtype=bool)
    if keep_only_points_with_data:
        valid &= inside_any
    if allowed_cells:
        allowed_cells = np.asarray([int(c) for c in allowed_cells], dtype=np.int32)
        valid &= np.isin(cell, allowed_cells)

    lon = lon[valid]
    lat = lat[valid]
    gpi = gpi[valid]
    cell = cell[valid]
    soil_all = {name: values[valid] for name, values in soil_all.items()}

    global_attrs = settings.get("attributes", {})
    global_attrs.update({
        "source_grid": grid_name,
        "source_tiffs": ",".join(Path(layer["file"]).name for layer in layers),
        "rule_common_points": "static raster sample at WARP5 GPI coordinates",
    })

    # merge layers end
    alg_logger.info(' ----> Merge datasets ... DONE')

    # info dump start
    alg_logger.info(' ----> Dump datasets ... ')
    written = []
    make_folder(output_folder)
    # iterate over cells
    for cell_id in np.unique(cell):

        # info cell start
        alg_logger.info(f' -----> Cell {cell_id} ... ')

        idx = cell == cell_id
        order = np.argsort(gpi[idx])

        gpi_cell = gpi[idx][order]
        lon_cell = lon[idx][order]
        lat_cell = lat[idx][order]
        soil_cell = {name: values[idx][order] for name, values in soil_all.items()}

        # define output file
        output_file = os.path.join(
            output_folder,
            filename_template.format(cell=int(cell_id), cell_n=str(int(cell_id)).zfill(4))
        )
        alg_logger.info(f" ::: Write -- cell: {int(cell_id)} - file: {output_file}")

        # write cell file
        write_cell_file(
            output_file=output_file,
            cell_id=int(cell_id),
            gpi=gpi_cell,
            lon=lon_cell,
            lat=lat_cell,
            soil_data=soil_cell,
            layer_meta=layer_meta,
            attrs=global_attrs,
            coordinate_scale=coordinate_scale,
        )
        written.append(output_file)

        # info cell end
        alg_logger.info(f' -----> Cell {cell_id} ... DONE')

    # info dump end
    alg_logger.info(' ----> Dump datasets ... DONE')

    # info end
    alg_logger.info(' ---> Convert grids to cell ... DONE')

    return written
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get algorithm args
    algorithm_args = get_args()
    # get algorithm settings
    algorithm_settings = read_file_json(algorithm_args.settings_file)
    # set logger
    setup_logger(algorithm_settings)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    # time algorithm
    start_time = time.time()
    #-------------------------------------------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------
    # method to convert ce
    algorithm_info = convert_grid2cell(algorithm_settings)
    #-------------------------------------------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------
    # info algorithm (end)
    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')
    #-------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# wrapper main
if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
#-----------------------------------------------------------------------------------------------------------------------
