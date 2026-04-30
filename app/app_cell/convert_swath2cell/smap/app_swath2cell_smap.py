#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SMAP SPL2SMP_E SWATH2CELL

Convert SMAP SPL2SMP_E swath files to cell NetCDF files using one of two grid modes:

1) grid_domain
   - Uses a raster reference grid, usually a domain mask GeoTIFF.
   - Computes gpi from raster row/col.
   - Splits the valid domain into configurable n_rows x n_cols cells.

2) grid_ascat
   - Uses an ASCAT grid NetCDF containing at least lon, lat, gpi, cell.
   - Assigns each SMAP point to the nearest ASCAT grid point.
   - Uses native ASCAT gpi and cell.

General command line:
python app_swath2cell_smap_grid_domain_ascat_flat.py \
    -settings_file app_swath2cell_smap.json \
    -time "YYYY-MM-DD HH:MM"

Example grid_domain JSON section:
"grid": {
  "mode": "grid_domain",
  "reference_grid": "/path/to/grid_reference_500m.tiff",
  "n_rows": 3,
  "n_cols": 3,
  "id_digits": 2
}

Example grid_ascat JSON section:
"grid": {
  "mode": "grid_ascat",
  "reference_grid": "/path/to/TUW_WARP5_grid_info_2_3.nc",
  "max_distance_km": 25,
  "use_land_flag": true,
  "output_smap_coordinates": false,
  "id_digits": 4
}
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import re
import sys
import json
import time
import logging
import argparse

from pathlib import Path
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from rasterio.transform import rowcol

try:
    from scipy.spatial import cKDTree
except ImportError:  # handled later only if grid_ascat is used
    cKDTree = None
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application for converting SMAP SPL2SMP_E swath to cells'
alg_type = 'Package'
alg_version = '1.2.0'
alg_release = '2026-04-29'

# algorithm globals
alg_logger = logging.getLogger('app_swath2cell')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# basic utils

def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    raise FileNotFoundError(f' ===> File "{file_name}" not found. Exit')


def get_args():
    parser = argparse.ArgumentParser(
        description='Convert SMAP SPL2SMP_E swath files to cell files using grid_domain or grid_ascat'
    )

    parser.add_argument(
        '-settings_file',
        dest='settings_file',
        required=True,
        help='Path to JSON settings file'
    )

    parser.add_argument(
        '-time',
        dest='time_run',
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM"'
    )

    return parser.parse_args()


def make_folder(folder_name):
    if folder_name is not None and folder_name != '':
        os.makedirs(folder_name, exist_ok=True)


def get_logger(settings, time_obj=None):
    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    alg_logger.handlers = []
    alg_logger.setLevel(level)
    alg_logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    alg_logger.addHandler(ch)

    log_folder_template = log_settings.get("folder", None)
    log_file_template = log_settings.get("filename", None)

    if log_folder_template and log_file_template and time_obj is not None:
        log_folder = fill_time_template(log_folder_template, time_obj)
        log_file = fill_time_template(log_file_template, time_obj)
        make_folder(log_folder)

        fh = logging.FileHandler(Path(log_folder) / log_file)
        fh.setLevel(level)
        fh.setFormatter(formatter)
        alg_logger.addHandler(fh)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# time utils

def parse_time_string(time_string):
    if time_string is None:
        return None
    return datetime.strptime(time_string, '%Y-%m-%d %H:%M')


def round_time_to_hour(time_obj):
    return time_obj.replace(minute=0, second=0, microsecond=0)


def round_time_to_midnight(time_obj):
    return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)


def get_reference_time(args):
    if args.time_run is not None:
        reference_time = parse_time_string(args.time_run)
    else:
        reference_time = datetime.now()
    return round_time_to_hour(reference_time)


def get_time_window(args, settings):
    time_settings = settings.get('time', {})

    period_mode = time_settings.get('period_mode', 'last_n_hours')
    time_start_raw = time_settings.get('time_start', None)
    time_end_raw = time_settings.get('time_end', None)

    if period_mode == 'time_range':
        if time_start_raw is None or time_end_raw is None:
            raise RuntimeError('time_start and time_end must be defined for period_mode time_range')

        date_start = round_time_to_hour(parse_time_string(time_start_raw))
        date_end = round_time_to_hour(parse_time_string(time_end_raw))

    else:
        if time_end_raw is None:
            date_end = get_reference_time(args)
        else:
            date_end = round_time_to_hour(parse_time_string(time_end_raw))

        if period_mode == 'last_n_hours':
            n_hours = int(time_settings.get('n_hours', 24))
            date_start = date_end - timedelta(hours=n_hours)

        elif period_mode == '00_12':
            date_start = round_time_to_midnight(date_end)
            date_end = date_start + timedelta(hours=12)

        elif period_mode == '12_24':
            date_start = round_time_to_midnight(date_end) + timedelta(hours=12)
            date_end = date_start + timedelta(hours=12)

        elif period_mode == 'daily':
            date_start = round_time_to_midnight(date_end)
            date_end = date_start + timedelta(hours=24)

        else:
            raise RuntimeError(f'Unsupported period_mode: {period_mode}')

        date_start = round_time_to_hour(date_start)
        date_end = round_time_to_hour(date_end)

    if date_start > date_end:
        raise RuntimeError('date_start must be less than or equal to date_end')

    return date_start, date_end


def hourly_range(date_start, date_end):
    time_step = date_start
    while time_step <= date_end:
        yield time_step
        time_step += timedelta(hours=1)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# template utils

def fill_time_template(path_template, time_obj):
    if path_template is None:
        return None

    if "{time_step" not in path_template:
        return path_template

    pattern = r"\{time_step:\s*(.*?)\}"

    def replace(match):
        time_format = match.group(1)
        return time_obj.strftime(time_format)

    return re.sub(pattern, replace, path_template)


def extract_datetime_from_filename(file_name, settings):
    source_settings = settings.get('source', {})

    file_regex = source_settings.get(
        'file_regex',
        r'^SMAP_L2_SM_P_E_[0-9]+_[AD]_(?P<time>[0-9]{8}T[0-9]{6})_R[0-9]+_[0-9]+\.h5$'
    )

    match = re.match(file_regex, file_name)
    if match is None:
        return None

    time_stamp = match.group('time')
    return datetime.strptime(time_stamp, "%Y%m%dT%H%M%S")


def define_cell_file(settings, cell_id, time_obj):
    destination_settings = settings.get('destination', {})
    grid_settings = settings.get('grid', {})

    folder_template = destination_settings['folder']
    file_template = destination_settings.get('filename', 'smap_{cell}.nc')

    cell_digits = int(grid_settings.get('id_digits', 2))
    cell_name = str(int(cell_id)).zfill(cell_digits)

    folder_name = fill_time_template(folder_template, time_obj)
    file_name = fill_time_template(file_template, time_obj)

    file_name = file_name.replace('{cell}', cell_name)
    file_name = file_name.replace('{cell_n}', cell_name)

    return Path(folder_name) / file_name
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# settings validation

def validate_settings(settings):
    grid_settings = settings.get('grid', {})
    mode = grid_settings.get('mode', 'grid_domain')

    if mode not in ['grid_domain', 'grid_ascat']:
        raise RuntimeError(f'Unsupported grid mode: {mode}. Supported modes: grid_domain, grid_ascat')

    if 'reference_grid' not in grid_settings:
        raise RuntimeError('grid.reference_grid must be defined')

    if mode == 'grid_domain':
        if 'n_rows' not in grid_settings or 'n_cols' not in grid_settings:
            raise RuntimeError('grid_domain mode requires grid.n_rows and grid.n_cols')

        n_rows = int(grid_settings.get('n_rows'))
        n_cols = int(grid_settings.get('n_cols'))

        if n_rows <= 0 or n_cols <= 0:
            raise RuntimeError('grid.n_rows and grid.n_cols must be positive integers')

    if mode == 'grid_ascat':
        if cKDTree is None:
            raise RuntimeError('grid_ascat mode requires scipy. Install it with: pip install scipy')

        max_distance_km = float(grid_settings.get('max_distance_km', 25))
        if max_distance_km <= 0:
            raise RuntimeError('grid.max_distance_km must be positive')

    return mode
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# file utils

def collect_smap_files(settings, date_start, date_end):
    source_settings = settings.get('source', {})

    folder_template = source_settings['folder']
    file_pattern = source_settings.get('filename', source_settings.get('file_pattern', '*.h5'))

    file_list = []

    for time_step in hourly_range(date_start, date_end):
        folder_name = fill_time_template(folder_template, time_step)
        folder_path = Path(folder_name)

        alg_logger.info(f' ----> Search source folder: {folder_path}')

        if not folder_path.exists():
            continue

        for file_path in folder_path.glob(file_pattern):
            file_time = extract_datetime_from_filename(file_path.name, settings)
            if file_time is None:
                continue

            file_time_hour = round_time_to_hour(file_time)
            if date_start <= file_time_hour <= date_end:
                file_list.append(file_path)

    return sorted(list(set(file_list)))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# smap utils

def read_smap_file(file_path, settings):
    dataset_settings = settings.get('datasets', {})

    lat_name = dataset_settings.get('latitude', '/Soil_Moisture_Retrieval_Data/latitude')
    lon_name = dataset_settings.get('longitude', '/Soil_Moisture_Retrieval_Data/longitude')
    sm_name = dataset_settings.get('soil_moisture', '/Soil_Moisture_Retrieval_Data/soil_moisture')

    with h5py.File(file_path, 'r') as file_handle:
        latitude = np.array(file_handle[lat_name])
        longitude = np.array(file_handle[lon_name])
        soil_moisture = np.array(file_handle[sm_name])

    latitude = latitude.astype(float)
    longitude = longitude.astype(float)
    soil_moisture = soil_moisture.astype(float)

    soil_moisture[soil_moisture < -9000] = np.nan

    return longitude, latitude, soil_moisture
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# grid_domain utils

def compute_domain_bounds(reference_data):
    valid_mask = np.isfinite(reference_data)
    valid_rows, valid_cols = np.where(valid_mask)

    if valid_rows.size == 0:
        raise RuntimeError('No valid cells found in reference grid')

    row_min = int(valid_rows.min())
    row_max = int(valid_rows.max())
    col_min = int(valid_cols.min())
    col_max = int(valid_cols.max())

    return row_min, row_max, col_min, col_max


def compute_gpi_from_row_col(row, col, n_cols_grid):
    return int(row) * int(n_cols_grid) + int(col)


def compute_domain_cell_id(row, col, domain_bounds, n_cell_rows, n_cell_cols):
    row_min, row_max, col_min, col_max = domain_bounds

    domain_height = row_max - row_min + 1
    domain_width = col_max - col_min + 1

    row_rel = int(row) - row_min
    col_rel = int(col) - col_min

    cell_height = float(domain_height) / float(n_cell_rows)
    cell_width = float(domain_width) / float(n_cell_cols)

    cell_row = int(float(row_rel) / cell_height)
    cell_col = int(float(col_rel) / cell_width)

    cell_row = min(max(cell_row, 0), n_cell_rows - 1)
    cell_col = min(max(cell_col, 0), n_cell_cols - 1)

    cell_id = cell_row * n_cell_cols + cell_col
    return int(cell_id)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# grid_ascat utils

def normalize_longitudes_to_grid(point_lon, grid_lon):
    """Normalize SMAP longitudes if ASCAT grid and SMAP use different longitude conventions."""
    point_lon = np.asarray(point_lon, dtype=float)
    grid_lon = np.asarray(grid_lon, dtype=float)

    grid_min = np.nanmin(grid_lon)
    grid_max = np.nanmax(grid_lon)

    if grid_min >= 0 and grid_max > 180:
        # ASCAT grid likely uses 0..360
        point_lon = np.where(point_lon < 0, point_lon + 360.0, point_lon)
    elif grid_min < 0 and grid_max <= 180:
        # ASCAT grid likely uses -180..180
        point_lon = np.where(point_lon > 180, point_lon - 360.0, point_lon)

    return point_lon


def read_ascat_grid(grid_file, use_land_flag=True):
    ds = xr.open_dataset(grid_file)

    required_vars = ['lon', 'lat', 'gpi', 'cell']
    for var_name in required_vars:
        if var_name not in ds:
            ds.close()
            raise RuntimeError(f'ASCAT grid file must contain variable: {var_name}')

    lon = ds['lon'].values.astype(float)
    lat = ds['lat'].values.astype(float)
    gpi = ds['gpi'].values.astype(np.int64)
    cell = ds['cell'].values.astype(np.int64)

    valid = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(gpi) & np.isfinite(cell)

    if use_land_flag and 'land_flag' in ds:
        land_flag = ds['land_flag'].values
        valid &= land_flag == 1

    grid = {
        'lon': lon[valid],
        'lat': lat[valid],
        'gpi': gpi[valid],
        'cell': cell[valid]
    }

    ds.close()

    if grid['lon'].size == 0:
        raise RuntimeError('No valid points found in ASCAT grid')

    return grid


def build_ascat_kdtree(ascat_grid):
    points = np.column_stack((ascat_grid['lon'], ascat_grid['lat']))
    tree = cKDTree(points)
    return tree


def match_points_to_ascat_grid(lon_flat, lat_flat, ascat_grid, ascat_tree, max_distance_km):
    lon_query = normalize_longitudes_to_grid(lon_flat, ascat_grid['lon'])
    lat_query = np.asarray(lat_flat, dtype=float)

    query_points = np.column_stack((lon_query, lat_query))
    dist_deg, idx = ascat_tree.query(query_points, k=1)

    # Approximate conversion. Good enough for filtering nearest regular grid candidates.
    max_distance_deg = float(max_distance_km) / 111.0
    valid_match = dist_deg <= max_distance_deg

    return valid_match, idx, lon_query, lat_query
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# cell utils

def update_cell_point(cell_data, gpi, row, col, lon, lat, soil_moisture, point_time):
    point_key = int(gpi)

    record = {
        'gpi': int(gpi),
        'row': int(row),
        'col': int(col),
        'longitude': float(lon),
        'latitude': float(lat),
        'soil_moisture': float(soil_moisture),
        'time': point_time
    }

    if point_key not in cell_data:
        cell_data[point_key] = record
    else:
        if point_time > cell_data[point_key]['time']:
            cell_data[point_key] = record

    return cell_data


def read_cell_netcdf(file_name):
    cell_data = {}

    if not os.path.exists(file_name):
        return cell_data

    ds = xr.open_dataset(file_name)

    gpis = ds['gpi'].values
    rows = ds['row'].values if 'row' in ds else np.full(gpis.shape, -1, dtype=np.int32)
    cols = ds['col'].values if 'col' in ds else np.full(gpis.shape, -1, dtype=np.int32)

    if 'longitude' in ds:
        lons = ds['longitude'].values
    else:
        lons = ds['lon'].values

    if 'latitude' in ds:
        lats = ds['latitude'].values
    else:
        lats = ds['lat'].values

    sms = ds['soil_moisture'].values
    times = pd.to_datetime(ds['time'].values).to_pydatetime()

    ds.close()

    for gpi, row, col, lon, lat, sm, point_time in zip(gpis, rows, cols, lons, lats, sms, times):
        cell_data[int(gpi)] = {
            'gpi': int(gpi),
            'row': int(row),
            'col': int(col),
            'longitude': float(lon),
            'latitude': float(lat),
            'soil_moisture': float(sm),
            'time': point_time
        }

    return cell_data


def write_cell_netcdf(file_name, cell_data, attrs=None):
    make_folder(os.path.dirname(file_name))

    records = list(cell_data.values())
    if len(records) == 0:
        return

    df = pd.DataFrame(records)
    df['time'] = pd.to_datetime(df['time'])

    # Sort by time and keep newest value for common gpi.
    df = df.sort_values('time')
    df = df.drop_duplicates(subset=['gpi'], keep='last')
    df = df.sort_values(['gpi'])

    if attrs is None:
        attrs = {}

    ds = xr.Dataset(
        data_vars={
            'gpi': ('locations', df['gpi'].values.astype('int64')),
            'soil_moisture': ('locations', df['soil_moisture'].values.astype('float32')),

            # ASCAT-compatible names
            'lon': ('locations', df['longitude'].values.astype('float32')),
            'lat': ('locations', df['latitude'].values.astype('float32')),

            # Optional aliases
            'longitude': ('locations', df['longitude'].values.astype('float32')),
            'latitude': ('locations', df['latitude'].values.astype('float32')),

            # Meaningful in grid_domain; set to -1 in grid_ascat.
            'row': ('locations', df['row'].values.astype('int32')),
            'col': ('locations', df['col'].values.astype('int32')),

            'time': ('locations', df['time'].values.astype('datetime64[ns]'))
        },
        coords={
            'locations': np.arange(len(df), dtype='int32')
        },
        attrs=attrs
    )

    ds['gpi'].attrs['long_name'] = 'grid point index'

    ds['soil_moisture'].attrs['long_name'] = 'SMAP soil moisture'
    ds['soil_moisture'].attrs['units'] = 'm3 m-3'

    ds['lon'].attrs['long_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees_east'
    ds['lat'].attrs['long_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees_north'

    ds['longitude'].attrs['long_name'] = 'longitude'
    ds['longitude'].attrs['units'] = 'degrees_east'
    ds['latitude'].attrs['long_name'] = 'latitude'
    ds['latitude'].attrs['units'] = 'degrees_north'

    ds['row'].attrs['long_name'] = 'reference grid row; -1 when grid mode is grid_ascat'
    ds['col'].attrs['long_name'] = 'reference grid column; -1 when grid mode is grid_ascat'

    ds.to_netcdf(file_name)
    ds.close()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# processing utils - grid_domain

def organize_smap_swath2cell_grid_domain(settings, smap_files, date_start, date_end):
    grid_settings = settings.get('grid', {})

    reference_grid = grid_settings['reference_grid']
    n_cell_rows = int(grid_settings.get('n_rows', 3))
    n_cell_cols = int(grid_settings.get('n_cols', 3))

    with rasterio.open(reference_grid) as ref:
        reference_transform = ref.transform
        reference_height = ref.height
        reference_width = ref.width
        reference_data = ref.read(1)

    domain_bounds = compute_domain_bounds(reference_data)

    alg_logger.info(
        f' ::: Domain bounds rows/cols: '
        f'{domain_bounds[0]}-{domain_bounds[1]} / {domain_bounds[2]}-{domain_bounds[3]}'
    )
    alg_logger.info(f' ::: Domain split: {n_cell_rows} x {n_cell_cols}')

    cells_runtime = {}

    for file_path in smap_files:
        file_time = extract_datetime_from_filename(file_path.name, settings)
        if file_time is None:
            alg_logger.warning(f' ===> Time not found in filename: {file_path.name}')
            continue

        alg_logger.info(f' ----> Processing: {file_path.name}')
        alg_logger.info(f' ----> File time : {file_time.strftime("%Y-%m-%d %H:%M:%S")}')

        lon, lat, sm = read_smap_file(file_path, settings)

        lon_flat = lon.ravel()
        lat_flat = lat.ravel()
        sm_flat = sm.ravel()

        valid = np.isfinite(lon_flat) & np.isfinite(lat_flat) & np.isfinite(sm_flat)
        lon_flat = lon_flat[valid]
        lat_flat = lat_flat[valid]
        sm_flat = sm_flat[valid]

        rows, cols = rowcol(reference_transform, lon_flat, lat_flat)
        rows = np.array(rows)
        cols = np.array(cols)

        inside = (
            (rows >= 0) &
            (rows < reference_height) &
            (cols >= 0) &
            (cols < reference_width)
        )

        rows = rows[inside]
        cols = cols[inside]
        lon_flat = lon_flat[inside]
        lat_flat = lat_flat[inside]
        sm_flat = sm_flat[inside]

        mask_values = reference_data[rows, cols]
        valid_mask = np.isfinite(mask_values)

        rows = rows[valid_mask]
        cols = cols[valid_mask]
        lon_flat = lon_flat[valid_mask]
        lat_flat = lat_flat[valid_mask]
        sm_flat = sm_flat[valid_mask]

        alg_logger.info(f' ----> Valid points over domain grid: {len(sm_flat)}')

        for row, col, point_lon, point_lat, point_sm in zip(rows, cols, lon_flat, lat_flat, sm_flat):
            gpi = compute_gpi_from_row_col(row, col, reference_width)

            cell_id = compute_domain_cell_id(
                row=row,
                col=col,
                domain_bounds=domain_bounds,
                n_cell_rows=n_cell_rows,
                n_cell_cols=n_cell_cols
            )

            if cell_id not in cells_runtime:
                cells_runtime[cell_id] = {}

            cells_runtime[cell_id] = update_cell_point(
                cell_data=cells_runtime[cell_id],
                gpi=gpi,
                row=row,
                col=col,
                lon=point_lon,
                lat=point_lat,
                soil_moisture=point_sm,
                point_time=file_time
            )

    written_files = write_runtime_cells(
        settings=settings,
        cells_runtime=cells_runtime,
        date_start=date_start,
        date_end=date_end,
        attrs_extra={
            'grid_mode': 'grid_domain',
            'reference_grid': reference_grid,
            'cell_rows': int(n_cell_rows),
            'cell_cols': int(n_cell_cols),
            'gpi_definition': 'row * grid_width + col'
        }
    )

    return cells_runtime, written_files
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# processing utils - grid_ascat

def organize_smap_swath2cell_grid_ascat(settings, smap_files, date_start, date_end):
    grid_settings = settings.get('grid', {})

    reference_grid = grid_settings['reference_grid']
    max_distance_km = float(grid_settings.get('max_distance_km', 25))
    use_land_flag = bool(grid_settings.get('use_land_flag', True))
    output_smap_coordinates = bool(grid_settings.get('output_smap_coordinates', False))

    alg_logger.info(f' ::: Read ASCAT grid: {reference_grid}')
    ascat_grid = read_ascat_grid(reference_grid, use_land_flag=use_land_flag)
    ascat_tree = build_ascat_kdtree(ascat_grid)

    alg_logger.info(f' ::: ASCAT valid grid points: {ascat_grid["gpi"].size}')
    alg_logger.info(f' ::: ASCAT max matching distance: {max_distance_km} km')

    cells_runtime = {}

    for file_path in smap_files:
        file_time = extract_datetime_from_filename(file_path.name, settings)
        if file_time is None:
            alg_logger.warning(f' ===> Time not found in filename: {file_path.name}')
            continue

        alg_logger.info(f' ----> Processing: {file_path.name}')
        alg_logger.info(f' ----> File time : {file_time.strftime("%Y-%m-%d %H:%M:%S")}')

        lon, lat, sm = read_smap_file(file_path, settings)

        lon_flat = lon.ravel()
        lat_flat = lat.ravel()
        sm_flat = sm.ravel()

        valid = np.isfinite(lon_flat) & np.isfinite(lat_flat) & np.isfinite(sm_flat)
        lon_flat = lon_flat[valid]
        lat_flat = lat_flat[valid]
        sm_flat = sm_flat[valid]

        valid_match, idx, lon_query, lat_query = match_points_to_ascat_grid(
            lon_flat=lon_flat,
            lat_flat=lat_flat,
            ascat_grid=ascat_grid,
            ascat_tree=ascat_tree,
            max_distance_km=max_distance_km
        )

        idx = idx[valid_match]
        sm_flat = sm_flat[valid_match]

        if output_smap_coordinates:
            out_lon = lon_query[valid_match]
            out_lat = lat_query[valid_match]
        else:
            out_lon = ascat_grid['lon'][idx]
            out_lat = ascat_grid['lat'][idx]

        gpis = ascat_grid['gpi'][idx]
        cells = ascat_grid['cell'][idx]

        alg_logger.info(f' ----> Valid points matched to ASCAT grid: {len(sm_flat)}')

        for gpi, cell_id, point_lon, point_lat, point_sm in zip(gpis, cells, out_lon, out_lat, sm_flat):
            cell_id = int(cell_id)

            if cell_id not in cells_runtime:
                cells_runtime[cell_id] = {}

            cells_runtime[cell_id] = update_cell_point(
                cell_data=cells_runtime[cell_id],
                gpi=gpi,
                row=-1,
                col=-1,
                lon=point_lon,
                lat=point_lat,
                soil_moisture=point_sm,
                point_time=file_time
            )

    written_files = write_runtime_cells(
        settings=settings,
        cells_runtime=cells_runtime,
        date_start=date_start,
        date_end=date_end,
        attrs_extra={
            'grid_mode': 'grid_ascat',
            'reference_grid': reference_grid,
            'max_distance_km': float(max_distance_km),
            'use_land_flag': int(use_land_flag),
            'output_smap_coordinates': int(output_smap_coordinates),
            'gpi_definition': 'native ASCAT gpi from reference grid',
            'cell_definition': 'native ASCAT cell from reference grid'
        }
    )

    return cells_runtime, written_files
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# shared writer

def write_runtime_cells(settings, cells_runtime, date_start, date_end, attrs_extra=None):
    if attrs_extra is None:
        attrs_extra = {}

    written_files = []

    for cell_id, new_cell_data in sorted(cells_runtime.items()):
        cell_file = define_cell_file(settings, cell_id, date_end)

        old_cell_data = read_cell_netcdf(cell_file)

        for point_data in new_cell_data.values():
            old_cell_data = update_cell_point(
                cell_data=old_cell_data,
                gpi=point_data['gpi'],
                row=point_data['row'],
                col=point_data['col'],
                lon=point_data['longitude'],
                lat=point_data['latitude'],
                soil_moisture=point_data['soil_moisture'],
                point_time=point_data['time']
            )

        attrs = {
            'source_product': 'SMAP_SPL2SMP_E',
            'processing': 'swath2cell',
            'time_start': date_start.strftime('%Y-%m-%d %H:%M:%S'),
            'time_end': date_end.strftime('%Y-%m-%d %H:%M:%S'),
            'cell_id': int(cell_id),
            'rule_common_points': 'keep latest observation by time for same gpi'
        }
        attrs.update(attrs_extra)

        write_cell_netcdf(
            file_name=cell_file,
            cell_data=old_cell_data,
            attrs=attrs
        )

        written_files.append(str(cell_file))
        alg_logger.info(f' ----> Written cell: {cell_file}')

    return written_files
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm main

def main():
    args = get_args()
    settings = read_file_json(args.settings_file)

    # Validate before processing.
    mode = validate_settings(settings)

    # Time window is needed also for file logger time tags.
    date_start, date_end = get_time_window(args, settings)

    get_logger(settings, time_obj=date_end)

    alg_logger.info(' ============================================================================ ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> START ... ')
    alg_logger.info(' ')

    start_time = time.time()

    alg_logger.info(f' ::: Grid mode : {mode}')
    alg_logger.info(f' ::: Time start: {date_start.strftime("%Y-%m-%d %H:%M:%S")}')
    alg_logger.info(f' ::: Time end  : {date_end.strftime("%Y-%m-%d %H:%M:%S")}')

    alg_logger.info(' ---> Collect SMAP files ... ')

    smap_files = collect_smap_files(
        settings=settings,
        date_start=date_start,
        date_end=date_end
    )

    alg_logger.info(f' ::: SMAP files selected: {len(smap_files)}')

    for file_path in smap_files:
        alg_logger.info(f' ----> Selected: {file_path}')

    if len(smap_files) == 0:
        alg_logger.warning(' ===> No SMAP files found in selected time window')
        alg_logger.info(' ---> Collect SMAP files ... DONE')
        return 0

    alg_logger.info(' ---> Collect SMAP files ... DONE')

    alg_logger.info(' ---> Organize SMAP swath to cell ... ')

    if mode == 'grid_domain':
        cells_runtime, written_files = organize_smap_swath2cell_grid_domain(
            settings=settings,
            smap_files=smap_files,
            date_start=date_start,
            date_end=date_end
        )
    elif mode == 'grid_ascat':
        cells_runtime, written_files = organize_smap_swath2cell_grid_ascat(
            settings=settings,
            smap_files=smap_files,
            date_start=date_start,
            date_end=date_end
        )
    else:
        raise RuntimeError(f'Unsupported grid mode: {mode}')

    alg_logger.info(f' ::: Cells written/updated: {len(cells_runtime)}')
    alg_logger.info(f' ::: Files written/updated: {len(written_files)}')
    alg_logger.info(' ---> Organize SMAP swath to cell ... DONE')

    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(' ')
    alg_logger.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    alg_logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    alg_logger.info(' ==> ... END')
    alg_logger.info(' ==> Bye, Bye')
    alg_logger.info(' ============================================================================ ')

    return 0
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# call script from external library
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(1)
# ----------------------------------------------------------------------------------------------------------------------