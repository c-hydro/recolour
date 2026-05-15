#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - SSM H122 CONVERT - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20260430'
__version__ = '1.4.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_convert_ssm_h122.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20260420 (1.3.0) --> Adapt configuration to parameters/source/destination sections
20260430 (1.4.0) --> Add grid section with grid_ascat_default/grid_ascat_file modes
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import time
import os
import re
import sys
import json
import argparse

from glob import glob
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Application to convert ssm h122 files to gpi/cell outputs'
alg_type = 'Package'
alg_version = '1.4.0'
alg_release = '2026-04-30'

logger = logging.getLogger("app_convert")
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# constants
TIME_FMT_CLI = "%Y-%m-%d %H:%M"
TIME_FMT_FILE = "%Y%m%d%H%M%S"

FILENAME_RE = re.compile(
    r"^W_IT-HSAF-ROME,SAT,SSM-ASCAT-(?P<sat>METOP[BC])-6\.25km-H122_C_LIIB_"
    r"(?P<created>\d{14})_(?P<start>\d{14})_(?P<end>\d{14})____\.nc$"
)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# dataclass
@dataclass(frozen=True)
class FileInfo:
    path: str
    created_time: datetime
    sensing_start: datetime
    sensing_end: datetime
    sat: str
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# basic utils
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    else:
        raise FileNotFoundError(f'File "{file_name}" not found. Exit')


def make_folder(path_name):
    if path_name is not None and path_name != "":
        os.makedirs(path_name, exist_ok=True)


def normalize_time_to_hour(time_obj):
    return time_obj.replace(minute=0, second=0, microsecond=0)


def normalize_time_to_midnight(time_obj):
    return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)


def parse_time(time_string):
    try:
        return datetime.strptime(time_string, TIME_FMT_CLI)
    except ValueError as exc:
        raise ValueError('Time must have format "YYYY-MM-DD HH:MM"') from exc


def parse_reference_time(time_string=None):
    if time_string is None:
        return normalize_time_to_hour(datetime.now())
    return normalize_time_to_hour(parse_time(time_string))


def overlaps(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or a_start > b_end)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# time methods
def parse_time_delta(time_string):
    """
    Supported formats:
    - H, D
    - 1H, 3H, 24H
    - 1D, 2D, 5D
    """
    if time_string is None:
        raise ValueError("Time string cannot be None")

    time_string = str(time_string).strip().upper()

    if time_string == "H":
        time_string = "1H"
    elif time_string == "D":
        time_string = "1D"

    match = re.match(r"^(?P<value>\d+)(?P<unit>[HD])$", time_string)
    if match is None:
        raise ValueError('Time string must have format like "H", "D", "1H", "24H", "2D"')

    time_value = int(match.group("value"))
    time_unit = match.group("unit")

    if time_unit == "H":
        return timedelta(hours=time_value)
    elif time_unit == "D":
        return timedelta(days=time_value)
    else:
        raise ValueError('Unsupported time unit. Use "H" or "D"')


def resolve_time_window(settings, reference_time):
    time_settings = settings.get("time", {})

    time_start_str = time_settings.get("time_start")
    time_end_str = time_settings.get("time_end")
    time_delta_str = time_settings.get("time_delta", time_settings.get("time_period", "2D"))

    floor_start_to_midnight = bool(time_settings.get("floor_start_to_midnight", False))
    floor_end_to_midnight = bool(time_settings.get("floor_end_to_midnight", False))

    if (time_start_str is not None) and (time_end_str is not None):
        time_start = normalize_time_to_hour(parse_time(time_start_str))
        time_end = normalize_time_to_hour(parse_time(time_end_str))

    elif (time_start_str is None) and (time_end_str is None):
        time_end = normalize_time_to_hour(reference_time)
        time_delta = parse_time_delta(time_delta_str)
        time_start = normalize_time_to_hour(time_end - time_delta)

    else:
        raise RuntimeError("time_start and time_end must be both provided, or neither provided")

    if floor_start_to_midnight:
        time_start = normalize_time_to_midnight(time_start)

    if floor_end_to_midnight:
        time_end = normalize_time_to_midnight(time_end)

    if time_end < time_start:
        raise RuntimeError(f"time_end {time_end} is earlier than time_start {time_start}")

    return time_start, time_end


def iter_time_steps(settings, time_start, time_end):
    time_settings = settings.get("time", {})
    time_frequency_str = time_settings.get("time_frequency", "H")
    time_frequency_delta = parse_time_delta(time_frequency_str)

    current_time = normalize_time_to_hour(time_start)
    final_time = normalize_time_to_hour(time_end)

    while current_time <= final_time:
        yield current_time
        current_time += time_frequency_delta
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# template methods
def resolve_time_tags(path_template, time_dict):
    path_resolved = str(path_template)

    for time_key, time_value in time_dict.items():
        if time_value is None:
            continue

        tag_pattern = r"\{" + re.escape(time_key) + r":([^}]+)\}"

        def replace_tag(match):
            time_format = match.group(1)
            return time_value.strftime(time_format)

        path_resolved = re.sub(tag_pattern, replace_tag, path_resolved)

    return path_resolved


def resolve_destination_filename(filename_template, tags):
    filename_resolved = resolve_time_tags(filename_template, tags)

    for tag_key, tag_value in tags.items():
        if isinstance(tag_value, (str, int, float)):
            filename_resolved = filename_resolved.replace("{" + str(tag_key) + "}", str(tag_value))

    return filename_resolved
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# file methods
def parse_file_info(file_path):
    file_name = os.path.basename(file_path)
    match = FILENAME_RE.match(file_name)
    if match is None:
        return None

    return FileInfo(
        path=file_path,
        created_time=datetime.strptime(match.group("created"), TIME_FMT_FILE),
        sensing_start=datetime.strptime(match.group("start"), TIME_FMT_FILE),
        sensing_end=datetime.strptime(match.group("end"), TIME_FMT_FILE),
        sat=match.group("sat"),
    )


def discover_files(settings, time_start, time_end, reference_time):
    source_settings = settings.get("source", {})
    source_folder_template = source_settings["folder"]
    source_filename_pattern = source_settings.get("filename", "*.nc")

    selected_files = []
    seen_files = set()

    for time_step in iter_time_steps(settings, time_start, time_end):
        src_folder = resolve_time_tags(
            source_folder_template,
            {
                "time_step": time_step,
                "time_start": time_start,
                "time_end": time_end,
                "time_run": reference_time
            }
        )

        logger.info(f" ----> Search source folder: {src_folder}")

        if not os.path.exists(src_folder):
            continue

        file_list = sorted(glob(os.path.join(src_folder, source_filename_pattern)))
        for file_path in file_list:
            if file_path in seen_files:
                continue

            file_info = parse_file_info(file_path)
            if file_info is None:
                continue

            seen_files.add(file_path)

            if overlaps(file_info.sensing_start, file_info.sensing_end, time_start, time_end):
                selected_files.append(file_info)

    return sorted(
        selected_files,
        key=lambda file_obj: (file_obj.sensing_start, file_obj.sensing_end, os.path.basename(file_obj.path))
    )
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# dataframe methods
def collect_to_dataframe(settings, file_list):
    source_settings = settings.get("source", {})

    value_col = source_settings.get("value_col", "surface_soil_moisture")
    agg_method = source_settings.get("agg", "first")
    include_all_vars = bool(source_settings.get("include_all_vars", False))
    dropna_value = bool(source_settings.get("dropna_value", False))

    frames = []

    for file_info in file_list:
        logger.info(
            f" ----> File ref: {file_info.created_time} | "
            f"{file_info.sensing_start} -> {file_info.sensing_end} | {file_info.path}"
        )

        with xr.open_dataset(file_info.path) as dataset:
            if "location_id" not in dataset:
                raise KeyError(f'location_id missing in "{file_info.path}"')
            if value_col not in dataset:
                raise KeyError(f'{value_col} missing in "{file_info.path}"')
            if "latitude" not in dataset or "longitude" not in dataset or "time" not in dataset:
                raise KeyError(f'latitude/longitude/time missing in "{file_info.path}"')

            cols = {
                "gpi": dataset["location_id"].values.astype(np.int64),
                "lat": dataset["latitude"].values.astype(np.float64),
                "lon": dataset["longitude"].values.astype(np.float64),
                "time": pd.to_datetime(dataset["time"].values),
                value_col: dataset[value_col].values,
                "created_time": np.full(
                    dataset["location_id"].shape[0], file_info.created_time, dtype="datetime64[ns]"
                ),
                "source_file": np.full(
                    dataset["location_id"].shape[0], os.path.basename(file_info.path), dtype=object
                ),
                "satellite": np.full(
                    dataset["location_id"].shape[0], file_info.sat, dtype=object
                ),
            }

            if include_all_vars:
                for var_name, data_array in dataset.data_vars.items():
                    if var_name in cols or var_name == "location_id":
                        continue
                    if data_array.ndim == 1 and data_array.dims == ("obs",):
                        cols[var_name] = data_array.values

            frames.append(pd.DataFrame(cols))

    if not frames:
        return pd.DataFrame(columns=["gpi", "lat", "lon", "time", value_col])

    df = pd.concat(frames, ignore_index=True)

    if dropna_value:
        df = df[df[value_col].notna()].copy()

    agg_map = {
        "lat": "first",
        "lon": "first",
        "time": agg_method,
        value_col: agg_method,
        "source_file": "first",
        "satellite": "first",
    }

    for col_name in df.columns:
        if col_name != "gpi" and col_name not in agg_map:
            if pd.api.types.is_numeric_dtype(df[col_name]) or np.issubdtype(df[col_name].dtype, np.datetime64):
                agg_map[col_name] = agg_method
            else:
                agg_map[col_name] = "first"

    df_out = df.groupby("gpi", sort=True, dropna=False).agg(agg_map).reset_index()
    return df_out
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# grid/domain methods
def validate_grid_settings(grid_settings):
    grid_mode = grid_settings.get("mode", "grid_ascat_default").lower()

    if grid_mode not in ["grid_ascat_default", "grid_ascat_file"]:
        raise RuntimeError(
            'Unsupported grid.mode. Use "grid_ascat_default" or "grid_ascat_file"'
        )

    if grid_mode == "grid_ascat_file" and grid_settings.get("reference_grid") is None:
        raise RuntimeError('grid.reference_grid must be defined for mode "grid_ascat_file"')

    return grid_mode

def get_grid_metadata(grid_settings):
    grid_mode = grid_settings.get("mode", "grid_ascat_default").lower()
    grid_name = grid_settings.get("name", "fibgrid_6.25")
    reference_grid = grid_settings.get("reference_grid", None)
    id_digits = int(grid_settings.get("id_digits", 4))

    return grid_mode, grid_name, reference_grid, id_digits

def add_tuw_cell_from_ascat_package(df, grid_name="fibgrid_6.25"):
    """
    Add TUW/FibGrid cell id using the tuw-geo/ascat package.

    For H122 use:
        grid_name = "fibgrid_6.25"

    For H29 use:
        grid_name = "fibgrid_12.5"
    """

    try:
        from ascat.grids.grid_registry import GridRegistry
    except Exception as exc:
        raise RuntimeError(
            'Missing package "ascat". Install it first, for example with:\n'
            'pip install ascat\n'
            'or use the conda environment from https://github.com/TUW-GEO/ascat'
        ) from exc

    registry = GridRegistry()
    grid = registry.get(grid_name)

    gpis = df["gpi"].to_numpy(np.int64)

    # Preferred accessor
    if hasattr(grid, "gpi2cell"):
        cells = np.array(
            [grid.gpi2cell(int(gpi)) for gpi in gpis],
            dtype=np.int32
        )

    # Fallback for pygeogrids-style grid objects
    elif hasattr(grid, "activegpis") and hasattr(grid, "activearrcell"):
        lookup = dict(
            zip(
                np.asarray(grid.activegpis, dtype=np.int64),
                np.asarray(grid.activearrcell, dtype=np.int32)
            )
        )

        missing = [int(gpi) for gpi in gpis if int(gpi) not in lookup]
        if missing:
            raise KeyError(
                f"{len(missing)} gpis not found in ASCAT grid {grid_name}. "
                f"First missing gpis: {missing[:20]}"
            )

        cells = np.array([lookup[int(gpi)] for gpi in gpis], dtype=np.int32)

    else:
        raise AttributeError(
            "Unsupported ASCAT grid object: no gpi2cell, activegpis, or activearrcell found"
        )

    df_out = df.copy()
    df_out.insert(1, "cell", cells)

    return df_out


def add_tuw_cell_from_reference_file(df, reference_grid, missing_policy="drop"):
    gpis = df["gpi"].to_numpy(np.int64)

    with xr.open_dataset(reference_grid) as grid_ds:
        if "gpi" not in grid_ds:
            raise KeyError(f'gpi missing in reference grid "{reference_grid}"')
        if "cell" not in grid_ds:
            raise KeyError(f'cell missing in reference grid "{reference_grid}"')

        grid_gpi = grid_ds["gpi"].values.astype(np.int64)
        grid_cell = grid_ds["cell"].values.astype(np.int32)

    lookup = dict(zip(grid_gpi, grid_cell))

    cells = np.array(
        [lookup.get(int(gpi), -1) for gpi in gpis],
        dtype=np.int32
    )

    missing_mask = cells == -1
    n_missing = int(missing_mask.sum())

    if n_missing > 0:
        missing_gpis = np.unique(gpis[missing_mask])
        logger.warning(
            f" ----> Missing {n_missing} rows / {len(missing_gpis)} unique gpis "
            f"in reference grid {reference_grid}. First missing gpis: "
            f"{missing_gpis[:20].tolist()}"
        )

        if missing_policy == "raise":
            raise KeyError(
                f"{len(missing_gpis)} gpis not found in reference grid. "
                f"First missing gpis: {missing_gpis[:20].tolist()}"
            )

        if missing_policy == "drop":
            df = df.loc[~missing_mask].copy()
            cells = cells[~missing_mask]

        elif missing_policy == "keep":
            df = df.copy()

        else:
            raise RuntimeError('grid.missing_policy must be "raise", "drop", or "keep"')

    df_out = df.copy()
    df_out.insert(1, "cell", cells)

    return df_out

def add_tuw_cell_from_default_grid(df, grid_name="fibgrid_6.25"):
    try:
        from ascat.grids.grid_registry import GridRegistry
    except Exception as exc:
        raise RuntimeError(
            'mode "grid_ascat_default" requires the ascat package. '
            'Use mode "grid_ascat_file" with grid.reference_grid otherwise.'
        ) from exc

    registry = GridRegistry()
    grid = registry.get(grid_name)

    gpis = df["gpi"].to_numpy(np.int64)

    if hasattr(grid, "gpi2cell"):
        cells = np.array([grid.gpi2cell(int(gpi)) for gpi in gpis], dtype=np.int32)

    elif hasattr(grid, "activegpis") and hasattr(grid, "activearrcell"):
        lookup = dict(
            zip(
                np.asarray(grid.activegpis, dtype=np.int64),
                np.asarray(grid.activearrcell, dtype=np.int32)
            )
        )
        try:
            cells = np.array([lookup[int(gpi)] for gpi in gpis], dtype=np.int32)
        except KeyError as exc:
            raise KeyError(f"gpi {exc.args[0]} not found in {grid_name}") from exc

    else:
        raise AttributeError("Unsupported ASCAT grid object: no known gpi->cell accessor found")

    df_out = df.copy()
    df_out.insert(1, "cell", cells)

    return df_out


def add_tuw_cell(df, grid_settings):
    grid_mode, grid_name, reference_grid, _ = get_grid_metadata(grid_settings)

    if grid_mode == "grid_ascat_default":
        logger.info(f" ----> Select grid mode: {grid_mode}")
        logger.info(f" ----> Select ASCAT grid: {grid_name}")
        return add_tuw_cell_from_ascat_package(
            df=df,
            grid_name=grid_name
        )

    if grid_mode == "grid_ascat_file":
        logger.info(f" ----> Select grid mode: {grid_mode}")
        logger.info(f" ----> Select grid file: {reference_grid}")
        return add_tuw_cell_from_reference_file(
            df=df,
            reference_grid=reference_grid
        )

    raise RuntimeError(
        'Unsupported grid.mode. Use "grid_ascat_default" or "grid_ascat_file"'
    )

def filter_by_domain_cells(df, allowed_cells):
    return df[df["cell"].isin(allowed_cells)].copy()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# output methods
def write_gpi_table(df, out_folder, filename_template, time_start, time_end, reference_time, overwrite=True):
    make_folder(out_folder)

    out_filename = resolve_destination_filename(
        filename_template,
        {
            "time_step": reference_time,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": reference_time
        }
    )

    out_path = os.path.join(out_folder, out_filename)

    if os.path.exists(out_path) and not overwrite:
        logger.info(f" ----> Skip existing gpi file: {out_path}")
        return out_path

    file_ext = os.path.splitext(out_path)[1].lower()

    if file_ext == ".parquet":
        df.to_parquet(out_path, index=False)
    elif file_ext == ".csv":
        df.to_csv(out_path, index=False)
    else:
        try:
            import pyarrow  # noqa: F401
            out_path = os.path.splitext(out_path)[0] + ".parquet"
            df.to_parquet(out_path, index=False)
        except Exception:
            out_path = os.path.splitext(out_path)[0] + ".csv"
            df.to_csv(out_path, index=False)

    return out_path


def build_dataset_from_part(part_df, attrs):
    data_vars = {}
    for col_name in part_df.columns:
        if col_name == "cell":
            continue

        values = part_df[col_name].to_numpy()
        if np.issubdtype(values.dtype, np.datetime64):
            data_vars[col_name] = ("locations", values.astype("datetime64[ns]"))
        else:
            data_vars[col_name] = ("locations", values)

    dataset = xr.Dataset(
        data_vars=data_vars,
        coords={"locations": np.arange(len(part_df), dtype=np.int32)},
        attrs=attrs,
    )
    return dataset


def write_cell_files(
        df, out_folder, filename_template, time_start, time_end, reference_time,
        grid_name="fibgrid_6.25", grid_mode="grid_ascat_default",
        reference_grid=None, id_digits=4, overwrite=True):

    make_folder(out_folder)

    written_files = []
    attrs = {
        "source_product": "H122",
        "grid_mode": grid_mode,
        "grid_mapping_name": grid_name,
        "reference_grid": "" if reference_grid is None else reference_grid,
        "time_start": time_start.strftime("%Y-%m-%d %H:%M:%S"),
        "time_end": time_end.strftime("%Y-%m-%d %H:%M:%S"),
    }

    for cell_id, part_df in df.groupby("cell", sort=True):
        out_filename = resolve_destination_filename(
            filename_template,
            {
                "cell_n": f"{int(cell_id):0{id_digits}d}",
                "cell": int(cell_id),
                "time_step": reference_time,
                "time_start": time_start,
                "time_end": time_end,
                "time_run": reference_time
            }
        )
        out_path = os.path.join(out_folder, out_filename)

        if os.path.exists(out_path) and not overwrite:
            logger.info(f" ----> Skip existing cell file: {out_path}")
            continue

        dataset = build_dataset_from_part(part_df, attrs)
        dataset.to_netcdf(out_path)
        written_files.append(out_path)

    return written_files
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# process methods
def process_data(settings, reference_time):
    parameters_settings = settings.get("parameters", {})
    grid_settings = settings.get("grid", {})
    destination_settings = settings.get("destination", {})

    output_mode = destination_settings.get("mode", "cell").lower()
    output_dry_run = bool(destination_settings.get("dry_run", False))
    output_overwrite = bool(destination_settings.get("overwrite", True))

    allowed_cells = set(parameters_settings.get("cells", []))

    grid_mode = validate_grid_settings(grid_settings)
    grid_mode, grid_name, reference_grid, id_digits = get_grid_metadata(grid_settings)

    stats = {
        "selected_files": 0,
        "written_files": 0,
        "rows_before_filter": 0,
        "rows_after_filter": 0,
        "errors": 0,
    }

    written_objects = []

    time_start, time_end = resolve_time_window(settings, reference_time)

    destination_folder = resolve_time_tags(
        destination_settings["folder"],
        {
            "time_step": reference_time,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": reference_time
        }
    )

    logger.info(f" ----> Time start: {time_start}")
    logger.info(f" ----> Time end:   {time_end}")
    logger.info(f" ----> Destination folder: {destination_folder}")
    logger.info(f" ----> Grid mode: {grid_mode}")
    logger.info(f" ----> Grid name: {grid_name}")
    logger.info(f" ----> Grid reference file: {reference_grid}")
    logger.info(f" ----> Cell id digits: {id_digits}")

    file_list = discover_files(settings, time_start, time_end, reference_time)
    stats["selected_files"] = len(file_list)

    logger.info(f" ----> Selected files: {len(file_list)}")
    for file_info in file_list:
        logger.info(
            f" ------> {file_info.sensing_start} -> {file_info.sensing_end} | {file_info.path}"
        )

    if output_dry_run:
        logger.info(" ----> Dry run active. Exit without writing output")
        return stats, written_objects, time_start, time_end, destination_folder

    if not file_list:
        logger.warning(" ----> No files selected. Nothing to do")
        return stats, written_objects, time_start, time_end, destination_folder

    df = collect_to_dataframe(settings, file_list)
    logger.info(f" ----> Unique gpis before domain filter: {len(df)}")

    df = add_tuw_cell(df=df, grid_settings=grid_settings)
    stats["rows_before_filter"] = len(df)

    if allowed_cells:
        df = filter_by_domain_cells(df, allowed_cells)

    stats["rows_after_filter"] = len(df)

    logger.info(
        f" ----> Rows after domain filter: {stats['rows_after_filter']} / {stats['rows_before_filter']}"
    )
    logger.info(f" ----> Allowed cells: {sorted(list(allowed_cells))}")

    if df.empty:
        logger.warning(" ----> No data left after domain filtering. Nothing to write")
        return stats, written_objects, time_start, time_end, destination_folder

    destination_filename = destination_settings.get("filename", "output.nc")

    if output_mode == "gpi":
        out_path = write_gpi_table(
            df=df,
            out_folder=destination_folder,
            filename_template=destination_filename,
            time_start=time_start,
            time_end=time_end,
            reference_time=reference_time,
            overwrite=output_overwrite
        )
        written_objects.append(out_path)
        stats["written_files"] = 1
        logger.info(f" ----> Wrote gpi table: {out_path}")

    elif output_mode == "cell":
        written_objects = write_cell_files(
            df=df,
            out_folder=destination_folder,
            filename_template=destination_filename,
            time_start=time_start,
            time_end=time_end,
            reference_time=reference_time,
            grid_name=grid_name,
            grid_mode=grid_mode,
            reference_grid=reference_grid,
            id_digits=id_digits,
            overwrite=output_overwrite
        )
        stats["written_files"] = len(written_objects)
        logger.info(f" ----> Wrote cell files: {len(written_objects)} in {destination_folder}")

    else:
        raise ValueError('destination.mode must be "gpi" or "cell"')

    return stats, written_objects, time_start, time_end, destination_folder
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# report utils
# helper to collect info report
def collect_report(
        settings, reference_time, time_start, time_end, destination_folder,
        stats, start_time, written_objects=None):

    if written_objects is None:
        written_objects = []

    elapsed = round(time.time() - start_time, 1)

    report = {
        "algorithm": {
            "name": alg_name,
            "version": alg_version,
            "release": alg_release
        },
        "run": {
            "reference_time": reference_time.strftime("%Y-%m-%d %H:%M:%S"),
            "time_start": time_start.strftime("%Y-%m-%d %H:%M:%S"),
            "time_end": time_end.strftime("%Y-%m-%d %H:%M:%S"),
            "elapsed_seconds": elapsed
        },
        "parameters": settings.get("parameters", {}),
        "grid": settings.get("grid", {}),
        "source": settings.get("source", {}),
        "destination": {
            **settings.get("destination", {}),
            "resolved_folder": destination_folder
        },
        "time": settings.get("time", {}),
        "stats": stats,
        "written_objects": written_objects
    }

    return report


# helper to write info report
def write_report(report, file_path):
    folder_name = os.path.dirname(file_path)
    if folder_name:
        make_folder(folder_name)

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=4)


# helper to save info report
def save_report(settings, report, reference_time):
    report_settings = settings.get("report", {})
    if not report_settings.get("enabled", False):
        return

    report_folder_template = report_settings.get("folder", "./report")
    report_filename_template = report_settings.get(
        "filename", "swath2cell_h122_report_%Y%m%d_%H%M.json"
    )

    report_folder = resolve_time_tags(
        report_folder_template,
        {
            "time_step": reference_time,
            "time_start": reference_time,
            "time_end": reference_time,
            "time_run": reference_time
        }
    )

    report_filename = reference_time.strftime(report_filename_template)
    report_path = os.path.join(report_folder, report_filename)

    write_report(report, report_path)
    logger.info(" ----> Report saved: " + report_path)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# execution utils
# helper to configure logger
def get_logger(settings, reference_time=None):
    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder_template = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "swath2cell_h122.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    if reference_time is None:
        reference_time = datetime.now()

    log_folder = resolve_time_tags(
        log_folder_template,
        {
            "time_step": reference_time,
            "time_start": reference_time,
            "time_end": reference_time,
            "time_run": reference_time
        }
    )

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = reference_time.strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    logger.handlers = []
    logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    )

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)


# helper to get arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="Convert H122 swath files to gpi/cell outputs"
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file",
    )

    parser.add_argument(
        "-time",
        dest="time_run",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM". Minutes are rounded down to 00.',
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get arguments
    args = get_args()
    # get settings
    settings = read_file_json(args.settings_file)

    # get reference time
    try:
        reference_time = parse_reference_time(args.time_run)
    except Exception as exc:
        print(f" ===> ERROR: parsing time: {exc}")
        sys.exit(1)

    # get logger
    get_logger(settings, reference_time=reference_time)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    logger.info(" ============================================================================ ")
    logger.info(" ==> " + alg_name + " (Version: " + alg_version + " Release_Date: " + alg_release + ")")
    logger.info(" ==> START ... ")
    logger.info(" ")

    logger.info(f" ---> Settings file:      {args.settings_file}")
    logger.info(f" ---> Source folder:      {settings.get('source', {}).get('folder')}")
    logger.info(f" ---> Source filename:    {settings.get('source', {}).get('filename')}")
    logger.info(f" ---> Destination folder: {settings.get('destination', {}).get('folder')}")
    logger.info(f" ---> Destination file:   {settings.get('destination', {}).get('filename')}")
    logger.info(f" ---> Grid mode:          {settings.get('grid', {}).get('mode')}")
    logger.info(f" ---> Grid name:          {settings.get('grid', {}).get('name')}")
    logger.info(f" ---> Grid reference:     {settings.get('grid', {}).get('reference_grid')}")

    start_time = time.time()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # set time
    logger.info(" ---> Set time ... ")
    try:
        logger.info(f" ----> Get reference time: {reference_time.strftime('%Y-%m-%d %H:%M:%S')}")
        if args.time_run:
            logger.info(f" ----> Set by user: {args.time_run} (rounded to hour)")
        else:
            logger.info(" ----> Set by system (rounded to hour)")
        logger.info(" ---> Set time ... DONE")

    except Exception as exc:
        logger.error(f" ===> ERROR: parsing time: {exc}")
        logger.info(" ---> Set time ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # execute conversion
    logger.info(" ---> Execute conversion ... ")
    try:
        stats, written_objects, time_start, time_end, destination_folder = process_data(
            settings, reference_time
        )
        logger.info(" ---> Execute conversion ... DONE")

    except Exception as exc:
        logger.error(f" ===> ERROR in executing conversion algorithm: {exc}")
        logger.info(" ---> Execute conversion ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # save report
    logger.info(" ---> Save report ...")
    try:
        report = collect_report(
            settings=settings,
            reference_time=reference_time,
            time_start=time_start,
            time_end=time_end,
            destination_folder=destination_folder,
            stats=stats,
            start_time=start_time,
            written_objects=written_objects
        )

        save_report(settings, report, reference_time=reference_time)
        logger.info(" ---> Save report ... DONE")

    except Exception as exc:
        logger.error(f" ===> ERROR saving report: {exc}")
        logger.info(" ---> Save report ... FAILED")
        sys.exit(1)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    alg_time_elapsed = round(time.time() - start_time, 1)

    logger.info(" ")
    logger.info(" ==> " + alg_name + " (Version: " + alg_version + " Release_Date: " + alg_release + ")")
    logger.info(" ==> TIME ELAPSED: " + str(alg_time_elapsed) + " seconds")
    logger.info(" ==> ... END")
    logger.info(" ==> Bye, Bye")
    logger.info(" ============================================================================ ")
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------