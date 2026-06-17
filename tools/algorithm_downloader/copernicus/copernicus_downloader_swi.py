#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - COPERNICUS CLMS SOIL WATER INDEX DOWNLOADER

Downloads selected variables from the Copernicus CLMS Daily Soil Water Index
Global 12.5 km V4 COG catalogue, clips them to a bounding box, and optionally
converts encoded UINT8 SWI/QFLAG values to physical percentages.

Example:
python copernicus_downloader_swi_filtered.py \
    -settings_file copernicus_downloader_swi_filtered.json \
    -time "2026-06-11 00:00"

Required S3 credentials:
export CDSE_S3_ACCESS_KEY="your_access_key"
export CDSE_S3_SECRET_KEY="your_secret_key"
"""

import os
import re
import sys
import json
import time
import atexit
import signal
import logging
import argparse
from pathlib import Path
from datetime import datetime, timedelta

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds


project_name = "recolour"
alg_name = "Application for downloading Copernicus CLMS Soil Water Index files"
alg_type = "Package"
alg_version = "1.1.0"
alg_release = "2026-06-15"

alg_logger = logging.getLogger("copernicus_swi_downloader")
acquired_lock = None

S3_ENDPOINT_URL = "https://eodata.dataspace.copernicus.eu"


# ----------------------------------------------------------------------------------------------------------------------
# basic utils
def read_file_json(file_name):
    if os.path.exists(file_name):
        with open(file_name, "r", encoding="utf-8") as file_handle:
            return json.load(file_handle)
    raise FileNotFoundError(f' ===> File "{file_name}" not found. Exit')


def get_args():
    parser = argparse.ArgumentParser(
        description="Download Copernicus CLMS Soil Water Index files"
    )

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        required=True,
        help="Path to JSON settings file",
    )

    parser.add_argument(
        "-f",
        "-force",
        dest="force_lock",
        action="store_true",
        help="Force removal of stale lock files",
    )

    parser.add_argument(
        "-time",
        dest="time_run",
        default=None,
        help='Reference time in format "YYYY-MM-DD HH:MM"',
    )

    return parser.parse_args()


def get_logger(settings):
    log_settings = settings.get("logging", {})

    level_str = log_settings.get("level", "INFO").upper()
    level = getattr(logging, level_str, logging.INFO)

    log_folder = log_settings.get("folder", "./log")
    log_filename = log_settings.get("filename", "copernicus_swi_downloader.log")
    rotate_daily = log_settings.get("rotate_daily", False)

    os.makedirs(log_folder, exist_ok=True)

    if rotate_daily:
        timestamp = datetime.now().strftime("%Y%m%d")
        log_filename = f"{timestamp}_{log_filename}"

    log_path = os.path.join(log_folder, log_filename)

    alg_logger.handlers = []
    alg_logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(log_path)
    fh.setLevel(level)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)

    alg_logger.addHandler(fh)
    alg_logger.addHandler(ch)


def make_folder(folder_name):
    if folder_name is not None and folder_name != "":
        os.makedirs(folder_name, exist_ok=True)


# ----------------------------------------------------------------------------------------------------------------------
# time utils
def parse_time_string(time_string, time_format="%Y-%m-%d %H:%M"):
    if time_string is None:
        return None

    try:
        return datetime.strptime(time_string, time_format)
    except ValueError as exc:
        raise ValueError(f'time "{time_string}" must have format "{time_format}"') from exc


def parse_date_string(date_string):
    if date_string is None:
        return None

    try:
        return datetime.strptime(date_string, "%Y-%m-%d")
    except ValueError as exc:
        raise ValueError('date must have format "YYYY-MM-DD"') from exc


def round_time_to_day(time_obj):
    return time_obj.replace(hour=0, minute=0, second=0, microsecond=0)


def get_reference_time(args, settings):
    downloader_settings = settings.get("downloader", {})
    round_to_day = downloader_settings.get("reference_time_round_to_day", True)

    if args.time_run is not None:
        reference_time = parse_time_string(args.time_run)
    else:
        reference_time = datetime.now()

    if round_to_day:
        reference_time = round_time_to_day(reference_time)

    return reference_time


def get_time_window(settings):
    time_settings = settings.get("time", {})

    time_start_raw = time_settings.get("time_start", None)
    time_end_raw = time_settings.get("time_end", None)

    time_start = parse_date_string(time_start_raw) if time_start_raw is not None else None
    time_end = parse_date_string(time_end_raw) if time_end_raw is not None else None

    return time_start, time_end


def daterange(start_date, end_date):
    current_date = start_date

    while current_date <= end_date:
        yield current_date
        current_date += timedelta(days=1)


def fill_time_template(path_template, time_obj, variable=None):
    if path_template is None:
        return None

    out = path_template

    if "{time_step" in out:
        pattern = r"\{time_step:\s*(.*?)\}"

        def replace(match):
            time_format = match.group(1)
            return time_obj.strftime(time_format)

        out = re.sub(pattern, replace, out)

    if variable is not None:
        out = out.replace("{variable}", variable)

    return out


# ----------------------------------------------------------------------------------------------------------------------
# run utils
def select_run_mode(args, settings):
    time_start, time_end = get_time_window(settings)
    reference_time = get_reference_time(args, settings)

    downloader_settings = settings.get("downloader", {})
    n_days = int(downloader_settings.get("n_days", 3))

    if (time_start is not None) and (time_end is not None):
        mode = "date_range"
    else:
        mode = "rolling_window"
        time_end = reference_time
        time_start = reference_time - timedelta(days=n_days)

    return mode, reference_time, time_start, time_end, n_days


def acquire_lock(settings, force_lock=False):
    global acquired_lock

    lock_settings = settings.get("lock", {})
    lock_dir = lock_settings["folder"]
    lock_basename = lock_settings["basename"]
    max_slots = int(lock_settings.get("max_slots", 1))

    make_folder(lock_dir)

    for lock_id in range(1, max_slots + 1):
        lock_file = os.path.join(lock_dir, f"{lock_basename}{lock_id}")

        if os.path.exists(lock_file) and force_lock:
            alg_logger.warning(f" ===> Forcing removal of stale lock: {lock_file}")
            try:
                os.remove(lock_file)
            except OSError as exc:
                alg_logger.error(f' ===> Error removing stale lock "{lock_file}": {exc}')
                continue

        if os.path.exists(lock_file):
            continue

        try:
            file_descriptor = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            with os.fdopen(file_descriptor, "w") as file_handle:
                file_handle.write(str(os.getpid()) + "\n")

            acquired_lock = lock_file
            alg_logger.info(f" ----> Acquired slot {lock_id} using {lock_file}")
            return True

        except FileExistsError:
            continue
        except OSError as exc:
            alg_logger.error(f' ===> Error creating lock "{lock_file}": {exc}')

    alg_logger.error(f" ===> All {max_slots} slots are currently in use. Exiting.")
    return False


def release_lock():
    global acquired_lock

    if acquired_lock is not None:
        try:
            if os.path.exists(acquired_lock):
                os.remove(acquired_lock)
                alg_logger.info(f" ----> Released slot ({acquired_lock})")
        except OSError as exc:
            alg_logger.error(f' ===> Error releasing lock "{acquired_lock}": {exc}')
        finally:
            acquired_lock = None


def handle_signal(signum, frame):
    alg_logger.warning(f" ===> Caught signal {signum}. Releasing lock and exiting.")
    release_lock()
    sys.exit(1)


# ----------------------------------------------------------------------------------------------------------------------
# s3 utils
def get_s3_client():
    access_key = os.environ.get("CDSE_S3_ACCESS_KEY")
    secret_key = os.environ.get("CDSE_S3_SECRET_KEY")

    if access_key is None or secret_key is None:
        raise RuntimeError(
            "Missing CDSE_S3_ACCESS_KEY/CDSE_S3_SECRET_KEY environment variables"
        )

    return boto3.client(
        "s3",
        endpoint_url=S3_ENDPOINT_URL,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        config=Config(signature_version="s3v4"),
    )


def parse_s3_path(s3_path):
    s3_path = str(s3_path)

    if not s3_path.startswith("s3://"):
        raise RuntimeError(f"Unsupported S3 path: {s3_path}")

    path = s3_path.replace("s3://", "", 1)
    bucket, key = path.split("/", 1)

    return bucket, key


# ----------------------------------------------------------------------------------------------------------------------
# copernicus swi utils
def extract_variable_from_path(path):
    """Extract variables like SWI040, QFLAG040, SSF from catalogue paths or filenames."""
    text = Path(str(path)).name.upper()
    patterns = [
        r"(?:^|[_-])(SWI\d{3})(?:[_\.-]|$)",
        r"(?:^|[_-])(QFLAG\d{3})(?:[_\.-]|$)",
        r"(?:^|[_-])(SSF)(?:[_\.-]|$)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return "UNKNOWN"


def read_catalogue(settings):
    product_settings = settings.get("product", {})
    catalog_csv = product_settings["catalog_csv"]

    alg_logger.info(" ::: Reading Copernicus SWI catalogue")
    alg_logger.info(f" ::: Catalogue CSV: {catalog_csv}")

    df = pd.read_csv(catalog_csv, sep=";")

    alg_logger.info(f" ::: Catalogue columns: {list(df.columns)}")
    alg_logger.info(f" ::: Catalogue rows: {len(df)}")

    url_column = "s3_path"
    date_column = "nominal_date"

    if url_column not in df.columns:
        raise RuntimeError(f'Column "{url_column}" not found in catalogue')

    if date_column not in df.columns:
        raise RuntimeError(f'Column "{date_column}" not found in catalogue')

    df["time_step"] = pd.to_datetime(df[date_column], errors="coerce").dt.floor("D")
    df = df.dropna(subset=["time_step"])
    df["variable"] = df[url_column].apply(extract_variable_from_path)

    variables_found = sorted(df["variable"].dropna().unique().tolist())
    alg_logger.info(f" ::: Catalogue variables found: {variables_found}")
    alg_logger.info(f" ::: Catalogue valid dates: {len(df)}")

    if df.empty:
        raise RuntimeError("No valid dates found in catalogue")

    alg_logger.info(
        " ::: Catalogue time coverage: "
        f"{df['time_step'].min().strftime('%Y-%m-%d')} -> "
        f"{df['time_step'].max().strftime('%Y-%m-%d')}"
    )

    return df, url_column


def select_catalogue_rows(df, url_column, day, settings):
    """Return one download row per requested SWI variable for the selected day.

    Important: the COG catalogue has one row per day and its s3_path often
    points to a product folder/prefix, not to a specific SWI/QFLAG file.
    Therefore variable filtering cannot be done on the catalogue row itself;
    the variable is resolved later inside the S3 prefix.
    """
    product_settings = settings.get("product", {})
    variables = product_settings.get("variables", None)

    if variables is None or len(variables) == 0:
        raise RuntimeError(
            'No product.variables configured. Example: ["SWI001", "SWI005", "SWI010", "SWI015", "SWI020", "SWI040"]'
        )

    variables_upper = [str(v).upper() for v in variables]

    day_start = pd.Timestamp(day.strftime("%Y-%m-%d"))
    day_end = day_start + pd.Timedelta(days=1)

    selected = df[
        (df["time_step"] >= day_start) &
        (df["time_step"] < day_end)
    ].copy()

    if selected.empty:
        return []

    rows = []

    for _, row in selected.iterrows():
        for variable in variables_upper:
            rows.append(
                {
                    "time_step": row["time_step"].to_pydatetime(),
                    "variable": variable,
                    "url": str(row[url_column]),
                    "source_path": str(row[url_column]),
                }
            )

    return rows

def define_data_folder(settings, time_obj, variable=None):
    downloader_settings = settings.get("downloader", {})
    data_template = downloader_settings["data_folder"]

    data_folder = fill_time_template(data_template, time_obj, variable=variable)

    output_folder = Path(data_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    return output_folder


def define_raw_folder(settings, time_obj, variable=None):
    downloader_settings = settings.get("downloader", {})

    raw_template = downloader_settings.get("raw_folder", None)

    if raw_template is None:
        raw_template = str(define_data_folder(settings, time_obj, variable=variable) / "raw")

    raw_folder = fill_time_template(raw_template, time_obj, variable=variable)

    output_folder = Path(raw_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    return output_folder


def define_output_filename(settings, time_obj, variable):
    downloader_settings = settings.get("downloader", {})

    filename_template = downloader_settings.get(
        "filename_template",
        "copernicus_swi_{variable}_{time_step:%Y%m%d}.tif",
    )

    return fill_time_template(filename_template, time_obj, variable=variable)


def define_raw_filename(s3_path):
    _, key = parse_s3_path(s3_path)

    filename = Path(key).name

    if filename == "":
        filename = "copernicus_swi_raw.tif"

    if not filename.lower().endswith((".tif", ".tiff")):
        filename = filename + ".tif"

    return filename


def download_raw_file(s3_path, raw_file, s3_client, variable=None, overwrite=False):
    if raw_file.exists() and raw_file.stat().st_size > 0 and not overwrite:
        alg_logger.info(f" ----> Raw file already exists: {raw_file}")
        return raw_file

    raw_file.parent.mkdir(parents=True, exist_ok=True)

    tmp_file = raw_file.with_suffix(raw_file.suffix + ".part")

    if tmp_file.exists():
        tmp_file.unlink()

    bucket, key = parse_s3_path(s3_path)

    key_resolved = resolve_s3_object(s3_client, bucket, key, variable=variable)

    alg_logger.info(f" ----> Download raw S3 file: s3://{bucket}/{key_resolved}")

    s3_client.download_file(bucket, key_resolved, str(tmp_file))

    tmp_file.rename(raw_file)

    alg_logger.info(f" ----> Raw file saved: {raw_file}")

    return raw_file


def resolve_s3_object(s3_client, bucket, key, variable=None):
    """Resolve a catalogue s3_path to the actual GeoTIFF object.

    For SWI v4 COG catalogue rows, s3_path commonly points to a daily
    product prefix. Inside that prefix there are separate COG files for
    SWI001, SWI005, ..., QFLAG001, ..., SSF.
    """
    variable_upper = str(variable).upper() if variable is not None else None

    # If no variable is requested, allow direct object download first.
    if variable_upper is None:
        try:
            s3_client.head_object(Bucket=bucket, Key=key)
            return key
        except ClientError as exc:
            error_code = exc.response.get("Error", {}).get("Code", "")
            if error_code not in ["404", "NoSuchKey", "NotFound"]:
                raise

    prefix = key.rstrip("/") + "/"

    alg_logger.info(f" ----> Listing S3 prefix: s3://{bucket}/{prefix}")

    paginator = s3_client.get_paginator("list_objects_v2")
    candidates = []

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            obj_key = obj["Key"]
            if obj_key.endswith("/"):
                continue
            if obj_key.lower().endswith((".tif", ".tiff")):
                candidates.append(obj_key)

    if not candidates:
        raise RuntimeError(f"No GeoTIFF objects found under s3://{bucket}/{prefix}")

    if variable_upper is not None:
        def matches_variable(obj_key):
            name = Path(obj_key).name.upper()
            # Match CLMS names such as:
            # c_gls_SWI-SWI040_202606131200_GLOBE_ASCAT_V4.0.1.tiff
            # c_gls_SWI_QFLAG040_202606131200_GLOBE_ASCAT_V4.0.1.tiff
            tokens = [
                f"_{variable_upper}_",
                f"-{variable_upper}_",
                f"_{variable_upper}-",
                f"-{variable_upper}-",
                f"_{variable_upper}.",
                f"-{variable_upper}.",
            ]
            return any(token in name for token in tokens) or name.startswith(variable_upper + "_") or name.startswith(variable_upper + "-")

        variable_candidates = [obj_key for obj_key in candidates if matches_variable(obj_key)]

        if not variable_candidates:
            sample = "\n".join(["      " + Path(c).name for c in candidates[:30]])
            raise RuntimeError(
                f"Variable {variable_upper} not found under s3://{bucket}/{prefix}. "
                f"First available GeoTIFFs:\n{sample}"
            )

        # Prefer exact SWI layer, not quality flag.
        variable_candidates = [
            c for c in variable_candidates
            if "QFLAG" not in Path(c).name.upper() and "SSF" not in Path(c).name.upper()
        ] or variable_candidates

        return sorted(variable_candidates)[0]

    return sorted(candidates)[0]

def clip_local_file_to_bbox(raw_file, output_file, bbox, settings, variable, overwrite=False):
    if output_file.exists() and output_file.stat().st_size > 0 and not overwrite:
        alg_logger.info(f" ----> Skip existing: {output_file.name}")
        return "skipped"

    west, south, east, north = bbox
    tmp_file = output_file.with_suffix(output_file.suffix + ".part")

    if tmp_file.exists():
        tmp_file.unlink()

    output_settings = settings.get("output", {})
    apply_scale = bool(output_settings.get("apply_scale", True))
    scale_factor = float(output_settings.get("scale_factor", 0.5))
    source_nodata = output_settings.get("source_nodata", 255)
    target_nodata = output_settings.get("target_nodata", -9999.0)
    output_dtype = output_settings.get("dtype", "float32")

    with rasterio.open(raw_file) as src:
        window = from_bounds(
            west,
            south,
            east,
            north,
            transform=src.transform,
        )
        window = window.round_offsets().round_lengths()

        data_raw = src.read(window=window)
        transform = src.window_transform(window)

        profile = src.profile.copy()
        profile.update(
            driver="GTiff",
            height=data_raw.shape[1],
            width=data_raw.shape[2],
            transform=transform,
            compress="deflate",
            tiled=True,
        )

        if apply_scale:
            data = data_raw.astype("float32")
            data[data_raw == source_nodata] = np.nan
            data = data * scale_factor
            data = np.where(np.isnan(data), target_nodata, data).astype(output_dtype)

            profile.update(
                dtype=output_dtype,
                nodata=target_nodata,
                count=data.shape[0],
            )
        else:
            data = data_raw

        output_file.parent.mkdir(parents=True, exist_ok=True)

        with rasterio.open(tmp_file, "w", **profile) as dst:
            dst.write(data)
            dst.update_tags(
                variable=variable,
                source_file=Path(str(raw_file)).name,
                source_scale_factor=str(scale_factor),
                physical_units="percent" if apply_scale else "encoded_uint8",
                processing="bbox_clip_and_scale" if apply_scale else "bbox_clip_only",
            )
            for band_id in range(1, data.shape[0] + 1):
                dst.set_band_description(band_id, f"{variable} physical percent" if apply_scale else variable)

    tmp_file.rename(output_file)

    return "downloaded"


def download_copernicus_swi_rows(settings, rows, s3_client):
    spatial_settings = settings.get("spatial", {})
    downloader_settings = settings.get("downloader", {})

    bbox = spatial_settings.get("bounding_box", [6.0, 36.0, 19.0, 47.5])
    overwrite = downloader_settings.get("overwrite", False)
    remove_raw_file = downloader_settings.get("remove_raw_file", True)

    downloaded_files = []
    skipped_files = []

    for row in rows:
        time_step = row["time_step"]
        variable = row.get("variable", "UNKNOWN")
        s3_path = row["url"]

        output_folder = define_data_folder(settings, time_step, variable=variable)
        output_filename = define_output_filename(settings, time_step, variable=variable)
        output_file = output_folder / output_filename

        raw_folder = define_raw_folder(settings, time_step, variable=variable)
        raw_filename = f"raw_{variable}_{time_step:%Y%m%d}.tif"
        raw_file = raw_folder / raw_filename

        alg_logger.info(f" ----> Variable: {variable}")
        alg_logger.info(f" ----> Source:   {s3_path}")
        alg_logger.info(f" ----> Raw:      {raw_file}")
        alg_logger.info(f" ----> Target:   {output_file}")

        status = None

        try:
            download_raw_file(
                s3_path=s3_path,
                raw_file=raw_file,
                s3_client=s3_client,
                variable=variable,
                overwrite=overwrite,
            )

            status = clip_local_file_to_bbox(
                raw_file=raw_file,
                output_file=output_file,
                bbox=bbox,
                settings=settings,
                variable=variable,
                overwrite=overwrite,
            )

            if remove_raw_file and raw_file.exists():
                raw_file.unlink()
                alg_logger.info(f" ----> Removed raw file: {raw_file}")

        except Exception as exc:
            alg_logger.error(f" ===> Error downloading {s3_path}: {exc}")
            continue

        if status == "skipped":
            skipped_files.append(str(output_file))
        else:
            downloaded_files.append(str(output_file))
            alg_logger.info(f" ----> Downloaded: {output_file.name}")

    return downloaded_files, skipped_files


def download_copernicus_swi_date_range(settings, time_start, time_end, s3_client):
    product_settings = settings.get("product", {})
    product_name = product_settings.get("name", "swi_global_12.5km_daily_v4")
    variables = product_settings.get("variables", [])

    total_found = 0
    total_downloaded = 0
    total_skipped = 0

    alg_logger.info(" ----> Running date-range mode")
    alg_logger.info(f" ::: Product: {product_name}")
    alg_logger.info(f" ::: Variables requested: {variables}")
    alg_logger.info(f' ::: Time start: {time_start.strftime("%Y-%m-%d")}')
    alg_logger.info(f' ::: Time end  : {time_end.strftime("%Y-%m-%d")}')

    df_catalogue, url_column = read_catalogue(settings)

    for day in daterange(time_start, time_end):
        day_str = day.strftime("%Y-%m-%d")

        alg_logger.info(" ")
        alg_logger.info(f" ::: DAY: {day_str}")
        alg_logger.info(f" ::: Searching {product_name} files ...")

        rows = select_catalogue_rows(df_catalogue, url_column, day, settings)

        n_found = len(rows)
        total_found += n_found

        alg_logger.info(f" ::: Files found after variable filter: {n_found}")

        if not rows:
            continue

        downloaded_files, skipped_files = download_copernicus_swi_rows(
            settings=settings,
            rows=rows,
            s3_client=s3_client,
        )

        total_downloaded += len(downloaded_files)
        total_skipped += len(skipped_files)

    alg_logger.info(" ")
    alg_logger.info(" ::: DOWNLOAD SUMMARY")
    alg_logger.info(f" ::: Total files found      : {total_found}")
    alg_logger.info(f" ::: Total downloaded       : {total_downloaded}")
    alg_logger.info(f" ::: Total skipped existing : {total_skipped}")


# ----------------------------------------------------------------------------------------------------------------------
# algorithm main
def main():
    global alg_name, alg_version, alg_release

    args = get_args()
    settings = read_file_json(args.settings_file)
    get_logger(settings)

    alg_info = settings.get("algorithm", {})
    alg_name = alg_info.get("name", alg_name)
    alg_version = alg_info.get("version", alg_version)
    alg_release = alg_info.get("release", alg_release)

    alg_logger.info(" ============================================================================ ")
    alg_logger.info(f" ==> {alg_name} (Version: {alg_version} Release_Date: {alg_release})")
    alg_logger.info(" ==> START ... ")
    alg_logger.info(" ")
    alg_logger.info(" ---> Settings file: " + args.settings_file)

    start_time = time.time()

    atexit.register(release_lock)
    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)

    alg_logger.info(" ---> Select download mode ... ")

    try:
        mode, reference_time, time_start, time_end, n_days = select_run_mode(
            args,
            settings,
        )
    except Exception as exc:
        alg_logger.error(f" ===> Error in selecting download mode: {exc}")
        alg_logger.info(" ---> Select download mode ... FAILED")
        sys.exit(1)

    alg_logger.info(" ::: Reference time: " + reference_time.strftime("%Y-%m-%d %H:%M:%S"))
    alg_logger.info(" ::: Download mode: " + mode)
    alg_logger.info(" ::: Time start: " + time_start.strftime("%Y-%m-%d"))
    alg_logger.info(" ::: Time end:   " + time_end.strftime("%Y-%m-%d"))

    if mode == "rolling_window":
        alg_logger.info(" ::: Days back: " + str(n_days))

    alg_logger.info(" ---> Select download mode ... DONE")

    alg_logger.info(" ---> Initialize download mode ... ")
    alg_logger.info(" ")

    if not acquire_lock(settings=settings, force_lock=args.force_lock):
        sys.exit(1)

    try:
        s3_client = get_s3_client()
        alg_logger.info(" ----> S3 client initialized")
    except Exception as exc:
        alg_logger.error(f" ===> S3 initialization failed: {exc}")
        release_lock()
        sys.exit(1)

    alg_logger.info(" ---> Initialize download mode ... DONE")

    alg_logger.info(" ---> Execute download mode ... ")
    alg_logger.info(" ::: Timestamp Start -- " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    try:
        if time_start > time_end:
            raise RuntimeError("time_start must be less than or equal to time_end")

        download_copernicus_swi_date_range(
            settings=settings,
            time_start=time_start,
            time_end=time_end,
            s3_client=s3_client,
        )

        alg_logger.info(" ::: Timestamp End -- " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        alg_logger.info(" ---> Execute download mode ... DONE")

    except Exception as exc:
        alg_logger.error(f" ===> ERROR in executing download algorithm: {exc}")
        alg_logger.info(" ---> Execute download mode ... FAILED")
        sys.exit(1)

    alg_time_elapsed = round(time.time() - start_time, 1)

    alg_logger.info(" ")
    alg_logger.info(f" ==> {alg_name} (Version: {alg_version} Release_Date: {alg_release})")
    alg_logger.info(" ==> TIME ELAPSED: " + str(alg_time_elapsed) + " seconds")
    alg_logger.info(" ==> ... END")
    alg_logger.info(" ==> Bye, Bye")
    alg_logger.info(" ============================================================================ ")


if __name__ == "__main__":
    main()
