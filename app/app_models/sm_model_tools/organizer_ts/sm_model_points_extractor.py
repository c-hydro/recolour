#!/usr/bin/env python3
"""
SM model point extractor.

Reads hourly NetCDF variables and daily GeoTIFF SWI/SM values over a backward time
period and writes one CSV time-series file per output variable.

Output format:
    time,ID1,ID2,ID3,...
    2026061012,value,value,value,...

Missing input file value: -9998
Invalid/no valid data value: -9999

Main features:
  - hourly NetCDF extraction by radius mean
  - daily GeoTIFF extraction by nearest neighbour or radius mean
  - optional radius extension if the initial radius has no valid value
  - detailed logs for file read, missing files, read failures and radius extension

Example:
    python sm_model_points_extractor.py \
        -settings_file sm_model_points_extractor.json \
        -time "2026-06-10 12:00"
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import xarray as xr
except Exception:  # pragma: no cover
    xr = None

try:
    import h5py
except Exception:  # pragma: no cover
    h5py = None

try:
    import rasterio
    from rasterio.warp import transform as rio_transform
except Exception:  # pragma: no cover
    rasterio = None
    rio_transform = None

try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover
    cKDTree = None

EARTH_RADIUS_M = 6371000.0
VALUE_MISSING_FILE = -9998.0
VALUE_INVALID = -9999.0


@dataclass
class Point:
    point_id: str
    lon: float
    lat: float
    attrs: Dict[str, Any]


def log_info(message: str) -> None:
    print(f"INFO  --> {message}")


def log_warning(message: str) -> None:
    print(f"WARN  --> {message}")


def parse_time(value: str, tz: timezone = timezone.utc) -> datetime:
    value = str(value).strip()
    if value.lower() == "now":
        return datetime.now(tz).replace(minute=0, second=0, microsecond=0)

    formats = (
        "%Y%m%d%H",
        "%Y%m%d%H%M",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M",
        "%Y-%m-%dT%H:%M:%S",
    )
    for fmt in formats:
        try:
            return datetime.strptime(value, fmt).replace(tzinfo=tz)
        except ValueError:
            pass

    dt = datetime.fromisoformat(value)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=tz)
    return dt.astimezone(tz)


def build_time_steps(time_start: datetime, time_end: datetime, step_hours: int, n_steps: Optional[int]) -> List[datetime]:
    if step_hours <= 0:
        raise ValueError("time.step_hours must be > 0")

    direction = -1 if time_start >= time_end else 1
    delta = timedelta(hours=step_hours * direction)
    steps: List[datetime] = []
    current = time_start

    while current >= time_end if direction < 0 else current <= time_end:
        steps.append(current)
        if n_steps is not None and len(steps) >= n_steps:
            break
        current = current + delta

    return steps


def read_points(path: str, cfg: Dict[str, Any]) -> List[Point]:
    delimiter = cfg.get("delimiter", "auto")
    lon_col = cfg.get("lon_col", "Lon")
    lat_col = cfg.get("lat_col", "Lat")
    id_col = cfg.get("id_col", "ID")
    encoding = cfg.get("encoding", "utf-8-sig")

    text = Path(path).read_text(encoding=encoding)
    first_line = next((line for line in text.splitlines() if line.strip()), "")

    if delimiter == "auto":
        candidates = [";", ",", "\t", " "]
        delimiter = max(candidates, key=lambda d: first_line.count(d))
        if delimiter == " ":
            delimiter = None

    if delimiter is None:
        reader = csv.DictReader(text.splitlines(), delimiter=" ", skipinitialspace=True)
    else:
        reader = csv.DictReader(text.splitlines(), delimiter=delimiter)

    points: List[Point] = []
    for idx, row in enumerate(reader):
        clean = {str(k).strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k is not None}
        if not clean:
            continue
        if lon_col not in clean or lat_col not in clean:
            raise KeyError(f"Point file must contain columns '{lon_col}' and '{lat_col}'. Found: {list(clean)}")
        point_id = str(clean.get(id_col, idx + 1)).strip()
        points.append(Point(point_id=point_id, lon=float(clean[lon_col]), lat=float(clean[lat_col]), attrs=clean))

    if not points:
        raise RuntimeError(f"No points found in {path}")
    return points


def format_path(template: str, t: datetime, time_now: datetime) -> str:
    return template.format(
        time=t,
        time_now=time_now,
        yyyy=t.strftime("%Y"),
        mm=t.strftime("%m"),
        dd=t.strftime("%d"),
        yyyymmdd=t.strftime("%Y%m%d"),
        yyyymmddhh=t.strftime("%Y%m%d%H"),
    )


def lonlat_to_unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])


def radius_to_chord(radius_m: float) -> float:
    return 2.0 * math.sin(radius_m / (2.0 * EARTH_RADIUS_M))


def great_circle_distance_m(lon0: float, lat0: float, lon1: np.ndarray, lat1: np.ndarray) -> np.ndarray:
    lon0r = math.radians(float(lon0))
    lat0r = math.radians(float(lat0))
    lon1r = np.deg2rad(lon1)
    lat1r = np.deg2rad(lat1)
    dlon = lon1r - lon0r
    dlat = lat1r - lat0r
    a = np.sin(dlat / 2.0) ** 2 + math.cos(lat0r) * np.cos(lat1r) * np.sin(dlon / 2.0) ** 2
    return EARTH_RADIUS_M * 2.0 * np.arcsin(np.sqrt(a))


def apply_valid_rules(values: np.ndarray, cfg: Dict[str, Any]) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64).copy()

    nodata_values = cfg.get("nodata_values", cfg.get("nodata", []))
    if nodata_values is None:
        nodata_values = []
    if not isinstance(nodata_values, list):
        nodata_values = [nodata_values]
    for nodata in nodata_values:
        out[np.isclose(out, float(nodata), equal_nan=False)] = np.nan

    valid_min = cfg.get("valid_min", None)
    valid_max = cfg.get("valid_max", None)
    if valid_min is not None:
        out[out < float(valid_min)] = np.nan
    if valid_max is not None:
        out[out > float(valid_max)] = np.nan

    return out


def get_radius_sequence(radius_m: float, cfg: Dict[str, Any]) -> List[float]:
    """Return radius attempts based on extend options."""
    radius_m = float(radius_m)
    if not cfg.get("extend_radius_if_empty", False):
        return [radius_m]

    factor = float(cfg.get("extend_radius_factor", 2.0))
    max_radius = float(cfg.get("extend_radius_max", radius_m))

    if factor <= 1.0:
        raise ValueError("extend_radius_factor must be > 1")
    if max_radius < radius_m:
        max_radius = radius_m

    radii = [radius_m]
    current = radius_m
    while current < max_radius:
        current = min(current * factor, max_radius)
        if current > radii[-1]:
            radii.append(current)
        else:
            break
    return radii


def mean_by_radius(
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    value_grid: np.ndarray,
    points: Sequence[Point],
    radius_m: float,
    var_cfg: Dict[str, Any],
    var_key: str = "variable",
) -> List[float]:
    values = apply_valid_rules(value_grid, var_cfg)
    lon_arr = np.asarray(lon_grid, dtype=np.float64)
    lat_arr = np.asarray(lat_grid, dtype=np.float64)

    if values.ndim > 2:
        values = np.squeeze(values).astype(np.float64, copy=False)

    if lon_arr.ndim == 1 and lat_arr.ndim == 1 and values.ndim >= 2:
        lon_arr, lat_arr = np.meshgrid(lon_arr, lat_arr)

    valid = np.isfinite(lon_arr) & np.isfinite(lat_arr) & np.isfinite(values)
    if not np.any(valid):
        return [VALUE_INVALID for _ in points]

    lon_flat = lon_arr[valid]
    lat_flat = lat_arr[valid]
    val_flat = values[valid].astype(np.float64, copy=False)
    radii = get_radius_sequence(radius_m, var_cfg)

    output: List[float] = []
    if cKDTree is not None:
        xyz = lonlat_to_unit_xyz(lon_flat, lat_flat)
        tree = cKDTree(xyz)
        for p in points:
            p_xyz = lonlat_to_unit_xyz(np.array([p.lon], dtype=np.float64), np.array([p.lat], dtype=np.float64))[0]
            found_value = VALUE_INVALID
            used_radius = None
            for radius_current in radii:
                idxs = tree.query_ball_point(p_xyz, radius_to_chord(radius_current))
                if idxs:
                    value = float(np.nanmean(val_flat[idxs]))
                    if np.isfinite(value):
                        found_value = value
                        used_radius = radius_current
                        break
            if used_radius is not None and used_radius > radius_m:
                log_warning(
                    f"POINT {p.point_id} | {var_key}: radius extended from {radius_m:.0f} to {used_radius:.0f} m"
                )
            output.append(found_value)
        return output

    for p in points:
        dist = great_circle_distance_m(p.lon, p.lat, lon_flat, lat_flat)
        found_value = VALUE_INVALID
        used_radius = None
        for radius_current in radii:
            idxs = np.where(dist <= radius_current)[0]
            if idxs.size > 0:
                value = float(np.nanmean(val_flat[idxs]))
                if np.isfinite(value):
                    found_value = value
                    used_radius = radius_current
                    break
        if used_radius is not None and used_radius > radius_m:
            log_warning(
                f"POINT {p.point_id} | {var_key}: radius extended from {radius_m:.0f} to {used_radius:.0f} m"
            )
        output.append(found_value)
    return output


def read_netcdf_variable(path: str, variable_name: str, lon_name: str, lat_name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read NetCDF arrays as float64 without automatic CF masking/scaling."""
    if xr is not None:
        try:
            with xr.open_dataset(path, decode_cf=False, mask_and_scale=False) as ds:
                if variable_name not in ds:
                    raise KeyError(f"Variable '{variable_name}' not found. Available: {list(ds.data_vars)}")
                if lon_name not in ds:
                    raise KeyError(f"Longitude '{lon_name}' not found. Available: {list(ds.variables)}")
                if lat_name not in ds:
                    raise KeyError(f"Latitude '{lat_name}' not found. Available: {list(ds.variables)}")

                data = np.asarray(ds[variable_name].squeeze().values, dtype=np.float64)
                lon = np.asarray(ds[lon_name].squeeze().values, dtype=np.float64)
                lat = np.asarray(ds[lat_name].squeeze().values, dtype=np.float64)
                return lon, lat, data
        except Exception as exc:
            log_warning(f"XARRAY READ FAILED, TRY H5PY: {exc}")

    if h5py is None:
        raise ImportError("Install xarray/netCDF4/h5netcdf or h5py to read NetCDF files")

    with h5py.File(path, "r") as fp:
        if variable_name not in fp:
            raise KeyError(f"Variable '{variable_name}' not found in {path}")
        data = np.asarray(fp[variable_name][()], dtype=np.float64).squeeze()
        lon = np.asarray(fp[lon_name][()], dtype=np.float64).squeeze()
        lat = np.asarray(fp[lat_name][()], dtype=np.float64).squeeze()
    return lon, lat, data


def extract_netcdf_values(path: str, var_cfg: Dict[str, Any], points: Sequence[Point], global_radius_m: float, var_key: str) -> List[float]:
    variable_name = var_cfg["variable_name"]
    lon_name = var_cfg.get("lon_name", "longitude")
    lat_name = var_cfg.get("lat_name", "latitude")
    radius_m = float(var_cfg.get("radius_m", global_radius_m))

    lon, lat, data = read_netcdf_variable(path, variable_name, lon_name, lat_name)
    return mean_by_radius(lon, lat, data, points, radius_m, var_cfg, var_key=var_key)


def tiff_sample_nn(src: Any, data: np.ndarray, point: Point, var_cfg: Dict[str, Any]) -> float:
    lon, lat = point.lon, point.lat
    if src.crs and str(src.crs).upper() not in ("EPSG:4326", "OGC:CRS84") and rio_transform is not None:
        xs, ys = rio_transform("EPSG:4326", src.crs, [lon], [lat])
        x0, y0 = xs[0], ys[0]
    else:
        x0, y0 = lon, lat

    try:
        row, col = src.index(x0, y0)
    except Exception:
        return VALUE_INVALID

    if row < 0 or row >= src.height or col < 0 or col >= src.width:
        return VALUE_INVALID

    value = float(data[row, col])
    checked = apply_valid_rules(np.array([value], dtype=np.float64), var_cfg)[0]
    if not np.isfinite(checked):
        return VALUE_INVALID
    return float(checked)


def tiff_mean_radius(src: Any, data: np.ndarray, point: Point, radius_m: float, var_cfg: Dict[str, Any], var_key: str) -> float:
    lon, lat = point.lon, point.lat
    if src.crs and str(src.crs).upper() not in ("EPSG:4326", "OGC:CRS84") and rio_transform is not None:
        xs, ys = rio_transform("EPSG:4326", src.crs, [lon], [lat])
        x0, y0 = xs[0], ys[0]
        base_pixel_size = max(abs(src.transform.a), abs(src.transform.e))
    else:
        x0, y0 = lon, lat
        deg_lat = radius_m / 111320.0
        deg_lon = radius_m / (111320.0 * max(math.cos(math.radians(lat)), 0.1))
        base_pixel_size = max(abs(src.transform.a), abs(src.transform.e), 1.0e-12)

    try:
        row0, col0 = src.index(x0, y0)
    except Exception:
        return VALUE_INVALID

    radii = get_radius_sequence(radius_m, var_cfg)
    for radius_current in radii:
        if src.crs and str(src.crs).upper() not in ("EPSG:4326", "OGC:CRS84") and rio_transform is not None:
            pixel_radius = max(1, int(math.ceil(radius_current / base_pixel_size)))
        else:
            deg_lat = radius_current / 111320.0
            deg_lon = radius_current / (111320.0 * max(math.cos(math.radians(lat)), 0.1))
            pixel_radius = max(1, int(math.ceil(max(deg_lon / abs(src.transform.a), deg_lat / abs(src.transform.e)))))

        rows = np.arange(max(0, row0 - pixel_radius), min(src.height, row0 + pixel_radius + 1))
        cols = np.arange(max(0, col0 - pixel_radius), min(src.width, col0 + pixel_radius + 1))
        if rows.size == 0 or cols.size == 0:
            continue

        rr, cc = np.meshgrid(rows, cols, indexing="ij")
        xs, ys = rasterio.transform.xy(src.transform, rr, cc, offset="center")
        xs = np.asarray(xs, dtype=np.float64)
        ys = np.asarray(ys, dtype=np.float64)

        if src.crs and str(src.crs).upper() not in ("EPSG:4326", "OGC:CRS84") and rio_transform is not None:
            lons, lats = rio_transform(src.crs, "EPSG:4326", xs.ravel().tolist(), ys.ravel().tolist())
            lons = np.asarray(lons, dtype=np.float64).reshape(xs.shape)
            lats = np.asarray(lats, dtype=np.float64).reshape(ys.shape)
        else:
            lons, lats = xs, ys
		
        dist = great_circle_distance_m(point.lon, point.lat, lons, lats)
        mask = np.asarray(dist <= radius_current, dtype=bool).reshape(rr.shape)
        #ask = np.asarray(dist <= radius_current, dtype=bool)

        window_vals = np.asarray(data[rr, cc], dtype=np.float64)

        if window_vals.shape != mask.shape:
            window_vals = np.squeeze(window_vals)
            mask = np.squeeze(mask)

        if window_vals.shape != mask.shape:
            raise RuntimeError(
                f"GeoTIFF window/mask shape mismatch: "
                f"data={data.shape}, rr={rr.shape}, cc={cc.shape}, "
                f"window={window_vals.shape}, mask={mask.shape}"
        )

        vals = window_vals[mask]
        vals = apply_valid_rules(vals, var_cfg)
        vals = vals[np.isfinite(vals)]
		
        if vals.size > 0:
            if radius_current > radius_m:
                log_warning(
                    f"POINT {point.point_id} | {var_key}: radius extended from {radius_m:.0f} to {radius_current:.0f} m"
                )
            return float(np.nanmean(vals))
        
    if var_cfg.get("nn_if_invalid", True):
        return tiff_sample_nn(src, data, point, var_cfg)
    return VALUE_INVALID


def extract_tiff_values(path: str, var_cfg: Dict[str, Any], points: Sequence[Point], global_radius_m: float, var_key: str) -> List[float]:
    if rasterio is None:
        raise ImportError("Install rasterio to read GeoTIFF files")

    radius_m = float(var_cfg.get("radius_m", global_radius_m))
    band = int(var_cfg.get("band", 1))
    method = str(var_cfg.get("method", "nn")).lower()

    with rasterio.open(path) as src:
        data = src.read(band).astype(np.float64, copy=False)
        nodata = src.nodata
        if nodata is not None:
            data = data.copy()
            data[np.isclose(data, float(nodata), equal_nan=False)] = np.nan
        if method == "radius_mean":
            return [tiff_mean_radius(src, data, p, radius_m, var_cfg, var_key=var_key) for p in points]
        return [tiff_sample_nn(src, data, p, var_cfg) for p in points]


def write_variable_csv(path: str, rows: Sequence[List[Any]], header: Sequence[str], delimiter: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp, delimiter=delimiter)
        writer.writerow(header)
        writer.writerows(rows)


def format_value(value: float, precision: int) -> Any:
    value = float(value)
    if np.isclose(value, VALUE_MISSING_FILE) or np.isclose(value, VALUE_INVALID):
        return int(value)
    if not np.isfinite(value):
        return int(VALUE_INVALID)
    return round(value, precision)


def extract_variable_series(
    var_key: str,
    var_cfg: Dict[str, Any],
    steps: Sequence[datetime],
    points: Sequence[Point],
    time_now: datetime,
    global_radius_m: float,
    delimiter: str,
    output_cfg: Dict[str, Any],
) -> None:
    if not var_cfg.get("enabled", True):
        log_info(f"VARIABLE {var_key}: disabled")
        return

    variable_type = var_cfg.get("type", "netcdf")
    file_template = var_cfg["file_template"]
    missing_file_value = float(var_cfg.get("missing_file_value", VALUE_MISSING_FILE))
    precision = int(var_cfg.get("precision", output_cfg.get("precision", 3)))

    point_ids = [p.point_id for p in points]
    header = ["time"] + point_ids
    rows: List[List[Any]] = []

    log_info(f"VARIABLE {var_key}: start extraction")
    for t in steps:
        t_tag = t.strftime("%Y%m%d%H")
        
        # Daily GeoTIFF variables: use only the default daily hour
        # when filename has no hourly token.
        if variable_type == "tiff" and var_cfg.get("frequency", "daily") == "daily":
            daily_hour = int(var_cfg.get("daily_hour", 0))
            if t.hour != daily_hour:
                continue
        
        file_path = format_path(file_template, t, time_now)
        log_info(f"TIME {t_tag} | {var_key} | file: {file_path}")

        if not os.path.exists(file_path):
            log_warning(f"TIME {t_tag} | {var_key} | FILE NOT FOUND: {file_path}")
            values = [missing_file_value for _ in points]
        else:
            try:
                if variable_type == "netcdf":
                    values = extract_netcdf_values(file_path, var_cfg, points, global_radius_m, var_key=var_key)
                elif variable_type == "tiff":
                    values = extract_tiff_values(file_path, var_cfg, points, global_radius_m, var_key=var_key)
                else:
                    raise ValueError(f"Unsupported variable type '{variable_type}' for {var_key}")
            except Exception as exc:
                if var_cfg.get("raise_on_error", False):
                    raise
                log_warning(f"TIME {t_tag} | {var_key} | READ FAILED: {exc}")
                values = [VALUE_INVALID for _ in points]

        rows.append([t_tag] + [format_value(v, precision) for v in values])

    output_folder = output_cfg.get("folder_name", "./output")
    output_folder = output_folder.format(time_now=time_now)
    file_template_out = output_cfg.get("file_template", "obs_db_{var_name}_{time_now:%Y%m%d%H}.csv")
    out_path = os.path.join(output_folder, file_template_out.format(var_name=var_key, time_now=time_now))
    write_variable_csv(out_path, rows, header, delimiter)
    log_info(f"VARIABLE {var_key}: wrote {len(rows)} rows to {out_path}")


def normalize_variables(cfg: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    if "variables" in cfg:
        return cfg["variables"]

    variables: Dict[str, Dict[str, Any]] = {}
    for key in ("rain", "t", "temperature", "sm", "swi"):
        if key in cfg:
            out_key = "t" if key == "temperature" else ("sm" if key == "swi" else key)
            variables[out_key] = cfg[key]
    return variables


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract point time series from NetCDF and GeoTIFF files")
    parser.add_argument("-settings_file", "--settings_file", required=True, help="JSON settings file")
    parser.add_argument("-time", "--time", dest="time_now", default=None, help="Reference time/time_now")
    parser.add_argument("-period_hours", "--period_hours", type=int, default=None, help="Backward period in hours")
    args = parser.parse_args()

    cfg = json.loads(Path(args.settings_file).read_text())
    tz = timezone.utc
    time_cfg = cfg.get("time", {})

    time_now = parse_time(args.time_now if args.time_now else cfg.get("time_now", time_cfg.get("time_now", "now")), tz=tz)
    time_start_raw = time_cfg.get("time_start", "now")
    time_start = time_now if str(time_start_raw).lower() == "now" else parse_time(str(time_start_raw), tz=tz)

    step_hours = int(time_cfg.get("step_hours", 1))
    period_hours = int(args.period_hours if args.period_hours is not None else time_cfg.get("period_hours", 24 * 10))
    if time_cfg.get("time_end") not in (None, ""):
        time_end = parse_time(str(time_cfg["time_end"]), tz=tz)
    else:
        time_end = time_start - timedelta(hours=period_hours)

    n_steps = time_cfg.get("n_steps")
    n_steps = int(n_steps) if n_steps is not None else None
    steps = build_time_steps(time_start, time_end, step_hours, n_steps)

    log_info("POINTS EXTRACTOR - START")
    log_info(f"SETTINGS FILE: {args.settings_file}")
    log_info(f"TIME NOW   : {time_now:%Y-%m-%d %H:%M}")
    log_info(f"TIME START : {time_start:%Y-%m-%d %H:%M}")
    log_info(f"TIME END   : {time_end:%Y-%m-%d %H:%M}")
    log_info(f"STEP HOURS : {step_hours}")
    log_info(f"N STEPS    : {len(steps)}")

    points = read_points(cfg["points"]["file"], cfg["points"])
    log_info(f"POINTS     : {len(points)} from {cfg['points']['file']}")

    output_cfg = cfg.get("output", {})
    delimiter = output_cfg.get("delimiter", ",")
    global_radius_m = float(cfg.get("radius_m", 1000))

    variables = normalize_variables(cfg)
    if not variables:
        raise RuntimeError("No variables configured. Use a 'variables' block or rain/t/sm blocks.")

    for var_key, var_cfg in variables.items():
        extract_variable_series(
            var_key=var_key,
            var_cfg=var_cfg,
            steps=steps,
            points=points,
            time_now=time_now,
            global_radius_m=global_radius_m,
            delimiter=delimiter,
            output_cfg=output_cfg,
        )

    log_info("POINTS EXTRACTOR - END")


if __name__ == "__main__":
    main()
