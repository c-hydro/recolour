#!/usr/bin/env python3
"""
Filter NetCDF files by whether they have any points inside a shapefile domain.
Now with detailed step-by-step prints for debugging and execution tracing.

New:
- Report file now includes a timestamp in its filename: filter_report_YYYYMMDDHHMM.csv
- Optional JSON key "report_dir" to control where the report is saved.
  If omitted, the report is saved under the destination root (dst).
  If "report_dir" is relative, it is resolved relative to dst.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from typing import List, Optional, Tuple
from datetime import datetime
import re
import numpy as np
import xarray as xr
import fiona
from shapely.geometry import shape, Point
from shapely.ops import unary_union, transform as shp_transform
from shapely.prepared import prep
from pyproj import CRS, Transformer
import pandas as pd


def path_has_date_in_range(path, start_time, end_time, date_format="%Y%m%d%H%M"):
    """
    Checks if the *path string* contains a date string between start_time and end_time (half-open).
    - path: string containing the date somewhere (e.g., ".../2025/07/01/file_202507010000.nc")
    - start_time, end_time: datetime-like (or pandas Timestamp)
    - date_format: strptime format matching the 12-digit date in the string

    Returns: True if date is within [start_time, end_time), else False.
    """
    if path is None:
        return False
    path = str(path)
    match = re.search(r"\d{12}", path)  # matches YYYYMMDDHHMM anywhere in the path
    if not match:
        return False
    file_date = datetime.strptime(match.group(), date_format)
    return start_time <= file_date < end_time


def filename_has_date_in_range(path, start_time, end_time, date_format="%Y%m%d%H%M"):
    """
    Like path_has_date_in_range, but only checks the basename (filename).
    """
    if path is None:
        return False
    fname = Path(str(path)).name
    match = re.search(r"\d{12}", fname)
    if not match:
        return False
    file_date = datetime.strptime(match.group(), date_format)
    return start_time <= file_date < end_time


def load_domain_geometry(shp_path: Path, target_crs: str = "EPSG:4326"):
    print(f"[INFO] Loading shapefile: {shp_path}")
    shp_path = Path(shp_path)
    if not shp_path.exists():
        raise FileNotFoundError(f"Shapefile not found: {shp_path}")
    with fiona.open(shp_path) as src:
        src_crs = None
        try:
            if src.crs:
                src_crs = CRS.from_user_input(src.crs)
                print(f"[INFO] Shapefile CRS: {src_crs}")
        except Exception:
            src_crs = None
        geoms = [shape(feat["geometry"]) for feat in src]
    print(f"[INFO] Loaded {len(geoms)} geometry(ies) from shapefile.")

    if not geoms:
        raise ValueError("No geometries found in shapefile.")

    target = CRS.from_user_input(target_crs)
    if src_crs is not None and src_crs != target:
        print(f"[INFO] Reprojecting geometries to target CRS: {target_crs}")
        transformer = Transformer.from_crs(src_crs, target, always_xy=True)

        def _proj(x, y, z=None):
            x2, y2 = transformer.transform(x, y)
            return (x2, y2) if z is None else (x2, y2, z)

        geoms = [shp_transform(_proj, g) for g in geoms]

    union_poly = unary_union(geoms)
    print("[INFO] Domain geometry prepared for fast point-in-polygon checks.")
    return prep(union_poly)


def find_lat_lon_vars(ds: xr.Dataset) -> Tuple[str, str]:
    print(f"[DEBUG] Searching for latitude/longitude variables in dataset.")
    cand_lat, cand_lon = None, None
    for v in ds.variables:
        std = str(ds[v].attrs.get("standard_name", "")).lower()
        if std == "latitude":
            cand_lat = v
        elif std == "longitude":
            cand_lon = v
    if cand_lat and cand_lon:
        print(f"[DEBUG] Found lat/lon via standard_name: {cand_lat}, {cand_lon}")
        return cand_lat, cand_lon

    def pick(cands: List[str]):
        for v in ds.variables:
            vn = v.lower()
            for c in cands:
                if c in vn and ds[v].ndim >= 1 and ds[v].size > 1:
                    return v
        return None

    lat_var = pick(["latitude", "lat"])
    lon_var = pick(["longitude", "lon", "long"])
    if not lat_var or not lon_var:
        raise KeyError("Could not find latitude/longitude variables in dataset.")
    print(f"[DEBUG] Found lat/lon via name pattern: {lat_var}, {lon_var}")
    return lat_var, lon_var


def dataset_time_window(ds: xr.Dataset):
    cands = []
    if "time" in ds.variables:
        cands.append("time")
    for coord in getattr(ds, "coords", {}):
        if coord == "time":
            cands.append(coord)
    if not cands:
        return None
    for tvar in dict.fromkeys(cands):
        try:
            tvals = pd.to_datetime(ds[tvar].values)
            if tvals.size == 0:
                continue
            print(f"[DEBUG] Dataset time range: {tvals.min()} to {tvals.max()}")
            return pd.to_datetime(np.min(tvals)), pd.to_datetime(np.max(tvals))
        except Exception:
            continue
    return None


def file_in_range(nc_path: Path, start_date: Optional[pd.Timestamp], end_date: Optional[pd.Timestamp], mode: str = "auto") -> bool:
    """
    Returns True if the file should be considered based on the time filter.

    Modes:
      - auto: try path/filename first, then open NetCDF time if needed
      - filename: only check the filename
      - netcdf: only check the dataset's time variable
      - off: ignore any time filtering (always True)
    """
    if start_date is None or end_date is None or mode == "off":
        return True

    mode = mode.lower()
    print(f"[INFO] Checking date range for file: {nc_path.name}")

    if mode == "auto" and path_has_date_in_range(nc_path, start_date, end_date):
        print("[DEBUG] Date matched by folder structure or path.")
        return True
    if mode == "filename" and filename_has_date_in_range(nc_path, start_date, end_date):
        print("[DEBUG] Date matched by filename.")
        return True
    if mode in {"auto", "netcdf"}:
        try:
            print("[DEBUG] Opening dataset to check time variable.")
            ds = xr.open_dataset(nc_path)
            win = dataset_time_window(ds)
        except Exception:
            win = None
        finally:
            try:
                ds.close()
            except Exception:
                pass
        if win is not None:
            tmin, tmax = win
            # Overlap test on normalized days (inclusive on both sides here)
            overlap = not (tmax.normalize() < start_date.normalize() or tmin.normalize() > end_date.normalize())
            print(f"[DEBUG] NetCDF time match: {overlap}")
            return overlap
        if mode == "netcdf":
            return False

    # Fallback to checking the path/filename if auto
    if path_has_date_in_range(nc_path, start_date, end_date):
        print("[DEBUG] Date matched by folder structure (fallback).")
        return True
    return filename_has_date_in_range(nc_path, start_date, end_date)


def points_inside_domain(ds: xr.Dataset, prepared_domain, lat_var: str, lon_var: str,
                         max_checks: Optional[int] = None) -> int:
    lat = ds[lat_var].values
    lon = ds[lon_var].values

    if lat.ndim == 1 and lon.ndim == 1 and lat.shape == lon.shape:
        lat_flat = lat
        lon_flat = lon
    elif lat.ndim == 2 and lon.ndim == 2 and lat.shape == lon.shape:
        lat_flat = lat.ravel()
        lon_flat = lon.ravel()
    else:
        if lat.ndim == 1 and lon.ndim == 1:
            lon_mesh, lat_mesh = np.meshgrid(lon, lat)
            lat_flat = lat_mesh.ravel()
            lon_flat = lon_mesh.ravel()
        else:
            raise ValueError(f"Unsupported lat/lon shapes: {lat.shape} vs {lon.shape}")

    n = lat_flat.size
    indices = np.arange(n)

    if max_checks is not None:
        try:
            max_checks = int(max_checks)
        except Exception:
            max_checks = None

    if max_checks is not None and max_checks < n:
        print(f"[DEBUG] Sampling {max_checks} points from {n} total grid points.")
        rng = np.random.default_rng(123)
        indices = rng.choice(n, size=max_checks, replace=False)

    count = 0
    for i in indices:
        pt = Point(float(lon_flat[i]), float(lat_flat[i]))
        if prepared_domain.contains(pt):
            count += 1
    print(f"[DEBUG] Found {count} points inside domain.")
    return count


def file_is_in_domain(nc_path: Path, prepared_domain, min_points_inside: int = 1,
                      sample_limit: Optional[int] = None) -> bool:
    print(f"[INFO] Checking domain inclusion for: {nc_path.name}")
    try:
        ds = xr.open_dataset(nc_path)
    except Exception as e:
        print(f"[WARN] Could not open {nc_path.name}: {e}")
        return False
    try:
        lat_var, lon_var = find_lat_lon_vars(ds)
        inside = points_inside_domain(ds, prepared_domain, lat_var, lon_var, max_checks=sample_limit)
        result = inside >= int(min_points_inside)
        print(f"[INFO] Domain check result: {result} ({inside} points inside)")
        return result
    except Exception as e:
        print(f"[WARN] {nc_path.name}: {e}")
        return False
    finally:
        try:
            ds.close()
        except Exception:
            pass


def merge_settings(json_cfg: dict, cli_args: dict) -> dict:
    out = dict(json_cfg) if json_cfg else {}
    for k, v in cli_args.items():
        if v is not None:
            out[k] = v
    print(f"[INFO] Merged settings: {out}")
    return out


def main():
    ap = argparse.ArgumentParser(description="Filter NetCDF files by shapefile domain")
    ap.add_argument("-settings_file", type=str, help="Path to JSON settings file")
    ap.add_argument("-shapefile", type=str, help="Path to .shp")
    ap.add_argument("-src", type=str, help="Source folder root (searched recursively)")
    ap.add_argument("-dst", type=str, help="Destination root (preserve relative structure)")
    ap.add_argument("-pattern", type=str, default=None, help='Glob pattern (e.g., "*.nc")')
    ap.add_argument("-min-points", type=int, dest="min_points", default=None)
    ap.add_argument("-sample-limit", type=int, dest="sample_limit", default=None)
    ap.add_argument("-time", type=str, dest="time_str", default=None,
                    help="End time (e.g., 2025-07-01T00:00) for half-open [start, end) window")
    ap.add_argument("-ndays", type=int, dest="ndays", default=None,
                    help="Number of days back from -time to include (default 0)")
    ap.add_argument("-time-mode", type=str, dest="time_mode", default=None,
                    choices=["auto", "netcdf", "filename", "off"])
    ap.add_argument("-dry-run", action="store_true")
    args = ap.parse_args()

    cfg = {}
    if args.settings_file:
        cfg_path = Path(args.settings_file)
        if not cfg_path.exists():
            print(f"[ERROR] settings_file not found: {cfg_path}")
            sys.exit(2)
        print(f"[INFO] Loading settings from: {cfg_path}")
        with open(cfg_path, "r") as fh:
            cfg = json.load(fh)

    merged = merge_settings(cfg, {
        "shapefile": args.shapefile,
        "src": args.src,
        "dst": args.dst,
        "pattern": args.pattern,
        "min_points": args.min_points,
        "sample_limit": args.sample_limit,
        "time_mode": args.time_mode,
        "ndays": args.ndays,
        "dry_run": args.dry_run
    })

    # Defaults
    merged.setdefault("pattern", "*.nc")
    merged.setdefault("min_points", 1)
    merged.setdefault("sample_limit", None)
    merged.setdefault("time_mode", "auto")
    merged.setdefault("ndays", 0)
    # NEW: optional report_dir setting (can be abs or relative to dst)
    merged.setdefault("report_dir", None)

    # Required
    for key in ("shapefile", "src", "dst"):
        if key not in merged or not merged[key]:
            print(f"[ERROR] Missing required setting: {key}")
            sys.exit(2)

    shapefile = Path(merged["shapefile"])
    src_root = Path(merged["src"]).resolve()
    dst_root = Path(merged["dst"]).resolve()
    pattern = merged["pattern"]
    min_points = int(merged["min_points"])
    sample_limit = merged["sample_limit"]
    time_mode = (merged.get("time_mode") or "auto").lower()
    dry_run = bool(merged.get("dry_run", False))
    ndays = int(merged.get("ndays", 0) or 0)

    start_date = end_date = None
    if args.time_str:
        try:
            # end_date is *exclusive* in the [start, end) convention
            end_date = pd.to_datetime(args.time_str)
            start_date = end_date - pd.Timedelta(days=ndays)
            print(f"[INFO] Scanning only from {start_date} (inclusive) to {end_date} (exclusive)")
        except Exception:
            print(f"[ERROR] Invalid --time value: {args.time_str}")
            sys.exit(2)
    else:
        print("[INFO] No time filter provided — scanning all files.")

    dst_root.mkdir(parents=True, exist_ok=True)
    prepared = load_domain_geometry(shapefile, target_crs="EPSG:4326")

    files = []
    if start_date is not None and end_date is not None:
        print("[INFO] Restricting search to matching date folders...")
        # Walk only day folders that overlap [start_date, end_date] in *calendar* sense
        for yyyy in range(start_date.year, end_date.year + 1):
            year_path = src_root / f"{yyyy:04d}"
            if not year_path.exists():
                continue
            for mm in range(1, 13):
                month_path = year_path / f"{mm:02d}"
                if not month_path.exists():
                    continue
                for dd in range(1, 32):
                    try:
                        dt = pd.Timestamp(f"{yyyy:04d}-{mm:02d}-{dd:02d}")
                    except ValueError:
                        continue
                    if start_date.normalize() <= dt.normalize() <= end_date.normalize():
                        day_path = month_path / f"{dd:02d}"
                        if day_path.exists():
                            print(f"[DEBUG] Scanning folder: {day_path}")
                            files.extend(day_path.rglob(pattern))
    else:
        print(f"[INFO] Scanning all files under {src_root}")
        files = sorted(src_root.rglob(pattern))

    if not files:
        print("[INFO] No files matched pattern.")
        sys.exit(0)

    print(f"[INFO] Found {len(files)} file(s). Starting checks...")
    kept, skipped = [], []

    for f in files:
        if not file_in_range(f, start_date, end_date, mode=time_mode):
            print(f"[TIME-SKIP] {f.relative_to(src_root)}")
            skipped.append(f)
            continue
        if file_is_in_domain(f, prepared, min_points_inside=min_points, sample_limit=sample_limit):
            kept.append(f)
            if not dry_run:
                rel = f.relative_to(src_root)
                out_path = (dst_root / rel).resolve()
                out_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, out_path)
                print(f"[COPY] {f} -> {out_path}")
            print(f"[KEEP] {f.relative_to(src_root)}")
        else:
            skipped.append(f)
            print(f"[SKIP] {f.relative_to(src_root)}")

    # === Report output ===
    # choose a timestamp: use end_date (if provided) else current time
    ts_for_name = (pd.to_datetime(end_date) if end_date is not None else pd.Timestamp.now()).strftime("%Y%m%d%H%M")
    report_filename = f"filter_report_{ts_for_name}.csv"

    # resolve report_dir:
    # - if JSON "report_dir" provided and absolute, use it
    # - if JSON "report_dir" provided and relative, resolve relative to dst
    # - otherwise, default to dst
    report_dir = merged.get("report_dir")
    if report_dir:
        rpt_dir = Path(report_dir)
        if not rpt_dir.is_absolute():
            rpt_dir = (dst_root / rpt_dir).resolve()
    else:
        rpt_dir = dst_root
    rpt_dir.mkdir(parents=True, exist_ok=True)
    report = rpt_dir / report_filename

    with open(report, "w", newline="") as fh:
        import csv
        w = csv.writer(fh)
        w.writerow(["relative_path", "kept"])
        for f in kept:
            w.writerow([str(f.relative_to(src_root)), True])
        for f in skipped:
            w.writerow([str(f.relative_to(src_root)), False])

    print(f"\n[INFO] Done. Kept {len(kept)} / {len(files)} files.")
    print(f"[INFO] Report saved to: {report}")


if __name__ == "__main__":
    main()


