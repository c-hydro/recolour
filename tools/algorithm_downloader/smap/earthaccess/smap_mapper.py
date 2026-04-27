#!/usr/bin/env python3
"""
Quick SMAP L2 soil moisture plotter with domain lookup.

Features
--------
- Read SMAP L2 HDF5 soil moisture, latitude, longitude
- Build bbox from:
    * predefined domain names (italy, europe, global, ...)
    * country names using Natural Earth through cartopy
- Plot zoomed scatter map
- Optionally save figure

Example
-------
python smap_plot_domain.py \
    --file SMAP_L2_SM_P_E_59818_A_20260413T153804_R19240_001.h5 \
    --domain italy \
    --save smap_italy.png

python smap_plot_domain.py \
    --file SMAP_L2_SM_P_E_59818_A_20260413T153804_R19240_001.h5 \
    --domain france

python smap_plot_domain.py \
    --file SMAP_L2_SM_P_E_59818_A_20260413T153804_R19240_001.h5 \
    --bbox 6 36 19 47.5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.io import shapereader
    HAS_CARTOPY = True
except Exception:
    HAS_CARTOPY = False


# -------------------------------------------------------------------------------------
# Predefined named domains
DOMAIN_BOXES = {
    "italy": (6.0, 36.0, 19.0, 47.5),
    "europe": (-12.0, 34.0, 35.0, 72.0),
    "mediterranean": (-7.0, 28.0, 40.0, 48.5),
    "global": (-180.0, -90.0, 180.0, 90.0),
}

DEFAULT_PARAMS = {
    # --- IO ---
    "group_name": "Soil_Moisture_Retrieval_Data",
    "output_figure": None,

    # --- Domain ---
    "domain": "italy",
    "bbox": None,
    "pad_deg": 0.5,

    # --- Data filtering ---
    "good_quality_only": False,
    "valid_min": 0.0,
    "valid_max": 1.0,

    # --- Plot ---
    "point_size": 7.0,
    "vmin": 0.0,
    "vmax": 0.6,
    "figsize": (9, 8),
    "cmap": None,  # use matplotlib default

    # --- Behaviour ---
    "show": True,
    "save": None,
}

# -------------------------------------------------------------------------------------
def read_smap_soil_moisture(file_name: str, group_name: str = "Soil_Moisture_Retrieval_Data"):
    """
    Read SMAP soil moisture, lat, lon and quality flag from the HDF5 file.
    """
    with h5py.File(file_name, "r") as f:
        grp = f[group_name]

        soil_moisture = grp["soil_moisture"][:]
        latitude = grp["latitude"][:]
        longitude = grp["longitude"][:]

        retrieval_qual_flag = None
        if "retrieval_qual_flag" in grp:
            retrieval_qual_flag = grp["retrieval_qual_flag"][:]

    return {
        "soil_moisture": soil_moisture,
        "latitude": latitude,
        "longitude": longitude,
        "retrieval_qual_flag": retrieval_qual_flag,
    }


def clean_smap_data(soil_moisture,
                    latitude, longitude,
                    retrieval_qual_flag=None, good_quality_only=False,
                    valid_min=0.0, valid_max=1.0):
    """
    Flatten arrays and mask invalid values.
    """
    sm = np.asarray(soil_moisture, dtype=float)
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)

    # Common physical mask for volumetric soil moisture
    sm[(sm < 0.0) | (sm > 1.0)] = np.nan

    if good_quality_only and retrieval_qual_flag is not None:
        qf = np.asarray(retrieval_qual_flag)
        # Common quick filter: keep only flag == 0
        sm[qf != 0] = np.nan

    sm = sm.ravel()
    lat = lat.ravel()
    lon = lon.ravel()

    sm[(sm < valid_min) | (sm > valid_max)] = np.nan

    valid = np.isfinite(sm) & np.isfinite(lat) & np.isfinite(lon)

    return lon[valid], lat[valid], sm[valid]


def get_bbox_from_country_name(country_name: str, pad_deg: float = 0.5):
    """
    Get bounding box from Natural Earth country polygons using cartopy.
    Returns (west, south, east, north).
    """
    if not HAS_CARTOPY:
        raise RuntimeError("Cartopy is required for country-name lookup")

    shp = shapereader.natural_earth(
        resolution="110m",
        category="cultural",
        name="admin_0_countries",
    )

    reader = shapereader.Reader(shp)

    target = country_name.strip().lower()

    for record in reader.records():
        attrs = record.attributes

        candidates = [
            attrs.get("NAME"),
            attrs.get("NAME_LONG"),
            attrs.get("ADMIN"),
            attrs.get("SOVEREIGNT"),
            attrs.get("FORMAL_EN"),
        ]
        candidates = [c for c in candidates if c]

        if any(str(c).strip().lower() == target for c in candidates):
            minx, miny, maxx, maxy = record.geometry.bounds
            return (minx - pad_deg, miny - pad_deg, maxx + pad_deg, maxy + pad_deg)

    raise ValueError(f"Country/domain '{country_name}' not found")


def get_bbox(area_name: str | None = None, bbox=None, pad_deg: float = 0.5):
    """
    Resolve bbox from:
    - explicit bbox
    - predefined domains
    - country name lookup
    """
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError("bbox must be (west, south, east, north)")
        return tuple(float(v) for v in bbox)

    if area_name is None:
        raise ValueError("Either area_name or bbox must be provided")

    key = area_name.strip().lower()

    if key in DOMAIN_BOXES:
        return DOMAIN_BOXES[key]

    return get_bbox_from_country_name(area_name, pad_deg=pad_deg)


def subset_to_bbox(lon, lat, sm, bbox):
    """
    Keep only points inside bbox.
    """
    west, south, east, north = bbox
    mask = (lon >= west) & (lon <= east) & (lat >= south) & (lat <= north)
    return lon[mask], lat[mask], sm[mask]


def plot_soil_moisture(
    lon,
    lat,
    sm,
    bbox,
    title="SMAP Soil Moisture",
    point_size=DEFAULT_PARAMS["point_size"],
    vmin=DEFAULT_PARAMS["vmin"],
    vmax=DEFAULT_PARAMS["vmax"],
    figsize=DEFAULT_PARAMS["figsize"],
    save_path=None,
    show=True,
):
    """
    Plot soil moisture over a bbox.
    """
    west, south, east, north = bbox

    if HAS_CARTOPY:
        fig = plt.figure(figsize=figsize)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([west, east, south, north], crs=ccrs.PlateCarree())

        ax.coastlines(resolution="10m", linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, alpha=0.2)
        ax.add_feature(cfeature.OCEAN, alpha=0.15)
        ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4)

        sc = ax.scatter(
            lon,
            lat,
            c=sm,
            s=point_size,
            transform=ccrs.PlateCarree(),
            vmin=0.0,
            vmax=0.6,
        )
        plt.colorbar(sc, ax=ax, shrink=0.8, label="Soil moisture [m3/m3]")
        ax.set_title(title)

    else:
        fig, ax = plt.subplots(figsize=(9, 8))
        sc = ax.scatter(lon, lat, c=sm, s=point_size, vmin=vmin, vmax=vmax)
        ax.set_xlim(west, east)
        ax.set_ylim(south, north)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        plt.colorbar(sc, ax=ax, shrink=0.8, label="Soil moisture [m3/m3]")

    plt.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"Saved figure: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


# -------------------------------------------------------------------------------------
def build_title(file_name: str, area_name: str | None, good_quality_only: bool):
    base = Path(file_name).name
    quality_txt = " | QF=0 only" if good_quality_only else ""
    area_txt = f" | area={area_name}" if area_name else ""
    return f"SMAP soil_moisture{area_txt}{quality_txt}\n{base}"


def parse_args():
    parser = argparse.ArgumentParser(description="Plot SMAP soil moisture over a named area or bbox")

    parser.add_argument("--file", required=True)

    parser.add_argument("--group", default=DEFAULT_PARAMS["group_name"])
    parser.add_argument("--domain", default=DEFAULT_PARAMS["domain"])
    parser.add_argument("--bbox", nargs=4, type=float)

    parser.add_argument("--pad", type=float, default=DEFAULT_PARAMS["pad_deg"])

    parser.add_argument(
        "--good-quality-only",
        action="store_true",
        default=DEFAULT_PARAMS["good_quality_only"],
    )

    parser.add_argument("--point-size", type=float, default=DEFAULT_PARAMS["point_size"])

    parser.add_argument("--vmin", type=float, default=DEFAULT_PARAMS["vmin"])
    parser.add_argument("--vmax", type=float, default=DEFAULT_PARAMS["vmax"])

    parser.add_argument("--save", default=DEFAULT_PARAMS["save"])

    parser.add_argument("--no-show", action="store_true")

    return parser.parse_args()


def main():
    args = parse_args()

    data = read_smap_soil_moisture(args.file, group_name=args.group)

    lon, lat, sm = clean_smap_data(
        soil_moisture=data["soil_moisture"],
        latitude=data["latitude"],
        longitude=data["longitude"],
        retrieval_qual_flag=data["retrieval_qual_flag"],
        good_quality_only=args.good_quality_only,
        valid_min=args.vmin,
        valid_max=args.vmax,
    )

    bbox = get_bbox(area_name=args.domain if args.bbox is None else None, bbox=args.bbox, pad_deg=args.pad)

    lon_sub, lat_sub, sm_sub = subset_to_bbox(lon, lat, sm, bbox)

    print(f"Resolved bbox: {bbox}")
    print(f"Points in bbox: {sm_sub.size}")

    if sm_sub.size == 0:
        print("No valid soil moisture points found inside the selected domain")
        return

    title = build_title(args.file, args.domain if args.bbox is None else "custom_bbox", args.good_quality_only)

    save_path = "/home/fabio/Desktop/recolour/exec/smap/smap_italy.png"

    plot_soil_moisture(
        lon=lon_sub,
        lat=lat_sub,
        sm=sm_sub,
        bbox=bbox,
        title=title,
        point_size=args.point_size,
        save_path=save_path,
        #save_path=args.save,
        show=not args.no_show,
    )

    print()

if __name__ == "__main__":
    main()