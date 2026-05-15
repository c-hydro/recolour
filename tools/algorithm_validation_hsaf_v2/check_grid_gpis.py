#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Map reference TUW GPIs to nearest GLDAS and CCI GPIs.

Run
---
python map_gpis.py \
    -settings settings_gpi_map.json \
    -out_json gpi_mapping.json

Settings JSON example
---------------------
{
    "cells": {
        "1394": [2128135, 2128136],
        "1395": [2129000]
    },

    "ref_grid": "TUW_WARP5_grid_info_2_2.nc",
    "k1_grid": "grid_gldas.nc",
    "k2_grid": "grid_cci.nc",

    "ref_name": "TUW",
    "k1_name": "GLDAS",
    "k2_name": "CCI",

    "ref_gpi_var": "gpi",
    "ref_lat_var": "lat",
    "ref_lon_var": "lon",
    "ref_cell_var": "cell",

    "k1_gpi_var": "gpi",
    "k1_lat_var": "lat",
    "k1_lon_var": "lon",
    "k1_cell_var": "cell",

    "k2_gpi_var": "gpi",
    "k2_lat_var": "lat",
    "k2_lon_var": "lon",
    "k2_cell_var": "cell"
}

Output JSON example
-------------------
{
    "1394": {
        "2128135": {
            "TUW": {
                "gpi": 2128135,
                "lat": 40.85997,
                "lon": 14.825693,
                "cell": 1394,
                "found": true
            },
            "GLDAS": {
                "gpi": 581099,
                "lat": 40.875,
                "lon": 14.875,
                "cell": 403,
                "distance_km": 4.51
            },
            "CCI": {
                "gpi": 753899,
                "lat": 40.875,
                "lon": 14.875,
                "cell": 1789,
                "distance_km": 4.51
            }
        }
    }
}
"""

import json
import argparse

import numpy as np
import xarray as xr


def to_json_value(value):

    if isinstance(value, np.integer):
        return int(value)

    if isinstance(value, np.floating):
        if np.isnan(value):
            return None
        return float(value)

    if isinstance(value, np.ndarray):
        return value.tolist()

    return value


def haversine_distance_km(lon1, lat1, lon2, lat2):
    """Compute great-circle distance in km."""

    radius_earth_km = 6371.0

    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    )

    c = 2.0 * np.arcsin(np.sqrt(a))

    return radius_earth_km * c


def read_grid(
        grid_file,
        gpi_var="gpi",
        lat_var="lat",
        lon_var="lon",
        cell_var="cell"):
    """Read grid information from NetCDF."""

    with xr.open_dataset(grid_file) as ds:

        for var_name in [gpi_var, lat_var, lon_var]:

            if var_name not in ds.variables:
                raise RuntimeError(
                    f'Variable "{var_name}" not found in "{grid_file}"'
                )

        grid = {
            "gpi": np.asarray(ds[gpi_var].values).squeeze(),
            "lat": np.asarray(ds[lat_var].values).squeeze(),
            "lon": np.asarray(ds[lon_var].values).squeeze(),
        }

        if cell_var in ds.variables:
            grid["cell"] = np.asarray(ds[cell_var].values).squeeze()
        else:
            grid["cell"] = None

    return grid


def get_point_by_gpi(grid, gpi):
    """Extract point info from GPI."""

    idx = np.where(grid["gpi"].astype(int) == int(gpi))[0]

    if len(idx) == 0:
        return None

    idx = int(idx[0])

    point = {
        "gpi": int(grid["gpi"][idx]),
        "lat": float(grid["lat"][idx]),
        "lon": float(grid["lon"][idx]),
    }

    if grid["cell"] is not None:
        point["cell"] = int(grid["cell"][idx])

    return point


def find_nearest_point(grid, lon, lat):
    """Find nearest point in target grid."""

    distances = haversine_distance_km(
        lon,
        lat,
        grid["lon"],
        grid["lat"]
    )

    idx = int(np.nanargmin(distances))

    point = {
        "gpi": int(grid["gpi"][idx]),
        "lat": float(grid["lat"][idx]),
        "lon": float(grid["lon"][idx]),
        "distance_km": float(distances[idx]),
    }

    if grid["cell"] is not None:
        point["cell"] = int(grid["cell"][idx])

    return point


def map_ref_gpis_to_target_grids(
        cells,
        ref_grid,
        k1_grid,
        k2_grid,
        ref_name="TUW",
        k1_name="GLDAS",
        k2_name="CCI"):
    """
    Map TUW GPIs to nearest GLDAS and CCI GPIs.
    """

    output = {}

    for cell, gpis in cells.items():

        output[str(cell)] = {}

        for gpi in gpis:

            ref_point = get_point_by_gpi(ref_grid, gpi)

            if ref_point is None:

                output[str(cell)][str(gpi)] = {
                    ref_name: {
                        "gpi": int(gpi),
                        "found": False
                    },
                    k1_name: None,
                    k2_name: None
                }

                continue

            ref_point["found"] = True

            k1_point = find_nearest_point(
                k1_grid,
                lon=ref_point["lon"],
                lat=ref_point["lat"]
            )

            k2_point = find_nearest_point(
                k2_grid,
                lon=ref_point["lon"],
                lat=ref_point["lat"]
            )

            output[str(cell)][str(gpi)] = {
                ref_name: ref_point,
                k1_name: k1_point,
                k2_name: k2_point
            }

    return output


def load_json(file_path):

    with open(file_path, "r", encoding="utf-8") as file_handle:
        return json.load(file_handle)


def save_json(file_path, data):

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(
            data,
            file_handle,
            indent=4,
            default=to_json_value
        )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Map TUW reference GPIs to nearest GLDAS and CCI GPIs."
        )
    )

    parser.add_argument(
        "-settings",
        dest="settings_file",
        type=str,
        required=True,
        help="Settings JSON file"
    )

    parser.add_argument(
        "-out_json",
        dest="out_json",
        type=str,
        default="gpi_mapping.json",
        help="Output JSON file"
    )

    args = parser.parse_args()

    settings = load_json(args.settings_file)

    cells = settings["cells"]

    ref_grid = read_grid(
        settings["ref_grid"],
        gpi_var=settings.get("ref_gpi_var", "gpi"),
        lat_var=settings.get("ref_lat_var", "lat"),
        lon_var=settings.get("ref_lon_var", "lon"),
        cell_var=settings.get("ref_cell_var", "cell")
    )

    k1_grid = read_grid(
        settings["k1_grid"],
        gpi_var=settings.get("k1_gpi_var", "gpi"),
        lat_var=settings.get("k1_lat_var", "lat"),
        lon_var=settings.get("k1_lon_var", "lon"),
        cell_var=settings.get("k1_cell_var", "cell")
    )

    k2_grid = read_grid(
        settings["k2_grid"],
        gpi_var=settings.get("k2_gpi_var", "gpi"),
        lat_var=settings.get("k2_lat_var", "lat"),
        lon_var=settings.get("k2_lon_var", "lon"),
        cell_var=settings.get("k2_cell_var", "cell")
    )

    mapping = map_ref_gpis_to_target_grids(
        cells=cells,
        ref_grid=ref_grid,
        k1_grid=k1_grid,
        k2_grid=k2_grid,
        ref_name=settings.get("ref_name", "TUW"),
        k1_name=settings.get("k1_name", "GLDAS"),
        k2_name=settings.get("k2_name", "CCI")
    )

    save_json(args.out_json, mapping)

    print(f' ===> Saved mapping file: {args.out_json}')


if __name__ == "__main__":
    main()