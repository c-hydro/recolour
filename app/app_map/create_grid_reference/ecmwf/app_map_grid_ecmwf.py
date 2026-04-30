#!/usr/bin/env python3
"""
Create a compact NetCDF grid file with:

    longitude(lon)  -> 1D
    latitude(lat)   -> 1D
    mask(lat, lon)  -> 2D

The mask is created from var40:
    mask = 1 where var40 is finite/valid
    mask = 0 where var40 is missing/NaN

The mask is also shifted along longitude to align original 0..360 data
with output longitude -180..180.

Example:
    python create_grid_mask_1d.py \
        -i h26_2026042600_R01.nc \
        -o grid_h26_1d.nc
"""

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


def pick_horizontal_dims(ds, y_dim=None, x_dim=None):
    """Detect latitude/y and longitude/x dimensions."""

    if y_dim is None:
        for candidate in ["lat", "latitude", "south_north", "y"]:
            if candidate in ds.dims:
                y_dim = candidate
                break

    if x_dim is None:
        for candidate in ["lon", "longitude", "west_east", "x"]:
            if candidate in ds.dims:
                x_dim = candidate
                break

    if y_dim is None or x_dim is None:
        dims = list(ds.dims)
        if len(dims) < 2:
            raise ValueError("Cannot detect horizontal dimensions.")
        y_dim = dims[-2]
        x_dim = dims[-1]

    return y_dim, x_dim


def squeeze_to_2d_grid(da, y_dim, x_dim):
    """
    Convert a variable to a 2D horizontal grid.

    Drops singleton dimensions.
    If other non-horizontal dimensions remain, selects index 0.
    """

    da = da.squeeze(drop=True)

    for dim in list(da.dims):
        if dim not in (y_dim, x_dim):
            da = da.isel({dim: 0})

    if set(da.dims) != {y_dim, x_dim}:
        raise ValueError(
            f"Variable cannot be reduced to ({y_dim}, {x_dim}). "
            f"Current dims: {da.dims}"
        )

    return da.transpose(y_dim, x_dim)


def create_1d_lat_lon(
    nx,
    ny,
    lon_min=-180.0,
    lon_max=180.0,
    lat_min=-90.0,
    lat_max=90.0,
):
    """
    Create 1D longitude and latitude arrays.

    endpoint=False avoids duplicating -180 and +180 on a cyclic grid.
    """

    longitude = np.linspace(
        lon_min,
        lon_max,
        nx,
        endpoint=False,
        dtype=np.float32,
    )

    latitude = np.linspace(
        lat_max,
        lat_min,
        ny,
        dtype=np.float32,
    )

    return latitude, longitude


def save_1d_grid_and_mask(
    file_in,
    file_out="grid_h26_1d.nc",
    mask_variable="var40",
    y_dim=None,
    x_dim=None,
    lon_min=-180.0,
    lon_max=180.0,
    lat_min=-90.0,
    lat_max=90.0,
    apply_lon_correction=True,
):
    """
    Save longitude(lon), latitude(lat), and mask(lat, lon).

    Parameters
    ----------
    file_in : str
        Input NetCDF file.
    file_out : str
        Output NetCDF file.
    mask_variable : str
        Variable used to create the mask. Default is var40.
    apply_lon_correction : bool
        If True, roll mask by nx//2 cells along longitude.
        This aligns source data on 0..360 longitude with output -180..180.
    """

    ds = xr.open_dataset(file_in)

    if mask_variable not in ds:
        raise KeyError(f"Variable '{mask_variable}' not found in input file.")

    y_dim, x_dim = pick_horizontal_dims(ds, y_dim=y_dim, x_dim=x_dim)

    ny = ds.sizes[y_dim]
    nx = ds.sizes[x_dim]

    latitude, longitude = create_1d_lat_lon(
        nx=nx,
        ny=ny,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
    )

    source = squeeze_to_2d_grid(ds[mask_variable], y_dim=y_dim, x_dim=x_dim)

    mask = xr.where(np.isfinite(source), 1, 0).astype("int8")
    mask = mask.transpose(y_dim, x_dim)

    if apply_lon_correction:
        lon_roll_cells = nx // 2
        mask = mask.roll({x_dim: lon_roll_cells}, roll_coords=False)
        correction_status = "active"
    else:
        lon_roll_cells = 0
        correction_status = "inactive"

    out_ds = xr.Dataset(
        data_vars={
            "mask": ((y_dim, x_dim), mask.values),
        },
        coords={
            "latitude": (y_dim, latitude),
            "longitude": (x_dim, longitude),
        },
        attrs={
            "title": "1D latitude longitude grid with aligned mask",
            "source_file": Path(file_in).name,
            "mask_source_variable": mask_variable,
            "source_longitude_convention": "0_to_360",
            "output_longitude_convention": "-180_to_180",
            "longitude_correction_status": correction_status,
            "longitude_roll_cells": int(lon_roll_cells),
            "geospatial_lon_min": float(lon_min),
            "geospatial_lon_max": float(lon_max),
            "geospatial_lat_min": float(lat_min),
            "geospatial_lat_max": float(lat_max),
        },
    )

    out_ds = out_ds.rename({
        y_dim: "lat",
        x_dim: "lon",
    })

    out_ds["latitude"].attrs = {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
        "axis": "Y",
    }

    out_ds["longitude"].attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "axis": "X",
    }

    out_ds["mask"].attrs = {
        "long_name": f"mask derived from {mask_variable}",
        "description": (
            f"1 where {mask_variable} is finite/valid; "
            f"0 where {mask_variable} is missing or NaN."
        ),
        "flag_values": np.array([0, 1], dtype=np.int8),
        "flag_meanings": "invalid valid",
        "source_variable": mask_variable,
        "longitude_correction_status": correction_status,
        "longitude_roll_cells": int(lon_roll_cells),
    }

    encoding = {
        "latitude": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "longitude": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "mask": {"zlib": True, "complevel": 4, "dtype": "int8"},
    }

    out_ds.to_netcdf(file_out, encoding=encoding)

    ds.close()
    out_ds.close()


def main():
    parser = argparse.ArgumentParser(
        description="Save 1D latitude, 1D longitude, and mask derived from var40."
    )

    parser.add_argument("-i", "--input", required=True, help="Input NetCDF file")
    parser.add_argument("-o", "--output", default="grid_h26.nc", help="Output NetCDF file")
    parser.add_argument("--mask-variable", default="var40", help="Variable used to create mask")
    parser.add_argument("--y-dim", default=None, help="Input latitude/y dimension name")
    parser.add_argument("--x-dim", default=None, help="Input longitude/x dimension name")

    parser.add_argument("--lon-min", type=float, default=-180.0)
    parser.add_argument("--lon-max", type=float, default=180.0)
    parser.add_argument("--lat-min", type=float, default=-90.0)
    parser.add_argument("--lat-max", type=float, default=90.0)

    parser.add_argument(
        "--no-lon-correction",
        action="store_true",
        help="Disable longitude rolling correction.",
    )

    args = parser.parse_args()

    save_1d_grid_and_mask(
        file_in=args.input,
        file_out=args.output,
        mask_variable=args.mask_variable,
        y_dim=args.y_dim,
        x_dim=args.x_dim,
        lon_min=args.lon_min,
        lon_max=args.lon_max,
        lat_min=args.lat_min,
        lat_max=args.lat_max,
        apply_lon_correction=not args.no_lon_correction,
    )


if __name__ == "__main__":
    main()