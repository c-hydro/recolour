#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract validation statistics from NetCDF files by GPI and save results
to JSON files named:

    {cell}_{gpi}.json

The cell is inferred from the NetCDF file name, for example:

    1394_2026.nc -> cell 1394

Author: Fabio Delogu

Run
---
python extract_stats.py \
    -settings settings_extract_stats.json

Settings JSON example
---------------------
{
    "datasets": [
        {
            "file": "1394_2026.nc",
            "gpis": [2128135, 2128136]
        },
        {
            "file": "1396_2026.nc",
            "gpis": [2130001, 2130002]
        }
    ],

    "out_folder": "results_json",

    "gpi_var": "gpi",
    "cell_var": "cell"
}
"""

import os
import json
import argparse

import numpy as np
import xarray as xr


def convert_numpy(obj):
    """Convert numpy objects to JSON-compatible Python objects."""

    if isinstance(obj, np.integer):
        return int(obj)

    if isinstance(obj, np.floating):
        if np.isnan(obj):
            return None
        return float(obj)

    if isinstance(obj, np.ndarray):
        return obj.tolist()

    if isinstance(obj, bytes):
        return obj.decode("utf-8")

    return obj


def load_json(file_path):
    """Load a JSON file."""

    with open(file_path, "r", encoding="utf-8") as file_handle:
        return json.load(file_handle)


def save_json(file_path, data):
    """Save dictionary to JSON file."""

    with open(file_path, "w", encoding="utf-8") as file_handle:
        json.dump(
            data,
            file_handle,
            indent=4,
            default=convert_numpy
        )


import re


def infer_cell_from_file(file_path):
    """
    Infer cell from file name and check file existence.

    Supported examples
    ------------------
    1394.nc
    1394_2026.nc
    cell_1394.nc
    metrics_cell_1394_2026.nc
    """

    if not os.path.exists(file_path):
        raise RuntimeError(
            f'File "{file_path}" does not exist'
        )

    file_name = os.path.basename(file_path)

    # extract all integer groups from filename
    numbers = re.findall(r'\d+', file_name)

    if len(numbers) == 0:
        raise RuntimeError(
            f'Cannot infer cell from file name "{file_name}"'
        )

    # first integer is assumed to be the cell
    cell = int(numbers[0])

    print(f' ===> File "{file_name}" -> inferred cell: {cell}')

    return cell

def extract_gpi_statistics(
        file_path,
        gpis,
        cell=None,
        gpi_var="gpi",
        cell_var="cell",
        decode_times=False):
    """
    Extract all statistics for selected GPIs from a NetCDF file.
    """

    if np.isscalar(gpis):
        gpis = [int(gpis)]
    else:
        gpis = [int(gpi) for gpi in gpis]

    results = {}

    with xr.open_dataset(file_path, decode_times=decode_times) as ds:

        if gpi_var not in ds.variables:
            raise RuntimeError(
                f'Variable "{gpi_var}" not found in file "{file_path}"'
            )

        file_gpis = np.asarray(ds[gpi_var].values).squeeze()

        if file_gpis.ndim != 1:
            raise RuntimeError(
                f'Variable "{gpi_var}" must be one-dimensional after squeeze. '
                f"Found shape: {file_gpis.shape}"
            )

        file_gpis_int = file_gpis.astype(int)
        available_gpis = set(file_gpis_int.tolist())

        missing_gpis = [gpi for gpi in gpis if gpi not in available_gpis]

        if missing_gpis:
            print(f" ===> WARNING: GPIs not found in file: {missing_gpis}")

        for gpi in gpis:

            idx = np.where(file_gpis_int == int(gpi))[0]

            if len(idx) == 0:
                continue

            idx = int(idx[0])

            gpi_data = {
                "file": os.path.basename(file_path),
                "cell": int(cell) if cell is not None else None,
                "gpi": int(gpi)
            }

            for var_name in ds.variables:

                values = np.asarray(ds[var_name].values).squeeze()

                try:
                    if values.ndim == 0:
                        gpi_data[var_name] = convert_numpy(values.item())

                    elif values.shape[0] == file_gpis.shape[0]:
                        value = values[idx]

                        if np.ndim(value) == 0:
                            gpi_data[var_name] = convert_numpy(value.item())
                        else:
                            gpi_data[var_name] = convert_numpy(value)

                except Exception as exc:
                    print(
                        f' ===> WARNING: Skip variable "{var_name}" '
                        f"for GPI {gpi}. Reason: {exc}"
                    )

            if cell is not None:
                gpi_data["requested_cell"] = int(cell)

                if cell_var in gpi_data:
                    file_cell = int(gpi_data[cell_var])

                    if file_cell != int(cell):
                        raise RuntimeError(
                            f"GPI {gpi} belongs to cell {file_cell}, "
                            f"not requested/inferred cell {cell}"
                        )
                else:
                    gpi_data["cell_check"] = (
                        f'Variable "{cell_var}" not available in file'
                    )

            results[int(gpi)] = gpi_data

    return results


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract all statistics from one or more NetCDF files for GPIs "
            "defined in a settings JSON file."
        )
    )

    parser.add_argument(
        "-settings",
        dest="settings_file",
        type=str,
        required=True,
        help="Settings JSON file"
    )

    args = parser.parse_args()

    settings = load_json(args.settings_file)

    datasets = settings["datasets"]

    out_folder = settings.get("out_folder", "results_json")
    gpi_var = settings.get("gpi_var", "gpi")
    cell_var = settings.get("cell_var", "cell")

    os.makedirs(out_folder, exist_ok=True)

    for dataset_info in datasets:

        file_path = dataset_info["file"]
        gpis = dataset_info["gpis"]

        cell = infer_cell_from_file(file_path)

        print(f" ===> Process file: {file_path}")
        print(f" ===> Inferred cell: {cell}")

        results = extract_gpi_statistics(
            file_path=file_path,
            gpis=gpis,
            cell=cell,
            gpi_var=gpi_var,
            cell_var=cell_var
        )

        if not results:
            print(f" ===> No results extracted for file {file_path}")
            continue

        for gpi, gpi_data in results.items():

            out_name = f"{cell}_{gpi}.json"
            out_path = os.path.join(out_folder, out_name)

            save_json(out_path, gpi_data)

            print(f" ===> Saved: {out_path}")


if __name__ == "__main__":
    main()