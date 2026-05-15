#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr


def inspect_gpi(reference_grid, gpi_target, max_print=10):

    print("\n============================================================")
    print(f" ---> Inspect GPI: {gpi_target}")
    print(f" ---> Grid file:   {reference_grid}")
    print("------------------------------------------------------------")

    with xr.open_dataset(reference_grid) as ds:

        if "gpi" not in ds:
            raise KeyError("Variable 'gpi' not found in grid file")

        gpis = ds["gpi"].values.astype(np.int64)

        print(f" ---> GPI range: {int(gpis.min())} - {int(gpis.max())}")

        mask = gpis == int(gpi_target)

        if mask.any():
            idx = np.where(mask)[0][0]

            print(" ---> STATUS: FOUND")
            print(f" ---> Index: {idx}")

            if "cell" in ds:
                print(f" ---> Cell: {int(ds['cell'].values[idx])}")

            if "lat" in ds:
                print(f" ---> Lat: {float(ds['lat'].values[idx])}")
            if "lon" in ds:
                print(f" ---> Lon: {float(ds['lon'].values[idx])}")

            if "latitude" in ds:
                print(f" ---> Lat: {float(ds['latitude'].values[idx])}")
            if "longitude" in ds:
                print(f" ---> Lon: {float(ds['longitude'].values[idx])}")

            return True

        print(" ---> STATUS: NOT FOUND")

        distances = np.abs(gpis - int(gpi_target))
        idx_sorted = np.argsort(distances)[:max_print]

        print(f" ---> Nearest GPIs: {gpis[idx_sorted].tolist()}")

        if "cell" in ds:
            print(f" ---> Nearest Cells: {ds['cell'].values[idx_sorted].tolist()}")

        if "latitude" in ds:
            print(f" ---> Nearest Lat: {ds['latitude'].values[idx_sorted].tolist()}")
        if "longitude" in ds:
            print(f" ---> Nearest Lon: {ds['longitude'].values[idx_sorted].tolist()}")

        if gpi_target < gpis.min():
            print(" ---> Position: BELOW grid range")
        elif gpi_target > gpis.max():
            print(" ---> Position: ABOVE grid range")
        else:
            print(" ---> Position: INSIDE range but missing")

        return False


# -------------------------------------------------------------------------------------
# OPTIONAL: helper for multiple GPIs
def inspect_multiple_gpis(reference_grid, gpi_list):

    print("\n==================== MULTI GPI CHECK ====================")

    with xr.open_dataset(reference_grid) as ds:
        gpis = set(ds["gpi"].values.astype(np.int64).tolist())

    found = []
    missing = []

    for gpi in gpi_list:
        if gpi in gpis:
            found.append(gpi)
        else:
            missing.append(gpi)

    print(f" ---> Total:   {len(gpi_list)}")
    print(f" ---> Found:   {len(found)}")
    print(f" ---> Missing: {len(missing)}")

    if missing:
        print(f" ---> Missing sample: {missing[:20]}")

    return found, missing


# -------------------------------------------------------------------------------------
def main():

    grid_file = (
        "/home/fabio/Desktop/recolour/dset/ascat_auxiliary/reference/ascat/grid/"
        "TUW_WARP5_grid_info_2_3.nc"
    )

    gpi_target = 3264395
    #gpi_target = [  3264395, 3264396,1000000, 2000000]
    inspect_gpi(grid_file, gpi_target)


if __name__ == "__main__":
    main()