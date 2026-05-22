"""
Library Features:

Name:           lib_utils_analysis_points
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""


# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd
import xarray as xr

from scipy.spatial import cKDTree

from config_info import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to check 1d arrays
def check_1d_locations_var(dataset: xr.Dataset, var_name: (str, list) = None):

    # check variable name
    if var_name is None:
        raise KeyError('Variable name must be defined. Found NoneType.')

    # allow single string or list/tuple of candidate names
    if isinstance(var_name, str):
        var_names = [var_name]
    else:
        var_names = list(var_name)

    # find first existing variable/coord
    selected_name = None
    for name in var_names:
        if name in dataset.variables or name in dataset.coords:
            selected_name = name
            break
    if selected_name is None:
        raise KeyError(f"None of the variables {var_names} found")

    # get data array
    data_array = dataset[selected_name]

    # check dimensions
    valid_dims = [("locations",), ("obs",)]
    if data_array.dims not in valid_dims:
        raise ValueError(
            f"Variable '{selected_name}' must have dims "
            f"{valid_dims}, got {data_array.dims}"
        )

    return selected_name
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to collect porosity collection by gpi
def collect_porosity_to_dataframe(parameters, file_list):

    value_var = parameters.get("variable_name", "porosity_0_5cm_mean")

    porosity_min = float(parameters.get("min_value", 0.05))
    porosity_max = float(parameters.get("max_value", 0.90))

    frames = []

    for file_path in file_list:

        logger.info(f" ----> Read porosity file: {file_path}")

        if not os.path.exists(file_path):
            raise RuntimeError(f"Porosity file not found: {file_path}")

        file_name = os.path.basename(file_path)

        try:
            cell = int(
                file_name
                .replace("soilgrids_cell_", "")
                .replace(".nc", "")
            )
        except Exception:
            raise RuntimeError(
                f'Cannot extract cell identifier from porosity filename "{file_name}"'
            )

        logger.info(f" :::: Porosity cell: {cell}")
        logger.info(f" ::: GET POROSITY POINTS ... ")

        with xr.open_dataset(file_path) as dataset:

            var_name_geo_x = check_1d_locations_var(dataset, ["lon", "longitude", "x", "X"])
            logger.info(f" :::: Variable GeoX: {var_name_geo_x}")
            var_name_geo_y = check_1d_locations_var(dataset, ["lat", "latitude", "y", "Y"])
            logger.info(f" :::: Variable GeoY: {var_name_geo_y}")
            var_name_loc = check_1d_locations_var(dataset, ["gpi", "location_id"])
            logger.info(f" :::: Variable GPI/Locations: {var_name_loc}")
            var_name_data = check_1d_locations_var(dataset, value_var)
            logger.info(f" :::: Variable Porosity: {var_name_data}")

            lon = np.asarray(dataset[var_name_geo_x].values, dtype=np.float64)
            lat = np.asarray(dataset[var_name_geo_y].values, dtype=np.float64)
            gpi = np.asarray(dataset[var_name_loc].values, dtype=np.int64)
            porosity = np.asarray(dataset[var_name_data].values, dtype=np.float64)

            if not (gpi.shape == porosity.shape):
                raise ValueError(
                    f'{file_path}: shape mismatch among gpi/{value_var}: '
                    f'{gpi.shape}, {porosity.shape}'
                )

            valid_mask = (
                np.isfinite(gpi) &
                np.isfinite(porosity) &
                np.isfinite(lon) &
                np.isfinite(lat) &
                (porosity > porosity_min) &
                (porosity < porosity_max)
            )

            if valid_mask.any():
                frame = pd.DataFrame(
                    {
                        "cell": np.full(np.sum(valid_mask), cell, dtype=np.int32),
                        "gpi": gpi[valid_mask],
                        "lon": lon[valid_mask],
                        "lat": lat[valid_mask],
                        "porosity": porosity[valid_mask],
                        "porosity_file": file_name
                    }
                )
                frames.append(frame)

        logger.info(f" ::: GET POROSITY POINTS ... DONE")

    # manage merge porosity points
    logger.info(f" ::: MERGE POROSITY POINTS ... ")

    if not frames:
        raise RuntimeError(
            "Porosity conversion is enabled, but no valid porosity points were found"
        )

    merge_frames = pd.concat(frames, ignore_index=True)
    merge_frames = merge_frames.drop_duplicates(subset=["cell", "gpi"], keep="last")

    logger.info(f" ::: MERGE POROSITY POINTS ... DONE")

    return merge_frames
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# point collection and dedup
def collect_points_to_dataframe(source_settings, file_list, grid, max_distance_km=25):

    value_var = source_settings.get("value_col", "surface_soil_moisture")

    logger.info(" ::: GET DATASETS POINTS ... ")
    frames = []
    for file_path in file_list:
        logger.info(f" :::: Read datasets file: {file_path}")

        with xr.open_dataset(file_path) as dataset:

            var_name_geo_x = check_1d_locations_var(dataset, ["lon", "longitude", "x", "X"])
            var_name_geo_y = check_1d_locations_var(dataset, ["lat", "latitude", "y", "Y"])
            var_name_time = check_1d_locations_var(dataset, "time")
            var_name_data = check_1d_locations_var(dataset, value_var)

            lon = np.asarray(dataset[var_name_geo_x].values, dtype=np.float64)
            lat = np.asarray(dataset[var_name_geo_y].values, dtype=np.float64)
            obs_time = pd.to_datetime(dataset[var_name_time].values)
            obs_value = np.asarray(dataset[var_name_data].values, dtype=np.float64)

            if not (lon.shape == lat.shape == obs_time.shape == obs_value.shape):
                raise ValueError(
                    f'{file_path}: shape mismatch among lon/lat/time/{value_var}: '
                    f'{lon.shape}, {lat.shape}, {obs_time.shape}, {obs_value.shape}'
                )

            valid_mask = (
                np.isfinite(lon) &
                np.isfinite(lat) &
                np.isfinite(obs_value)
            )

            if valid_mask.any():

                valid_lon = lon[valid_mask]
                valid_lat = lat[valid_mask]

                nearest_gpi, distance = grid.find_nearest_gpi(
                    valid_lon,
                    valid_lat,
                    max_dist=max_distance_km * 1000.0
                )

                nearest_gpi = np.asarray(nearest_gpi, dtype=np.int64)
                distance = np.asarray(distance, dtype=np.float64)

                valid_grid = nearest_gpi >= 0

                if valid_grid.any():

                    valid_gpi = nearest_gpi[valid_grid]
                    valid_cell = grid.gpi2cell(valid_gpi).astype(np.int32)

                    frame = pd.DataFrame(
                        {
                            "cell": valid_cell,
                            "gpi": valid_gpi,
                            "lon": valid_lon[valid_grid],
                            "lat": valid_lat[valid_grid],
                            "time": obs_time[valid_mask][valid_grid],
                            value_var: obs_value[valid_mask][valid_grid],
                            "distance_grid_m": distance[valid_grid],
                            "source_file": os.path.basename(file_path)
                        }
                    )

                    frames.append(frame)

    logger.info(" ::: GET DATASETS POINTS ... DONE")

    logger.info(" ::: MERGE DATASETS POINTS ... ")
    if not frames:
        merge_frames = pd.DataFrame(
            columns=[
                "cell", "gpi", "lon", "lat", "time",
                value_var, "distance_grid_m", "source_file"
            ]
        )
    else:
        merge_frames = pd.concat(frames, ignore_index=True)

    logger.info(" ::: MERGE DATASETS POINTS ... DONE")

    return merge_frames
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# remove duplicates from points
def deduplicate_latest_points(df, value_var, reference_time=None, reference_flag=False):

    info = {
        "initial_rows": 0,
        "reference_active": bool(reference_flag),
        "reference_time": reference_time,
        "filtered_by_reference_time": 0,
        "duplicates_by_value": 0,
        "duplicates_by_gpi": 0,
        "final_rows": 0,
    }

    if df.empty:
        logger.warning(" ===> DataFrame is empty")
        return df.copy(), info

    info["initial_rows"] = len(df)

    logger.info(" ::: DEDUPLICATE POINTS ... ")
    logger.info(f" :::: Initial rows: {info['initial_rows']}")

    df = df.copy()

    # filter by reference time if active
    if reference_flag:
        if reference_time is None:
            logger.warning(" ===> Reference flag is active but reference_time is None")
        else:
            rows_before_filter = len(df)

            df = df.loc[df["time"] <= reference_time].copy()

            info["filtered_by_reference_time"] = rows_before_filter - len(df)

            logger.info(f" :::: Reference time active: {reference_time}")
            logger.info(
                f" :::: Rows filtered by reference time: "
                f"{info['filtered_by_reference_time']}"
            )

    if df.empty:
        logger.warning(" ===> DataFrame is empty after reference time filtering")
        info["final_rows"] = 0
        return df.copy(), info

    # sort by gpi and time
    df = df.sort_values(["gpi", "time"]).copy()

    # duplicates based on ["gpi", value_var]
    dup_mask_value = df.duplicated(subset=["gpi", value_var], keep="last")
    dup_points_value = df.loc[dup_mask_value, ["gpi", "time", value_var]]

    info["duplicates_by_value"] = len(dup_points_value)
    logger.info(f" :::: Duplicates by value: {info['duplicates_by_value']}")

    # remove duplicates by value, keeping latest
    df = df.drop_duplicates(subset=["gpi", value_var], keep="last")

    # duplicates based on ["gpi"]
    dup_mask_gpi = df.duplicated(subset=["gpi"], keep="last")
    dup_points_gpi = df.loc[dup_mask_gpi, ["gpi", "time", value_var]]

    info["duplicates_by_gpi"] = len(dup_points_gpi)
    logger.info(f" :::: Duplicates by gpi/locations: {info['duplicates_by_gpi']}")

    # keep latest point per gpi
    df = (
        df.sort_values(["gpi", "time"])
          .drop_duplicates(subset=["gpi"], keep="last")
    )

    # final sorting
    df = df.sort_values("gpi").reset_index(drop=True)

    info["final_rows"] = len(df)
    logger.info(f" :::: Final rows: {info['final_rows']}")

    logger.info(" ::: DEDUPLICATE POINTS ... DONE")

    return df, info
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# convert point values from SMAP VWC to ASCAT-like SSM using porosity
def convert_points_vwc_to_ssm(points_df, porosity_df, value_var, parameters):

    fill_value = float(parameters.get("fill_value", -9999.0))
    clip_min = float(parameters.get("clip_min", 0.0))
    clip_max = float(parameters.get("clip_max", 100.0))

    radius_list = parameters.get("radius_km", [5, 10, 15, 20, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
    lat_col = parameters.get("lat_col", "lat")
    lon_col = parameters.get("lon_col", "lon")

    if "gpi" not in points_df.columns:
        raise RuntimeError('Column "gpi" is needed in points dataframe')
    if "gpi" not in porosity_df.columns:
        raise RuntimeError('Column "gpi" is needed in porosity dataframe')
    if "porosity" not in porosity_df.columns:
        raise RuntimeError('Column "porosity" is needed in porosity dataframe')
    if value_var not in points_df.columns:
        raise RuntimeError(f'Column "{value_var}" is needed in points dataframe')

    for col in [lat_col, lon_col]:
        if col not in points_df.columns:
            raise RuntimeError(f'Column "{col}" is needed in points dataframe')
        if col not in porosity_df.columns:
            raise RuntimeError(f'Column "{col}" is needed in porosity dataframe')

    points_df = points_df.copy()
    porosity_df = porosity_df.copy()

    assign_info = {
        "n_points_initial": int(points_df.shape[0]),
        "n_assigned_by_gpi": 0,
        "n_assigned_by_radius": 0,
        "n_not_assigned": 0,
        "n_valid_ssm": 0,
        "n_removed_invalid_ssm": 0,
        "radius": {}
    }

    points_df["assign_method"] = None
    points_df["assign_id"] = -1
    points_df["assign_params"] = np.nan

    # First assignment: by gpi
    points_df = points_df.merge(
        porosity_df[["gpi", "porosity"]],
        on="gpi",
        how="left"
    )

    by_gpi = np.isfinite(points_df["porosity"])

    points_df.loc[by_gpi, "assign_method"] = "by_gpi"
    points_df.loc[by_gpi, "assign_id"] = 1
    points_df.loc[by_gpi, "assign_params"] = np.nan

    assign_info["n_assigned_by_gpi"] = int(by_gpi.sum())

    # Second assignment: nearest porosity point by radius
    missing = ~by_gpi

    if missing.any():

        valid_porosity = porosity_df[
            np.isfinite(porosity_df["porosity"]) &
            np.isfinite(porosity_df[lat_col]) &
            np.isfinite(porosity_df[lon_col])
        ].copy()

        if not valid_porosity.empty:

            lat0 = np.nanmean(points_df[lat_col])
            km_per_deg_lat = 111.32
            km_per_deg_lon = 111.32 * np.cos(np.deg2rad(lat0))

            porosity_xy = np.column_stack([
                valid_porosity[lon_col].values * km_per_deg_lon,
                valid_porosity[lat_col].values * km_per_deg_lat,
            ])

            tree = cKDTree(porosity_xy)

            missing_idx = points_df.index[
                missing &
                np.isfinite(points_df[lat_col]) &
                np.isfinite(points_df[lon_col])
            ]

            points_xy = np.column_stack([
                points_df.loc[missing_idx, lon_col].values * km_per_deg_lon,
                points_df.loc[missing_idx, lat_col].values * km_per_deg_lat,
            ])

            for radius_km in radius_list:

                still_missing = points_df.loc[missing_idx, "porosity"].isna().values

                if not still_missing.any():
                    break

                idx_to_check = missing_idx[still_missing]
                xy_to_check = points_xy[still_missing]

                dist, nn = tree.query(
                    xy_to_check,
                    distance_upper_bound=float(radius_km)
                )

                found = np.isfinite(dist)

                n_found_radius = int(found.sum())
                assign_info["radius"][float(radius_km)] = n_found_radius

                if n_found_radius > 0:
                    found_idx = idx_to_check[found]
                    found_nn = nn[found]

                    points_df.loc[found_idx, "porosity"] = (
                        valid_porosity.iloc[found_nn]["porosity"].values
                    )

                    points_df.loc[found_idx, "assign_method"] = "by_radius"
                    points_df.loc[found_idx, "assign_id"] = 2
                    points_df.loc[found_idx, "assign_params"] = float(radius_km)

    assign_info["n_assigned_by_radius"] = int(
        (points_df["assign_method"] == "by_radius").sum()
    )

    assign_info["n_not_assigned"] = int(
        points_df["porosity"].isna().sum()
    )

    valid = (
        np.isfinite(points_df[value_var]) &
        np.isfinite(points_df["porosity"]) &
        (points_df["porosity"] > 0)
    )

    points_df["surface_soil_moisture"] = fill_value

    points_df.loc[valid, "surface_soil_moisture"] = (
        100.0 *
        points_df.loc[valid, value_var] /
        points_df.loc[valid, "porosity"]
    )

    points_df.loc[valid, "surface_soil_moisture"] = np.clip(
        points_df.loc[valid, "surface_soil_moisture"],
        clip_min,
        clip_max
    )

    valid_ssm = (
        np.isfinite(points_df["surface_soil_moisture"]) &
        (points_df["surface_soil_moisture"] != fill_value)
    )

    assign_info["n_valid_ssm"] = int(valid_ssm.sum())
    assign_info["n_removed_invalid_ssm"] = int((~valid_ssm).sum())

    points_df = points_df[valid_ssm]

    return points_df, "surface_soil_moisture", assign_info
# ----------------------------------------------------------------------------------------------------------------------
