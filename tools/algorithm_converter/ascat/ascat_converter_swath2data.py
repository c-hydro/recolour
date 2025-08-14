
# libraries
import argparse
import os
import logging
import pandas as pd
from pathlib import Path
import geopandas as gpd

import numpy as np

from lib_utils_generic import setup_logging, setup_config
from lib_utils_grid import resample_points2grid, filter_grid
from lib_utils_plot import plot_grid, plot_points, plot_grid_and_points
from lib_utils_io import convert_meters_to_degrees

from lib_io_shapefile import shapefile_to_mask
from lib_io_bufr import decode_bufr
from lib_io_nc import decode_netcdf
from lib_io_tiff import write_tiff
from lib_io_csv import write_csv

# method to process SWATH BUFR/NETCDF files and convert to GeoTIFF
def process_swath_to_tiff(config_path, date_str):

    # load configuration
    config = setup_config(config_path)
    info_cfg = config["info"]
    input_cfg = config["input"]
    output_grids_cfg = config["output_grids"]
    output_points_cfg = config["output_points"]
    input_aux_cfg = config["input_auxiliary"]
    output_aux_cfg = config["output_auxiliary"]
    vars_cfg = config["variables"]
    flags_cfg = config["flags"]
    parameters_cfg = config["parameters"]
    log_cfg = config["log"]

    # product information
    prd_name = info_cfg["product_name"]
    prd_format, prd_ext = info_cfg["product_format"], info_cfg['product_extension']
    prd_tag = ' ::: '.join([prd_name, prd_format])

    # log start information
    logger_obj = setup_logging(log_cfg, date_str)
    # info start of processing
    logger_obj.info(f" --> {prd_tag} processing for {date_str} ... ")

    # organize mask and grid data
    (mask_values,
     mask_transform, mask_bbox,
     mask_lons, mask_lats, land_values, count_values, percentage_values) = shapefile_to_mask(
        shapefile_path_in=input_aux_cfg["shapefile"], shapefile_path_out=output_aux_cfg["shapefile"],
        buffer_km_fine=parameters_cfg["grid_buffer_km_fine"],
        resolution_km_fine=parameters_cfg["grid_spacing_km_fine"],
        resolution_km_coarse=parameters_cfg["grid_spacing_km_coarse"],
        extend_area_flag=parameters_cfg["grid_extend_area_flag"],
        extend_area_buffer_km_land=parameters_cfg["grid_extend_area_km_land"],
        extend_area_buffer_km_sea=parameters_cfg["grid_extend_area_km_sea"],
        update=flags_cfg["update_auxiliary"],
        plot=parameters_cfg["plot_debug"])

    # iterate through time windows
    year, month, day = date_str[:4], date_str[4:6], date_str[6:]
    for window in input_cfg["time_windows"]:

        # info start of processing for this time window
        logger_obj.info(f" ---> Time Window " + str(window) + " ... ")

        # get file parameters
        start_hour, end_hour = window["start_hour"], window["end_hour"]
        label = window["label"]

        # compose file path source
        file_path_src = Path(input_cfg["base_path"]) / year / month / day

        # compose file path destination points
        folder_name_out_pts = Path(output_points_cfg["base_path"]) / year / month / day
        folder_name_out_pts.mkdir(parents=True, exist_ok=True)
        file_name_out_pts = output_points_cfg["filename_pattern"].format(date=date_str, label=label)
        # compose file path destination grids
        folder_name_out_grids = Path(output_grids_cfg["base_path"]) / year / month / day
        folder_name_out_grids.mkdir(parents=True, exist_ok=True)
        file_name_out_grids = output_points_cfg["filename_pattern"].format(date=date_str, label=label)

        # update flags for output
        if flags_cfg["update_output"]:
            if os.path.exists(file_name_out_pts):
                os.remove(file_name_out_pts)
            if os.path.exists(file_name_out_grids):
                os.remove(file_name_out_grids)

        # check if output files already exist
        if not os.path.exists(file_name_out_pts) or not os.path.exists(file_name_out_grids):

            # info start of processing for this time window
            logger_obj.info(
                f" ----> Search {prd_tag} files in: {file_path_src} for window {label} ({start_hour}-{end_hour}) ... ")

            # collect NETCDF files
            list_files = []
            for h in range(start_hour, end_hour):
                hour_str = f"{h:02d}"
                list_files.extend((file_path_src / hour_str).glob("*.nc"))

            # check if any BUFR files were found
            if not list_files:
                # info end of processing for this time window (no files found)
                logger_obj.warning(f" ===> No {prd_tag} files for {label} on {date_str}")

                logger_obj.info(
                    f" ----> Search {prd_tag} files in: {file_path_src} for window {label} ({start_hour}-{end_hour}) ... SKIPPED")

                # info end of processing for this time window
                logger_obj.info(f" ---> Time Window " + str(window) + " ... SKIPPED. Input files are not available.")

                continue

            # log number of NETCDF files found
            logger_obj.info(f" ::: Found {len(list_files)} {prd_tag} files for {label}")

            # read and decode NETCDF files
            all_data = {
                "longitude": [], "latitude": [],
                "ssm": [], "flag_processing": [], "flag_corrections": [], "quality": [],
                'jd': [], 'time': []}
            for file in list_files:
                logging.info(f" ::: Decode file: {file}")

                if file.suffix == prd_ext:
                    logging.info(f" ::: Expected file format: {prd_ext}. Proceeding with decoding.")
                else:
                    logging.error(f" ===> Unsupported file format: {file.suffix}. Expected {prd_ext}. Skipping file.")
                    raise ValueError(f"Unsupported file format: {file.suffix}. Expected {prd_ext}.")

                if file.suffix == ".buf":
                    data = decode_bufr(file)
                elif file.suffix == ".nc":
                    data = decode_netcdf(file)
                else:
                    logging.error(f" ===> Unsupported file format: {file.suffix}. Only .buf/.nc files are supported.")
                    continue
                for k in all_data:
                    all_data[k].extend(data[k])

            values_dict = {
                "ssm": np.array(all_data["ssm"]),
                "flag_processing": np.array(all_data["flag_processing"]),
                "flag_corrections": np.array(all_data["flag_corrections"]),
                "quality": np.array(all_data["quality"]),
                "jd": np.array(all_data["jd"]),
                "time": pd.DatetimeIndex(all_data["time"])
            }

            # info end of processing for this time window
            logger_obj.info(
                f" ----> Search {prd_tag} files in: {file_path_src} for window {label} ({start_hour}-{end_hour}) ... DONE")

            # info start of processing for defining geographical variables
            logger_obj.info(f" ----> Organize geographical variable(s) ... ")

            # get point coordinates
            points_lons = np.array(all_data["longitude"])
            points_lats = np.array(all_data["latitude"])

            # define buffer (in degrees)
            # example ===> buffer_deg = 0.5  # for example, 0.5° latitude/longitude
            buffer_deg_lon, buffer_deg_lat = convert_meters_to_degrees(
                mask_lons, mask_lats, buffer_km=parameters_cfg['points_buffer_km'])

            # extract bounding box (extent) from grid_lons and grid_lats
            min_lon, max_lon = np.min(mask_lons), np.max(mask_lons)
            min_lat, max_lat = np.min(mask_lats), np.max(mask_lats)
            # expand bounding box
            buffered_min_lon = min_lon - buffer_deg_lon
            buffered_max_lon = max_lon + buffer_deg_lon
            buffered_min_lat = min_lat - buffer_deg_lat
            buffered_max_lat = max_lat + buffer_deg_lat

            # create mask with buffer zone
            mask = (
                    (points_lons >= buffered_min_lon) & (points_lons <= buffered_max_lon) &
                    (points_lats >= buffered_min_lat) & (points_lats <= buffered_max_lat)
            )
            # apply mask to points
            selected_lons = points_lons[mask]
            selected_lats = points_lats[mask]

            # info end of processing for defining geographical variables
            logger_obj.info(f" ----> Organize geographical variable(s) ... DONE")

            # iterate through variables and resample to grid
            grid_bands, points_bands = {}, {}

            # info start of processing for defining datasets variable(s)
            logger_obj.info(f" ----> Organize datasets variable(s) ... ")

            # add common grid information
            grid_bands['land'] = land_values
            grid_bands['count'] = count_values
            grid_bands['percentage'] = percentage_values
            for var_name, var_cfg in vars_cfg.items():

                # info start of processing for each variable
                logger_obj.info(f" -----> Variable " + var_name + " ... ")

                # get variable values
                points_values = values_dict[var_name]
                var_min, var_max, var_no_data = var_cfg["min"], var_cfg["max"], var_cfg["no_data"]

                # Apply the mask
                try:
                    selected_values = points_values[mask]
                except BaseException as base_exp:
                    logger_obj.error(f"Error applying mask for variable {var_name}: {base_exp}")
                    selected_values = np.array([])

                # check if selected values are ndarray
                if isinstance(selected_values, np.ndarray):

                    # compute grid resampled and filled
                    grid_resampled, grid_filled = resample_points2grid(
                        selected_lons, selected_lats, selected_values,
                        mask_lons, mask_lats, mask_values,
                        fill_value=var_no_data, max_value=var_max, min_value=var_min,
                        radius_of_influence=parameters_cfg['resample_radius_of_influence_km'],
                        neighbours_min=parameters_cfg['resample_neighbours_min'],
                        neighbours_n=parameters_cfg['resample_neighbours_n'],
                        mask=False,
                        plot=parameters_cfg['plot_debug'])

                    # compute grid filtered
                    grid_filtered = filter_grid(
                        grid_filled, percentage_values,
                        threshold=parameters_cfg['filter_threshold'])

                    # plot if enabled
                    if parameters_cfg['plot_output']:
                        plot_grid_and_points(
                            mask_lons, mask_lats, grid_filtered,
                            selected_lons, selected_lats, selected_values)

                        plot_grid_and_points(
                            mask_lons, mask_lats, grid_filled,
                            selected_lons, selected_lats, selected_values)

                        plot_grid_and_points(
                            mask_lons, mask_lats, grid_resampled,
                            selected_lons, selected_lats, selected_values)

                    # save in common grid bands
                    grid_bands[f"{var_name}_resampled"] = grid_resampled
                    grid_bands[f"{var_name}_filled"] = grid_filled
                    grid_bands[f"{var_name}_filtered"] = grid_filtered

                # save in common points bands
                if 'longitude' not in list(points_bands.keys()):
                    points_bands['longitude'] = selected_lons
                if 'latitude' not in list(points_bands.keys()):
                    points_bands['latitude'] = selected_lats
                points_bands[var_name] = selected_values

                # info start of processing for each variable
                logger_obj.info(f" -----> Variable " + var_name + " ... DONE")

            # info end of processing for defining datasets variable(s)
            logger_obj.info(f" ----> Organize datasets variable(s) ... DONE")

            # info start of dumping for defining datasets variable(s)
            logger_obj.info(f" ----> Dump datasets variable(s) ... ")

            # info start of saving CSV
            logger_obj.info(f" -----> Save CSV ... ")
            folder_name_out = Path(output_points_cfg["base_path"]) / year / month / day
            folder_name_out.mkdir(parents=True, exist_ok=True)
            file_name_out = output_points_cfg["filename_pattern"].format(date=date_str, label=label)
            path_name_out = os.path.join(folder_name_out, file_name_out)

            # save grids datasets as CSV
            logger_obj.info(f" ::: Saved CSV file: {path_name_out} ")
            # dump points bands to CSV
            write_csv(path_name_out, points_bands, time_column='time')
            # info end of saving CSV
            logger_obj.info(f" -----> Save CSV ... DONE")

            # info start of saving GeoTIFF
            logger_obj.info(f" -----> Save GeoTiff ... ")

            # compose file path destination
            folder_name_out = Path(output_grids_cfg["base_path"]) / year / month / day
            folder_name_out.mkdir(parents=True, exist_ok=True)
            file_name_out = output_grids_cfg["filename_pattern"].format(date=date_str, label=label)
            path_name_out = os.path.join(folder_name_out, file_name_out)

            # save grids datasets as GeoTIFF
            logger_obj.info(f" ::: Saved GeoTiff file: {path_name_out} ")
            # dump grid bands to GeoTIFF
            write_tiff(path_name_out, mask_lons, mask_lats, grid_bands)
            # info end of aving GeoTIFF
            logger_obj.info(f" -----> Save GeoTiff ... DONE")

            # info start of dumping for defining datasets variable(s)
            logger_obj.info(f" ----> Dump datasets variable(s) ... DONE ")

            # info end of processing for this time window
            logger_obj.info(f" ---> Time Window " + str(window) + " ... DONE")

        else:

            # info end of processing for this time window
            logger_obj.info(f" ---> Time Window " + str(window) + " ... SKIPPED. Output files already exist.")

    # info end of processing
    logger_obj.info(f" --> {prd_tag} processing for {date_str} ... DONE")

# main method to run the script
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description="Convert ASCAT SWATH BUFR/NETCDF to GeoTIFF")
    parser.add_argument("config", help="Path to config.json")
    parser.add_argument("date", help="Date in YYYYMMDD format")
    args = parser.parse_args()

    # method to process SWATH BUFR/NETCDF files and convert to GeoTIFF
    process_swath_to_tiff(args.config, args.date)

