"""
RESHUFFLE SOILGRIDS
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import sys
import shutil
import logging
import argparse
import numpy as np

from copy import deepcopy

from lib_grid_hmc import load_grid, subgrid4bbox
from lib_interface_soilgrids import soilgrids_ds
from lib_img2cell_soilgrids import Img2Cell
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to reset folder
def reset_folder(folder_name):

    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)

    os.makedirs(folder_name, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert string to boolean
def str2bool(val):

    if val in ["True", "true", "t", "T", "1"]:
        return True
    else:
        return False
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to reshuffle data
def reshuffle(
        product,
        dataset_grid_root,
        dataset_ts_root,
        grid_path,
        parameters,
        file_name_tmpl_grid="{variable}_{depth}_{statistic}.tif",
        file_name_tmpl_ts="{cell_n}.nc",
        depth="0-5cm",
        statistic="mean",
        variable_mapping=None,
        cell_format_ts="%04d",
        reset_ts=True,
        input_grid=None,
        target_grid=None,
        bbox=None,
        img_buffer=1):

    # ------------------------------------------------------------------------------------------------------------------
    # Section 0 - variable(s) initialization
    if reset_ts:
        reset_folder(dataset_ts_root)

    os.makedirs(dataset_ts_root, exist_ok=True)

    if variable_mapping is None:
        variable_mapping = {}
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 1 - datasets type and format
    logging.info(' -----> 1) Define datasets type and format ... ')

    dataset_grid_driver = soilgrids_ds(
        data_path=dataset_grid_root,
        file_name_tmpl=file_name_tmpl_grid,
        depth=depth,
        statistic=statistic,
        parameter=parameters,
        grid_path=grid_path,
        subgrid=input_grid,
        array_1D=True
    )

    logging.info(' -----> 1) Define datasets type and format ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 2 - datasets metadata
    logging.info(' -----> 2) Define datasets metadata ... ')

    global_attr = {
        "product": product,
        "source": "SoilGrids",
        "depth": depth,
        "statistic": statistic,
        "static": "true"
    }

    cell_attributes = {}

    for variable_name in parameters:

        var_out_name = variable_mapping.get(
            variable_name,
            variable_name
        )

        cell_attributes[var_out_name] = {
            "source_variable": variable_name,
            "depth": depth,
            "statistic": statistic,
            "units": "-"
        }

    if input_grid is None:
        input_grid = load_grid(
            land_points=True,
            grid_path=grid_path,
            bbox=None
        )

    if bbox is not None:

        min_lon = bbox[0]
        min_lat = bbox[1]
        max_lon = bbox[2]
        max_lat = bbox[3]

        target_grid = subgrid4bbox(
            input_grid,
            min_lon=min_lon,
            min_lat=min_lat,
            max_lon=max_lon,
            max_lat=max_lat
        )

    else:
        target_grid = deepcopy(input_grid)

    logging.info(' -----> 2) Define datasets metadata ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 3 - datasets conversion from grid to cell
    logging.info(' -----> 3) Convert datasets from grid to cell ... ')

    drv_reshuffle = Img2Cell(
        input_dataset=dataset_grid_driver,
        outputpath=dataset_ts_root,
        input_grid=input_grid,
        target_grid=target_grid,
        variable_rename=variable_mapping,
        cellsize_lat=5.0,
        cellsize_lon=5.0,
        filename_templ=file_name_tmpl_ts,
        cell_templ=cell_format_ts,
        gridname="grid.nc",
        global_attr=global_attr,
        cell_attributes=cell_attributes,
        cell_dtypes=np.float32,
        zlib=True
    )

    drv_reshuffle.calc(parameters)

    logging.info(' -----> 3) Convert datasets from grid to cell ... DONE')
    # ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to parse argument(s)
def parse_args(args):

    parser = argparse.ArgumentParser(
        description="Convert static SoilGrids data to cell netCDF format."
    )

    parser.add_argument(
        "product",
        default="soilgrids",
        help="Product name"
    )

    parser.add_argument(
        "flags",
        nargs=2,
        type=str2bool,
        help="Flags for reset static grid and output cell datasets."
    )

    parser.add_argument(
        "bbox",
        type=str,
        default="",
        help="Boundary box 'lon_min, lat_min, lon_max, lat_max'"
    )

    parser.add_argument(
        "image_buffer",
        type=str,
        default=1,
        help="Image buffer parameter"
    )

    parser.add_argument(
        "dataset_grid_root",
        help="Root folder where static SoilGrids maps are stored"
    )

    parser.add_argument(
        "dataset_ts_root",
        help="Root folder where cell netCDF files should be stored"
    )

    parser.add_argument(
        "dataset_stack_root",
        help="Root folder for stack files; unused for static SoilGrids"
    )

    parser.add_argument(
        "grid_path",
        type=str,
        default=os.path.dirname(__file__),
        help="Path where the grid reference is stored"
    )

    parser.add_argument(
        "templates_src",
        metavar="templates_src",
        type=str,
        nargs=1,
        help="Source filename template"
    )

    parser.add_argument(
        "templates_dst",
        metavar="templates_dst",
        type=str,
        nargs=1,
        help="Destination filename template"
    )

    parser.add_argument(
        "depth",
        type=str,
        default="0-5cm",
        help="SoilGrids depth, for example 0-5cm"
    )

    parser.add_argument(
        "statistic",
        type=str,
        default="mean",
        help="SoilGrids statistic, for example mean"
    )

    parser.add_argument(
        "parameters",
        metavar="parameters",
        nargs="+",
        help="SoilGrids variables to reshuffle"
    )

    parser.add_argument(
        "--land_points",
        type=str2bool,
        default="True",
        help="Set True to convert only land points"
    )

    parser.add_argument(
        "--bbox",
        type=float,
        default=None,
        nargs=4,
        help="min_lon min_lat max_lon max_lat"
    )

    args = parser.parse_args(args)

    return args
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main(args):

    product = args["product"]["name"]

    bbox = args["product"]["bbox"]

    reset_ts = args["flags"]["reset_dynamic"]

    dataset_grid_root = args["data"]["path_grid"]
    dataset_ts_root = args["data"]["path_ts"]

    grid_path = os.path.join(
        args["grid"]["folder_name"],
        args["grid"]["file_name"]
    )

    parameters = args["variables"]

    file_name_tmpl_grid = args["template"]["file"]["file_name_source"]
    file_name_tmpl_ts = args["template"]["file"]["file_name_destination"]

    depth = args["soilgrids"].get("depth", "0-5cm")
    statistic = args["soilgrids"].get("statistic", "mean")

    variable_mapping = args["soilgrids"].get("variable_mapping", {})

    cell_format_ts = args["template"]["dataset"].get(
        "cell_n",
        "%04d"
    )

    # method to load grid
    logging.info(' ----> Get reference grid ... ')

    input_grid = load_grid(
        land_points=True,
        grid_path=grid_path,
        bbox=None
    )

    logging.info(' ----> Get reference grid ... DONE')

    # method to reshuffle data from grid to cell format
    logging.info(' ----> Run reshuffle algorithm to convert grid to cell ... ')

    reshuffle(
        product,
        dataset_grid_root,
        dataset_ts_root,
        grid_path,
        parameters,
        file_name_tmpl_grid=file_name_tmpl_grid,
        file_name_tmpl_ts=file_name_tmpl_ts,
        depth=depth,
        statistic=statistic,
        variable_mapping=variable_mapping,
        cell_format_ts=cell_format_ts,
        reset_ts=reset_ts,
        input_grid=input_grid,
        target_grid=None,
        bbox=bbox,
        img_buffer=int(args["product"].get("image_buffer", 1))
    )

    logging.info(' ----> Run reshuffle algorithm to convert grid to cell ... DONE')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to run from command-line
def run():
    main(sys.argv[1:])
# ----------------------------------------------------------------------------------------------------------------------