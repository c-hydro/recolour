"""
RESHUFFLE

Module for a command line interface to convert the HMC data into a
time series format using the repurpose package

__date__ = '20230627'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_grid2ts.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20230627 (1.0.0) --> First development
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import sys
import shutil
import logging
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from copy import deepcopy
from pygeogrids import BasicGrid

from lib_utils_cci import read_obj, write_obj
from lib_img2cell_cci import Img2Ts
from lib_interface_cci import cci_ds
from lib_grid_cci import load_grid, cell_grid, subgrid4bbox
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to reset folder
def reset_folder(folder_name):
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name, exist_ok=True)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to find file type based on extension
def get_filetype(grid_path):

    # get folder level(s)
    one_down = os.path.join(grid_path, os.listdir(grid_path)[0])
    if os.path.isdir(one_down):
        level_down = deepcopy(one_down)
    else:
        # check if one down is enough by checking if file or directory
        level_down = os.path.join(one_down, os.listdir(one_down)[0]) \
            if os.path.isdir(one_down) else os.path.dirname(one_down)

    extension = None
    for path, dirs, files in os.walk(level_down):
        if extension is not None:
            break
        for name in files:
            filename, extension = os.path.splitext(name)
            break

    if extension == '.nc' or extension == '.nc4':
        return 'netcdf'
    elif extension == '.tiff' or extension == '.tif':
        return 'tiff'
    elif extension == '.grb' or extension == '.grib':
        return 'grib'
    else:
        logging.error(' ===> File format for grid datasets "' + str(extension) + '" is not supported')
        raise NotImplemented('Case not implemented yet')

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to convert sub path from string to list
def convert_sub_path(sub_path_string='', sep_string='/'):
    sub_path_parts = sub_path_string.split(sep=sep_string)
    sub_path_list = []
    for sub_path_element in sub_path_parts:
        if sub_path_element:
            sub_path_list.append(sub_path_element)
    return sub_path_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compose sub path with time
def compose_sub_path(root_path, time_path, template_path=None):
    tmp_path = ''
    if template_path is not None:
        for s in template_path:
            tmp_path = os.path.join(tmp_path, time_path.strftime(s))
    sub_path = os.path.join(root_path, tmp_path)
    return sub_path
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compose file name
def compose_file_name(file_name_tmpl,
                      file_time, time_placeholder='datetime_source', time_format='%Y%m%d',
                      file_cell=None, cell_placeholder='cell_n'):

    if file_cell is None:
        file_cell = '{cell_n}'

    if not isinstance(file_time, str):
        file_time = file_time.strftime(time_format)

    file_tags = {time_placeholder: file_time, cell_placeholder: file_cell}
    file_name_filled = file_name_tmpl.format(**file_tags)
    return file_name_filled
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create date string
def mkdate(datestring):
    """
    Create date string.

    Parameters
    ----------
    datestring : str
        Date string.

    Returns
    -------
    datestr : datetime
        Date string as datetime.
    """
    if len(datestring) == 10:
        return datetime.strptime(datestring, "%Y-%m-%d")
    if len(datestring) == 16 and 'T' in datestring:
        return datetime.strptime(datestring, "%Y-%m-%dT%H:%M")
    if len(datestring) == 16 and 'T' not in datestring:
        return datetime.strptime(datestring, "%Y-%m-%d %H:%M")
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
def reshuffle(product,
              dataset_grid_root, dataset_ts_root,
              start_date, end_date, run_date,
              grid_path,
              parameters,
              file_name_tmpl_grid="ESACCI-SOILMOISTURE-L3S-SSMV-PASSIVE-{datetime}-fv08.1.nc",
              datetime_format_grid="%Y%m%d000000",
              data_sub_path_grid=None,
              file_name_tmpl_ts=None,
              datetime_format_ts=None,
              cell_format_ts='%04d',
              data_sub_path_ts=None,
              reset_ts=True,
              input_grid=None, target_grid=None, bbox=None,
              img_buffer=50,
              dataset_stack_root=None,
              stack_flag=True,
              stack_file_name_tmpl='{cell_n}.stack',
              stack_reset=True):

    # ------------------------------------------------------------------------------------------------------------------
    # Section 0 - variable(s) initialization (if needed)
    if run_date is None:
        logging.error(' ===> The variable "run_date" is defined by NoneType')
        raise RuntimeError('The variable must be defined by finite element. Check your time settings.')

    # parse the path(s) element
    if data_sub_path_grid is None:
        data_sub_path_grid = ['%Y', '%m', '%d']
    if isinstance(data_sub_path_grid, str):
        data_sub_path_grid = convert_sub_path(data_sub_path_grid)
    if data_sub_path_ts is None:
        data_sub_path_ts = ['%Y', '%m']
    if isinstance(data_sub_path_ts, str):
        data_sub_path_ts = convert_sub_path(data_sub_path_ts)

    if dataset_stack_root is None:
        dataset_stack_root = deepcopy(dataset_ts_root)

    # set ts
    dataset_ts_root = compose_sub_path(dataset_ts_root, run_date, data_sub_path_ts)
    file_name_tmpl_ts = compose_file_name(
        file_name_tmpl_ts, file_time=run_date, time_placeholder='datetime_destination')
    os.makedirs(dataset_ts_root, exist_ok=True)
    if reset_ts:
        reset_folder(dataset_ts_root)
    # set stack
    dataset_stack_root = compose_sub_path(dataset_stack_root, run_date, data_sub_path_ts)
    # set start and end hour
    start_step, end_step = start_date.hour, end_date.hour
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 1 - datasets type and format
    logging.info(' -----> 1) Define datasets type and format ... ')

    if get_filetype(dataset_grid_root) == "netcdf":
        if input_grid is not None:
            logging.info(" ------> Land Grid is fit to CCI grid netCDF data")
            dataset_grid_driver = cci_ds(
                product=product,
                data_path=dataset_grid_root,
                grid_path=grid_path,
                data_sub_path=data_sub_path_grid,
                file_name_tmpl=file_name_tmpl_grid,
                datetime_format=datetime_format_grid,
                parameter=parameters, subgrid=input_grid, array_1D=True
            )
            logging.info(' -----> 1) Define file type and format ... DONE')
        else:
            logging.info(' -----> 1) Define file type and format ... FAILED')
            logging.error(' ===> Grid is not defined')
            raise RuntimeError('Grid must be defined to correctly run the algorithm')

    else:
        logging.info(' -----> 1) Define file type and format ... FAILED')
        logging.error(' ===> File type is not supported by the procedure')
        raise NotImplemented('Case not implemented yet')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 2 - datasets metadata
    logging.info(' -----> 2) Define datasets metadata ... ')

    # define global attributes
    global_attr = {'product': product}

    # iterate over selected time stamps
    data = None
    tstamp_list = dataset_grid_driver.tstamps_for_daterange(start_date, end_date)
    for tstamp in tstamp_list:
        logging.info(' -------> Analyze time "' + tstamp.isoformat() + '" ... ')

        # retrieve data information to initialize datasets
        try:
            data = dataset_grid_driver.read(tstamp)
        except BaseException as base_exp:
            logging.warning(' ===> Datasets are not available (' + str(base_exp) + ')')
            logging.info(' -------> Analyze time "' + tstamp.isoformat() + '" ... SKIPPED')

        if data is not None:
            logging.info(' -------> Analyze time "' + tstamp.isoformat() + '" ... DONE')
            break

    # data is not available for all steps
    if data is None:
        logging.error(' ===> All datasets are not available')
        logging.info(' -----> 2) Define datasets metadata ... FAILED. No data for the day. Exit ...')
        raise SystemExit('Algorithm will exit')

    # define ts attributes
    ts_attributes = {}
    if hasattr(data, 'metadata'):
        ts_attributes = data.metadata

    # get grid information
    if input_grid is None:
        if hasattr(data, 'lon') and hasattr(data, 'lat'):
            input_grid = BasicGrid(data.lon, data.lat)
        else:
            logging.error(' ===> Data obj attributes "lon" and "lat" are not available')
            raise RuntimeError('Attributes "lon" and "lat" needed by the procedure')

    # select grid based on boundary box
    if bbox is not None:
        min_lon, min_lat, max_lon, max_lat = bbox[0], bbox[1], bbox[2], bbox[3]
        target_grid = subgrid4bbox(input_grid, min_lon=min_lon, min_lat=min_lat, max_lon=max_lon, max_lat=max_lat)
    else:
        target_grid = deepcopy(input_grid)

    logging.info(' -----> 2) Define datasets metadata ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Section 3 - datasets conversion from grid to time-series
    logging.info(' -----> 3) Convert datasets from grid to time-series ... ')

    #'''
    # method to initialize reshuffle to convert grid to time-series
    drv_reshuffle = Img2Ts(
        input_dataset=dataset_grid_driver,
        outputpath=dataset_ts_root,
        startdate=start_date,
        enddate=end_date,
        input_grid=input_grid,
        target_grid=target_grid,
        imgbuffer=img_buffer,
        cellsize_lat=5.0,
        cellsize_lon=5.0,
        global_attr=global_attr,
        zlib=True,
        unlim_chunksize=None,
        filename_templ=file_name_tmpl_ts,
        cell_templ=cell_format_ts,
        gridname='grid.nc',
        ts_attributes=ts_attributes,
        stack_history=stack_flag,
        filestack_templ='stack_{:}_{:}.workspace',
        stack_id_templ='%04d',
        stackpath=dataset_stack_root
    )
    # convert grid to ts format
    stack_obj = drv_reshuffle.calc()

    # stack datasets
    logging.info(' ------> Organize stack datasets ... ')
    if stack_obj is not None:

        # iterate over cell object(s)
        for stack_cell, stack_file_list in stack_obj.items():

            # cell string
            string_cell = cell_format_ts % stack_cell
            # info cell start
            logging.info(' -------> Dump cell "' + string_cell + '" ... ')

            # join stacks (over chunks)
            stack_file_collections = None
            for stack_file_path in stack_file_list:
                stack_file_data = read_obj(stack_file_path)
                if stack_file_collections is None:
                    if stack_file_data is not None:
                        stack_file_collections = deepcopy(stack_file_data)
                    else:
                        logging.warning(' ===> Stack file "' + stack_file_path + '" is defined by NoneType')
                else:
                    stack_file_collections = pd.merge(stack_file_collections, stack_file_data,
                                                      left_index=True, right_index=True)
                # reset stacks tmp file(s)
                if stack_reset:
                    if os.path.exists(stack_file_path):
                        os.remove(stack_file_path)

            # save stacks evaluating all the chunks selected by the time steps
            if stack_file_collections is not None:
                stack_index = stack_file_collections.index.values
                stack_values = stack_file_collections.values
                stack_idx_finite_any = np.where(np.any(np.isfinite(stack_values), axis=1))[0]
                stack_idx_nan_all = np.where(np.all(np.isnan(stack_values), axis=1))[0]

                stack_summary_arr = np.zeros(shape=(stack_index.shape[0]))
                stack_summary_arr[:] = np.nan
                stack_summary_arr[stack_idx_finite_any] = 1
                stack_summary_arr[stack_idx_nan_all] = 0

                stack_summary_df = pd.DataFrame(data={'stack': stack_summary_arr}, index=stack_index)

                stack_file_name_summary = stack_file_name_tmpl.format(cell_n=string_cell)
                stack_file_path_summary = os.path.join(dataset_stack_root, stack_file_name_summary)

                if os.path.exists(stack_file_path_summary):
                    os.remove(stack_file_path_summary)

                write_obj(stack_file_path_summary, stack_summary_df)

                # info cell end
                logging.info(' -------> Dump cell "' + string_cell + '" ... DONE')
            else:
                logging.info(' -------> Dump cell "' + string_cell + '" ... SKIPPED. All stacks are defined by NoneType')

        logging.info(' ------> Organize stack datasets ... DONE')
    else:
        logging.info(' ------> Organize stack datasets ... SKIPPED. Flag is not activated')

    logging.info(' -----> 3) Convert datasets from grid to time-series ... DONE')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to parse argument(s)
def parse_args(args):
    """
    Parse command line parameters.

    Parameters
    ----------
    args : list of str
        Command line parameters as list of strings.

    Returns
    -------
    args : argparse.Namespace
        Command line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Convert CCI data to time series format."
    )

    parser.add_argument(
        "product",
        default="product",
        help="Product name",
    )

    parser.add_argument(
        "flags",
        nargs=2,
        type=str2bool,
        help="Flags for reset static grid and output ts."
    )

    parser.add_argument(
        "bbox",
        type=str,
        default='',
        help="Boundary box 'lon_min, lat_min, lon_max, lat_max]' "
    )

    parser.add_argument(
        "dataset_grid_root",
        help="Root of local filesystem where the data is stored",
    )

    parser.add_argument(
        "dataset_ts_root",
        help="Root of local filesystem where the timeseries should be stored"
    )

    parser.add_argument(
        "dataset_stack_root",
        help="Root of local filesystem where the stacks should be stored"
    )

    parser.add_argument(
        "start",
        type=mkdate,
        help="time_start. Either in format 'YYYY-MM-DD', 'YYYY-MM-DDTHH:MM' or 'YYYY-MM-DD HH:MM'"
    )

    parser.add_argument(
        "end",
        type=mkdate,
        help="time_end. Either in format 'YYYY-MM-DD', 'YYYY-MM-DDTHH:MM' or 'YYYY-MM-DD HH:MM'"
    )

    parser.add_argument(
        "run",
        type=mkdate,
        help="time_run. Either in format 'YYYY-MM-DD', 'YYYY-MM-DDTHH:MM' or 'YYYY-MM-DD HH:MM'"
    )

    parser.add_argument(
        "grid_path",
        type=str,
        default=os.path.dirname(__file__),
        help="Path of local filesystem where the grid reference is stored."
    )

    parser.add_argument(
        "templates_src",
        metavar="templates_src",
        type=str,
        nargs=3,
        help="Filename template, datetime template and sub path template"
    )

    parser.add_argument(
        "templates_dst",
        metavar="templates_dst",
        type=str,
        nargs=3,
        help="Filename template, datetime template and sub path template"
    )

    parser.add_argument(
        "parameters",
        metavar="parameters",
        nargs="+",
        help="Parameters to reshuffle into time series format."
    )

    parser.add_argument(
        "--land_points",
        type=str2bool,
        default="True",
        help=(
            "Set True to convert only land points as defined"
            " in the HMC land mask (faster and less/smaller files)"
        ),
    )

    parser.add_argument(
        "--bbox",
        type=float,
        default=None,
        nargs=4,
        help=(
            "min_lon min_lat max_lon max_lat. "
            "Bounding Box (lower left and upper right corner) "
            "of area to reshuffle (WGS84)"
        ),
    )

    parser.add_argument(
        "--imgbuffer",
        type=int,
        default=200,
        help=(
            "How many images to read at once. Bigger "
            "numbers make the conversion faster but "
            "consume more memory."
        ),
    )

    # set defaults that can not be handled by argparse
    args = parser.parse_args(args)

    return args
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main(args):
    """
    Main routine used for command line interface.

    Parameters
    ----------
    args : list of str
        Command line arguments.
    """

    # method to parse arguments
    args = parse_args(args)

    # info start
    logging.info(
        ' ---> Convert CCI datasets from {} to {} into folder {} ... '.format(
            args.start.isoformat(), args.end.isoformat(), args.dataset_ts_root
        )
    )

    # method to load grid
    logging.info(' ----> Get reference grid ... ')
    input_grid = load_grid(
        land_points=args.land_points, grid_path=args.grid_path,
        bbox=None
        # bbox=tuple(args.bbox) if args.bbox is not None else None,
    )
    logging.info(' ----> Get reference grid ... DONE')

    # method to get bbox
    logging.info(' ----> Get reference bbox ... ')
    bbox_list_str, bbox_list_arr = args.bbox.split(','), None
    for bbox_elem in bbox_list_str:
        if (isinstance(bbox_elem, str)) and (bbox_elem != ''):
            if bbox_list_arr is None:
                bbox_list_arr = []
            bbox_list_arr.append(float(bbox_elem))
    logging.info(' ----> Get reference bbox ... DONE')

    # method to reshuffle data from grid to time-series
    logging.info(' ----> Run reshuffle algorithm to convert grid to time-series ... ')

    reshuffle(
        args.product,
        args.dataset_grid_root,
        args.dataset_ts_root,
        args.start,
        args.end,
        args.run,
        args.grid_path,
        args.parameters,
        input_grid=input_grid,
        img_buffer=args.imgbuffer,
        file_name_tmpl_grid=args.templates_src[0],
        datetime_format_grid=args.templates_src[1],
        data_sub_path_grid=args.templates_src[2],
        file_name_tmpl_ts=args.templates_dst[0],
        datetime_format_ts=args.templates_dst[1],
        data_sub_path_ts=args.templates_dst[2],
        bbox=bbox_list_arr,
        target_grid=None,
        reset_ts=args.flags[1],
        dataset_stack_root=args.dataset_stack_root,
        stack_flag=True,
        stack_reset=True
    )

    logging.info(' ----> Run reshuffle algorithm to convert grid to time-series ... DONE')

    # info end
    logging.info(
        ' ---> Convert CCI datasets from {} to {} into folder {} ... DONE'.format(
            args.start.isoformat(), args.end.isoformat(), args.dataset_ts_root
        )
    )
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to run from command-line
def run():
    main(sys.argv[1:])
# ----------------------------------------------------------------------------------------------------------------------
