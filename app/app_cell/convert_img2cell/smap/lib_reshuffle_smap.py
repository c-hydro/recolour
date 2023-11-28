"""
RESHUFFLE

Module for a command line interface to convert the SMAP data into a
time series format using the repurpose package

__date__ = '20230711'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
    'Martina Natali (martina01.natali@edu.unife.it)'
__library__ = 'recolour'

General command line:
python app_grid2ts_smap.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20230711 (1.0.0) --> First development
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# libraries
import os
import sys
import shutil
import logging
import argparse

from copy import deepcopy
from datetime import datetime, timedelta
from pygeogrids import BasicGrid

from lib_img2cell_smap import Img2Ts
from lib_interface_smap import smap_ds
from lib_grid_smap import subgrid4bbox, cell_grid
# -------------------------------------------------------------------------------------


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

    one_down = os.path.join(grid_path, os.listdir(grid_path)[0])
    two_down = os.path.join(one_down, os.listdir(one_down)[0])

    extension = None
    for path, dirs, files in os.walk(two_down):
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


# -------------------------------------------------------------------------------------
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
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to convert string to boolean
def str2bool(val):
    if val in ["True", "true", "t", "T", "1"]:
        return True
    else:
        return False
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to reshuffle data
def reshuffle(product,
              dataset_grid_root, dataset_ts_root,
              start_date, end_date, run_date,
              grid_path,
              parameters,
              file_name_tmpl_grid="SoilMoistureItaly_{datetime}.tif",
              datetime_format_grid="%Y%m%d%H0000",
              data_sub_path_grid=None,
              file_name_tmpl_ts=None,
              datetime_format_ts=None,
              cell_format_ts='%04d',
              data_sub_path_ts=None,
              reset_ts=True,
              input_grid=None, target_grid=None, bbox=None,
              img_buffer=100,):

    # -------------------------------------------------------------------------------------
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

    # set ts
    dataset_ts_root = compose_sub_path(dataset_ts_root, run_date, data_sub_path_ts)
    file_name_tmpl_ts = compose_file_name(
        file_name_tmpl_ts, file_time=run_date, time_placeholder='datetime_destination')
    os.makedirs(dataset_ts_root, exist_ok=True)
    if reset_ts:
        reset_folder(dataset_ts_root)
    # set start and end hour
    start_step, end_step = start_date.hour, end_date.hour
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Section 1 - datasets type and format
    logging.info(' -----> 1) Define datasets type and format ... ')

    if get_filetype(dataset_grid_root) == "tiff":
        if input_grid is not None:
            dataset_grid_driver = smap_ds(
                product=product,
                data_path=dataset_grid_root,
                grid_path=grid_path,
                data_sub_path=data_sub_path_grid,
                file_name_tmpl=file_name_tmpl_grid,
                datetime_tmpl=datetime_format_grid,
                parameter=parameters, subgrid=input_grid, array_1D=True,
                start_step = 0, end_step = 23
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
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Section 2 - datasets metadata
    logging.info(' -----> 2) Define datasets metadata ... ')

    # create ts folder
    global_attr = {"product": product}

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
        target_grid = subgrid4bbox(input_grid, *bbox)
    else:
        target_grid = deepcopy(input_grid)

    logging.info(' -----> 2) Define datasets metadata ... DONE')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Section 3 - datasets conversion from grid to time-series
    logging.info(' -----> 3) Convert datasets from grid to time-series ... ')

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
    )
    # convert grid to ts format
    drv_reshuffle.calc()

    logging.info(' -----> 3) Convert datasets from grid to time-series ... DONE')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to parse argument(s)
def parse_args(args):
    """Parse command line parameters.

    WARNING TO DEVELOPER: if you add any parameter remember that this method
    parses args in order, so if you have an unknown number of arguments in
    your configuration file, then you should put it at the end.
    Always check with the debugger how args is built and sorted to prevent
    any errors.

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
        description="Convert SMAP data to time series format."
    )

    parser.add_argument(
        "product",
        default="product",
        help="Product name",
    )

    parser.add_argument(
        "flags",
        nargs=2,
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
        type=str,
        help="Root of local filesystem where the data is stored.",
    )

    parser.add_argument(
        "dataset_ts_root",
        type=str,
        help="Root of local filesystem where the timeseries should be stored.",
    )

    parser.add_argument(
        "start",
        type=mkdate,
        help="time_start. Either in format YYYY-MM-DD or " "YYYY-MM-DDTHH:MM.",
    )

    parser.add_argument(
        "end",
        type=mkdate,
        help=(
            "time_end. Either in format YYYY-MM-DD or " "YYYY-MM-DDTHH:MM."),
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
            " in the SMAP land mask (faster and less/smaller files)"
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
        default=4,
        help=(
            "How many images to read at once. Bigger "
            "numbers make the conversion faster but "
            "consume more memory."
        ),
    )

    # set defaults that can not be handled by argparse
    args = parser.parse_args(args)

    return args
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
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
        ' ---> Convert SMAP datasets from {} to {} into folder {} ... '.format(
            args.start.isoformat(), args.end.isoformat(), args.dataset_ts_root
        )
    )

    # method to load grid
    logging.info(' ----> Get reference grid ... ')
    input_grid = cell_grid(grid_path=args.grid_path)
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
        reset_ts=args.flags[1]
    )

    logging.info(' ----> Run reshuffle algorithm to convert grid to time-series ... DONE')

    # info end
    logging.info(
        ' ---> Convert SMAP datasets from {} to {} into folder {} ... DONE'.format(
            args.start.isoformat(), args.end.isoformat(), args.dataset_ts_root
        )
    )
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to run from command-line
def run():
    main(sys.argv[1:])
# -------------------------------------------------------------------------------------
