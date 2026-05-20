"""
Library Features:

Name:           lib_img2cell_soilgrids
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260521'
Version:        '1.0.0'
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# libraries
import os
import logging
import numpy as np
import xarray as xr

from copy import deepcopy
from datetime import datetime

import pygeogrids.netcdf as grid2nc
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to handle img2cell errors
class Img2CellError(Exception):
    pass
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to convert static image datasets to cell datasets
class Img2Cell(object):

    def __init__(
            self,
            input_dataset,
            outputpath,
            input_grid=None,
            target_grid=None,
            variable_rename=None,
            cellsize_lat=5.0,
            cellsize_lon=5.0,
            filename_templ='{cell_n}.nc',
            cell_templ='%04d',
            gridname='grid.nc',
            global_attr=None,
            cell_attributes=None,
            cell_dtypes=None,
            zlib=True):

        self.imgin = input_dataset
        self.zlib = zlib

        if not hasattr(self.imgin, 'grid'):
            self.input_grid = input_grid
        else:
            self.input_grid = self.imgin.grid

        if self.input_grid is None and target_grid is None:
            raise ValueError(
                "Either the input dataset has to have a grid, "
                "input_grid has to be specified or "
                "target_grid has to be set"
            )

        self.target_grid = target_grid

        if self.target_grid is None:
            self.target_grid = self.input_grid

        if not hasattr(self.target_grid, "grid_points_for_cell"):
            self.target_grid = self.target_grid.to_cell_grid(
                cellsize_lat=cellsize_lat,
                cellsize_lon=cellsize_lon
            )

        self.outputpath = outputpath
        self.variable_rename = variable_rename
        self.filename_templ = filename_templ
        self.cell_tmpl = cell_templ
        self.gridname = gridname
        self.global_attr = global_attr
        self.cell_attributes = cell_attributes
        self.cell_dtypes = cell_dtypes

    def calc(self, parameters):

        logging.info(' ------> Organize static cell datasets ... ')

        start_time = datetime.now()

        os.makedirs(self.outputpath, exist_ok=True)

        # save grid information in file
        grid2nc.save_grid(
            os.path.join(self.outputpath, self.gridname),
            self.target_grid
        )

        # read all static images
        img_data, img_metadata = self.img_bulk(parameters)

        # iterate over cell(s) in target grid
        for cell in self.target_grid.get_cells():

            logging.info(' -------> Dump cell "' + str(cell) + '" ... ')

            # get cell info
            cell_gpis, cell_lons, cell_lats = self.target_grid.grid_points_for_cell(cell)

            # look where in the subset the data is
            cell_index = np.where(
                cell == self.target_grid.activearrcell
            )[0]

            if cell_index.size == 0:
                raise Img2CellError('cell not found in grid subset')

            data = {}

            for key in img_data:

                # rename variable in output dataset
                if self.variable_rename is None:
                    var_new_name = str(key)
                else:
                    var_new_name = self.variable_rename[key]

                output_array = deepcopy(img_data[key][cell_index])

                # change dtype of output cell variable
                if self.cell_dtypes is not None:
                    if type(self.cell_dtypes) == dict:
                        output_dtype = self.cell_dtypes[key]
                    else:
                        output_dtype = self.cell_dtypes

                    output_array = output_array.astype(output_dtype)

                data[var_new_name] = (
                    ("locations",),
                    output_array.astype(np.float32)
                )

            cell_string = self.cell_tmpl % cell
            filename_cell = self.filename_templ.format(cell_n=cell_string)
            file_path = os.path.join(self.outputpath, filename_cell)

            dataset = xr.Dataset(
                data_vars=data,
                coords={
                    "gpi": (
                        ("locations",),
                        cell_gpis.astype(np.int32)
                    ),
                    "longitude": (
                        ("locations",),
                        cell_lons.astype(np.float32)
                    ),
                    "latitude": (
                        ("locations",),
                        cell_lats.astype(np.float32)
                    )
                },
                attrs={
                    **(self.global_attr or {}),
                    "cell": int(cell),
                    "geospatial_lat_min": float(np.min(cell_lats)),
                    "geospatial_lat_max": float(np.max(cell_lats)),
                    "geospatial_lon_min": float(np.min(cell_lons)),
                    "geospatial_lon_max": float(np.max(cell_lons)),
                }
            )

            for var_name, var_attrs in (self.cell_attributes or {}).items():
                if var_name in dataset:
                    dataset[var_name].attrs.update(var_attrs)

            encoding = {}

            for var_name in data:
                encoding[var_name] = {
                    "zlib": self.zlib,
                    "complevel": 4,
                    "_FillValue": -9999.0,
                    "dtype": "float32"
                }

            dataset.to_netcdf(
                file_path,
                encoding=encoding
            )

            dataset.close()

            logging.info(' -------> Dump cell "' + str(cell) + '" ... DONE')

        elapsed_time = datetime.now() - start_time

        logging.info(
            ' ------> Organize static cell datasets ... DONE '
            '(Elapsed_Time: ' + str(elapsed_time) + ')'
        )

    def img_bulk(self, parameters):

        img_data = {}
        img_metadata = {}

        for parameter in parameters:

            logging.info(' -------> Read static variable "' + parameter + '" ... ')

            img = self.imgin.read(parameter)

            for key in img.data:

                if key not in img_data:
                    img_data[key] = img.data[key]

                if hasattr(img, "metadata"):
                    img_metadata[key] = img.metadata.get(key, {})

            logging.info(' -------> Read static variable "' + parameter + '" ... DONE')

        return img_data, img_metadata
# -------------------------------------------------------------------------------------