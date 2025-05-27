"""
Library Features:

Name:          lib_utils_hmc
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230628'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import logging
import rasterio

import pickle
import numpy as np
import netCDF4 as nc

logging.getLogger('rasterio').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to select grid data origin
def select_origin_grid(data_path: str, data_origin: (str, None) =None, default_origin: str ='.nc') -> str:

    if data_origin is None:
        ext_origin = default_origin
    elif data_origin == 'source_datasets':

        extension = None
        for root, dirs, files in os.walk(data_path):
            for file in files:
                extension = os.path.splitext(file)[1]
                if extension in ['.nc', '.tif', '.tiff', '.nc4', '.netcdf']:
                    break

        if extension is not None:
            ext_origin = extension
        else:
            logging.error(' ===> Grid data origin is not selected')
            raise RuntimeError('Check the supported origin datasets to correctly set the grid data origin')

    else:
        logging.error(' ===> Grid data origin "' + data_origin + '" is not supported')
        raise NotImplemented('Case not supported for grid file')


    return ext_origin
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create file grid
def create_file_grid(grid_path, data_path, data_ext='.tif', grid_update=True):

    logging.info(' ---> Create grid reference for datasets ... ')

    if grid_update:
        if os.path.exists(grid_path):
            os.remove(grid_path)

    if not os.path.exists(grid_path):

        grid_folder, grid_name = os.path.split(grid_path)
        os.makedirs(grid_folder, exist_ok=True)

        file_path = None
        for root, dirs, files in os.walk(data_path):
            for file in files:

                if file.endswith(data_ext):
                    file_path = os.path.join(root, file)
                    break

        if file_path is None:
            logging.error(' ===> File for getting grid is not available')
            raise RuntimeError('Grid is needed by the procedure to convert files from grid to time-series format')

        # method to read file tiff
        # method definition below
        file_values, file_geo_x_1d, file_geo_y_1d = read_file_tiff(file_path)
        # method to save file netcdf
        # method definition below
        save_file_nc(grid_path, file_values, file_geo_x_1d, file_geo_y_1d)

        logging.info(' ---> Create grid reference for datasets ... DONE')

    else:
        logging.info(' ---> Create grid reference for datasets ... PREVIOUSLY DONE')

# ----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------
# method to read file tiff
def read_file_tiff(file_name):

    file_dset = rasterio.open(file_name)

    file_bounds, file_res, file_transform = file_dset.bounds, file_dset.res, file_dset.transform
    file_data = file_dset.read()

    file_values = np.float32(file_data[0, :, :])

    file_center_right = file_bounds.right - (file_res[0] / 2)
    file_center_left = file_bounds.left + (file_res[0] / 2)
    file_center_top = file_bounds.top - (file_res[1] / 2)
    file_center_bottom = file_bounds.bottom + (file_res[1] / 2)

    file_geo_x_1d = np.arange(file_center_left, file_center_right + np.abs(file_res[0] / 2),
                              np.abs(file_res[0]), float)
    file_geo_y_1d = np.flip(np.arange(file_center_bottom, file_center_top + np.abs(file_res[1] / 2),
                                      np.abs(file_res[1]), float), axis=0)
    file_geo_x_2d, file_geo_y_2d = np.meshgrid(file_geo_x_1d, file_geo_y_1d)

    file_geo_y_upper, file_geo_y_lower = file_geo_y_2d[0, 0], file_geo_y_2d[-1, 0]
    if file_geo_y_lower > file_geo_y_upper:
        file_geo_y_2d = np.flipud(file_geo_y_2d)
        file_values = np.flipud(file_values)

    return file_values, file_geo_x_1d, file_geo_y_1d
# ------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------
# method to create file netcdf
def save_file_nc(file_name, file_values, file_geo_x_1d, file_geo_y_1d,
                 var_name='soil_moisture', var_geo_x='longitude', var_geo_y='latitude', var_time='time'):

    dim_geo_x, dim_geo_y = file_geo_x_1d.shape[0], file_geo_y_1d.shape[0]

    file_dset = nc.Dataset(file_name, 'w', format='NETCDF4')

    file_dim_time = file_dset.createDimension(var_time, None)
    file_dim_geo_y = file_dset.createDimension(var_geo_y, dim_geo_y)
    file_dim_geo_x = file_dset.createDimension(var_geo_x, dim_geo_x)

    file_var_time = file_dset.createVariable(var_time, 'f4', (var_time,))
    file_var_geo_y = file_dset.createVariable(var_geo_y, 'f4', (var_geo_y,))
    file_var_geo_x = file_dset.createVariable(var_geo_x, 'f4', (var_geo_x,))
    file_var_values = file_dset.createVariable(var_name, 'f4', (var_time, var_geo_y, var_geo_x,))
    file_var_values.units = '-'

    file_var_geo_y[:] = file_geo_y_1d
    file_var_geo_x[:] = file_geo_x_1d
    file_var_values[0, :, :] = file_values

    file_dset.close()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read data obj
def read_obj(file_name):
    data = None
    if os.path.exists(file_name):
        data = pickle.load(open(file_name, "rb"))
    return data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write data obj
def write_obj(file_name, data):
    if os.path.exists(file_name):
        os.remove(file_name)
    with open(file_name, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
# ----------------------------------------------------------------------------------------------------------------------
