"""
Library Features:

Name:          lib_utils_smap
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230628'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import os
import logging
import rasterio
import numpy as np
import netCDF4 as nc

# debugging
# from matplotlib import pyplot as plt

logging.getLogger('rasterio').setLevel(logging.WARNING)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to create file grid
def create_file_grid(grid_path=None, data_path=None, data_ext='.tif', grid_update=True):
    """Creates nc grid file from first image in collection

    Takes the first image which is found by walking data_path and which matches
    the extension data_ext and produces a proper nc file containing the grid
    in the formatting wanted by TUW codes.

    :param grid_path: str, default None
        Leave None if your goal is to produce a nc grid file from a tif/tiff
    :param data_path: str, default None
        Path to folder of tif images
    :param data_ext: str, default '.tif'
    :param grid_update: bool, default True
        If True and grid_path is None, creates a nc grid file in that location
    :return:
        None, creates grid file at the location grid_path if not present
    """

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
            if file_path: break

        if file_path is None:
            logging.error(' ===> File for getting grid is not available')
            raise RuntimeError('Grid is needed by the procedure to convert files from grid to time-series format')

        # method to read file tiff
        file_values, file_geo_x_1d, file_geo_y_1d = read_file_tiff(file_path)

        # method to save file netcdf
        save_file_nc(grid_path, file_values, file_geo_x_1d, file_geo_y_1d)

        logging.info(' ---> Create grid reference for datasets ... DONE')

    else:
        logging.info(' ---> Create grid reference for datasets ... PREVIOUSLY DONE')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
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
    # file_geo_y_1d = np.arange(file_center_bottom, file_center_top + np.abs(file_res[1] / 2),
    #                           np.abs(file_res[1]), float)

    file_geo_x_2d, file_geo_y_2d = np.meshgrid(file_geo_x_1d, file_geo_y_1d)

    file_geo_y_upper, file_geo_y_lower = file_geo_y_2d[0, 0], file_geo_y_2d[-1, 0]
    if file_geo_y_lower > file_geo_y_upper:
        # file_geo_y_2d = np.flipud(file_geo_y_2d) # original
        file_values = np.flipud(file_values)

    ''' debug
    plt.figure()
    plt.imshow(file_dset.read(1), cmap='pink')
    plt.savefig('sample_grid.png', dpi=300)
    plt.close()
    '''
    return file_values, file_geo_x_1d, file_geo_y_1d
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
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

    file_var_geo_x[:] = file_geo_x_1d
    file_var_geo_y[:] = file_geo_y_1d
    file_var_values[0, :, :] = file_values

    ''' debug
    data = file_dset[var_name][:]
    plt.figure()
    plt.imshow(file_values, cmap='jet')
    plt.colorbar(label="Colorbar Label")
    plt.savefig('sample_nc.png')
    plt.close()
    '''

    file_dset.close()
# -------------------------------------------------------------------------------------
