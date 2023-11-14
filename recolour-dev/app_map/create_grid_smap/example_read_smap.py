#!/usr/bin/python3

import pyresample
import os
import pyproj
import numpy as np
import xarray as xr
import cartopy
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from cartopy.mpl.geoaxes import GeoAxes
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import AxesGrid
from pygeogrids.grids import BasicGrid, genreg_grid

from example_read_tools import read_smap_l2, read_grid_data

import matplotlib.pylab as plt

# -------------------------------------------------------------------------------------
# main
def main(file_name, grid_name):

    # GRID CASE ITALY
    grid_data, grid_attrs = read_grid_data(grid_name)
    grid_lons_1d, grid_lats_1d = grid_data['longitude'].values, grid_data['latitude'].values

    # TEST CASE
    grid_res = 0.009
    grid_x_min, grid_x_max = 110, 150
    grid_y_min, grid_y_max = 10, 50
    grid_lons_1d = np.arange(grid_x_min, grid_x_max + grid_res, grid_res)
    grid_lats_1d = np.arange(grid_y_min, grid_y_max + grid_res, grid_res)

    # define grid mesh
    grid_lons_2d, grid_lats_2d = np.meshgrid(grid_lons_1d, grid_lats_1d)
    # grid obj using pyresample definition
    grid_obj = pyresample.geometry.GridDefinition(lats=grid_lats_2d, lons=grid_lons_2d)
    # define grid data 2d
    grid_sm_2d = np.zeros(shape=(grid_lons_2d.shape[0], grid_lats_2d.shape[1]))
    grid_sm_2d[:, :] = np.nan

    # SMAP DATA
    file_data, file_geo_x, file_geo_y, file_col, file_row = read_smap_l2(file_name)

    file_reference = {
        'epsg': 6933, 'x_min': -17367530.45, 'y_max': 7314540.83, 'res': 9008.05, 'n_cols': 3856, 'n_rows': 1624}
    file_proj = pyproj.Proj(file_reference['epsg'])

    file_x = file_reference['x_min'] + file_col * file_reference['res'] + file_reference['res'] / 2
    file_y = file_reference['y_max'] - file_row * file_reference['res'] - file_reference['res'] / 2

    # get data, lons and lats
    values_raw = file_data.data
    lons_raw, lats_raw = file_proj(file_x, file_y, inverse=True)
    idx_finite = np.argwhere(np.isfinite(values_raw))[:, 0]

    # select data finite
    values_finite = values_raw[idx_finite]
    lons_finite = lons_raw[idx_finite]
    lats_finite = lats_raw[idx_finite]
    # get geographical window
    lon_finite_min, lon_finite_max = np.min(lons_finite), np.max(lons_finite)
    lat_finite_min, lat_finite_max = np.min(lats_finite), np.max(lats_finite)

    file_obj = pyresample.geometry.SwathDefinition(lons=lons_finite, lats=lats_finite)

    # join point index to grid index
    valid_input_index, valid_output_index, index_array, distance_array = pyresample.kd_tree.get_neighbour_info(
        source_geo_def=grid_obj, target_geo_def=file_obj, radius_of_influence=9000,
        neighbours=1)

    # select data valid
    values_valid = values_finite[valid_output_index]
    lons_valid = lons_finite[valid_output_index]
    lats_valid = lats_finite[valid_output_index]
    # convert 1d index to 2d index
    index_array_2d = np.unravel_index(index_array, grid_obj.shape)

    # join point data to grid data
    grid_sm_2d[index_array_2d[0], index_array_2d[1]] = values_valid

    # plot grid
    plot_crs = cartopy.crs.Mercator()
    data_crs = cartopy.crs.PlateCarree()

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=plot_crs)
    ax.set_title('sm grid')

    ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAND, facecolor='#aaaaaa')
    ax.set_extent([lon_finite_min, lon_finite_max, lat_finite_max, lat_finite_min])

    sc = ax.pcolormesh(grid_lons_2d, grid_lats_2d, grid_sm_2d, zorder=3,
                       transform=data_crs, vmin=0, vmax=1)

    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0,
                        0.02, ax.get_position().height])

    cbar = fig.colorbar(sc, ax=ax, cax=cax)
    cbar.set_label('Degree of Saturation (%)')

    plt.show()

    # plot swath
    plot_crs = cartopy.crs.Mercator()
    data_crs = cartopy.crs.PlateCarree()

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=plot_crs)
    ax.set_title('sm swaths')

    ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    ax.add_feature(cartopy.feature.LAND, facecolor='#aaaaaa')
    ax.set_extent([lon_finite_min, lon_finite_max, lat_finite_max, lat_finite_min])

    sc = ax.scatter(lons_valid, lats_valid,
                    c=values_valid, zorder=3, marker='s', s=2,
                    transform=data_crs, vmin=0, vmax=1)

    cax = fig.add_axes([ax.get_position().x1 + 0.01, ax.get_position().y0,
                        0.02, ax.get_position().height])

    cbar = fig.colorbar(sc, ax=ax, cax=cax)
    cbar.set_label('Degree of Saturation (%)')

    plt.show()

    print('')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# call script from external library
if __name__ == "__main__":

    file_folder = '/home/fabio/Desktop/Recolour_Workspace/recolour-ws/smap_test_georef/'
    file_name = 'SMAP_L2_SM_P_E_46118_D_20230919T222213_R2023_005.h5'

    file_path = os.path.join(file_folder, file_name)

    grid_folder = '/home/fabio/Desktop/Recolour_Workspace/recolour-ws/domain/'
    grid_name = 'italy_sm_smap.tiff'

    grid_path = os.path.join(grid_folder, grid_name)

    main(file_path, grid_path)
# -------------------------------------------------------------------------------------
