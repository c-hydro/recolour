"""
Library Features:

Name:          lib_figure_fx_results_data
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np

import cartopy

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pylab as plt

from lib_figure_cartography import CartoMap, oversample
from lib_figure_fx_tools import get_field, configure_plot_kwargs

coast = cartopy.feature.NaturalEarthFeature(
    category='physical', scale='50m', name='coastline',
    facecolor='none', edgecolor='#777777')

land = cartopy.feature.NaturalEarthFeature(
    category='physical', scale='50m', name='land',
    facecolor='#aaaaaa')

logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fiona').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results data
def plot_results_data(fig_file_name, fig_dframe, fig_kwargs, fig_committed_area=False):

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    plot_kwargs = configure_plot_kwargs(fig_kwargs)

    # get figure lons, lats and data
    fig_lons = fig_dframe['lon'].values
    fig_lats = fig_dframe['lat'].values
    fig_data = fig_dframe[fig_kwargs['parameter']].values

    tmp_dframe = fig_dframe.loc[fig_dframe['cell'] == 808]

    # clip data (if activated)
    fig_clip = get_field(fig_kwargs, field_name='clip')
    if fig_clip is not None:
        fig_data = np.clip(fig_data, fig_clip[0], fig_clip[1])
    # organize nan(s)
    idx_nan = np.isnan(fig_data)
    fig_data[idx_nan] = 100.5

    if fig_committed_area:
        fig_data[fig_dframe['committed_area'].values == 0] = 100.9

    # define map
    fig_size = get_field(fig_kwargs, field_name='fig_size')
    fig_map_pos = get_field(fig_kwargs, field_name='map_pos')
    m = CartoMap(figsize=fig_size, map_pos=fig_map_pos)

    # method to oversample data
    fig_data_extent = get_field(fig_kwargs, field_name='data_extent')
    fig_grid_sampling = get_field(fig_kwargs, field_name='grid_sampling')
    fig_max_distance = get_field(fig_kwargs, field_name='max_dist')

    img, reg_grid = oversample(
        fig_lons, fig_lats, fig_data, fig_data_extent,
        grid_sampling=fig_grid_sampling,
        max_dist=fig_max_distance)

    # masked should be not shown
    img[img.mask] = 101.5

    # compute figure bin(s)
    fig_vmin = get_field(fig_kwargs, field_name='vmin')
    fig_vmax = get_field(fig_kwargs, field_name='vmax')
    fig_cmap = get_field(fig_kwargs, field_name='cmap')

    figure_bins = np.linspace(fig_vmin, fig_vmax, len(fig_cmap.colors) + 1)
    figure_bins = np.concatenate((figure_bins, np.array([100.8, 101, np.inf])))

    x = np.digitize(img, figure_bins)
    y = [*fig_cmap.colors, *[np.array([0.4, 0.4, 0.4, 1]),
                             # np.array([0.9098, 0.8253, 0.5176, 1]),
                             # np.array([0.866, 0.235, 0.235, 1]),
                             np.array([0.8, 0.8, 0.8, 1]),
                             np.array([0, 0, 0, 0])]]
    color_arr = np.vstack(y)

    fig_data_extent = get_field(fig_kwargs, field_name='data_extent')

    m.plot_img(color_arr[x - 1], fig_data_extent, **plot_kwargs)
    # m.ax.add_feature(land, lw=0.3, zorder=0.5)
    m.ax.add_feature(coast, lw=0.3, zorder=0.5)

    # define patch(es) for nan and land
    color = (0.4, 0.4, 0.4)
    patch_result_nan = mpatches.Patch(color=color, label='NaN')
    color = '#aaaaaa'
    patch_land = mpatches.Patch(color=color, label='Land')

    if fig_committed_area:
        # color = (0.9098, 0.8253, 0.5176)
        # color = (0.866, 0.235, 0.235)
        color = (0.8, 0.8, 0.8)
        patch_not_committed = mpatches.Patch(
            color=color, label='Non-committed area')

        # handles = [patch_land, patch_result_nan, patch_not_committed]
        handles = [patch_land, patch_result_nan, patch_not_committed]
    else:
        # handles = [patch_land, patch_result_nan]
        handles = [patch_land, patch_result_nan]

    # add legend to map
    m.ax.legend(handles=handles, loc='lower center',
                borderaxespad=0., bbox_to_anchor=(0, -.03, 1., .102),
                ncol=4, frameon=False, prop={'size': 4})
    # save figure
    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')
    m.export_fig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------
