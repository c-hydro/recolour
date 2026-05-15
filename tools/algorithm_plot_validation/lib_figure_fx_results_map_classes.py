"""
Library Features:

Name:          lib_figure_fx_results_map_classes
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260514'
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
# method to plot results classes pearson
def plot_results_classes_pearson(fig_file_name, fig_dframe, fig_kwargs, fig_committed_area=False):

    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    plot_kwargs = configure_plot_kwargs(fig_kwargs)

    fig_lons = fig_dframe['lon'].values
    fig_lats = fig_dframe['lat'].values

    fig_data = fig_dframe[fig_kwargs['parameter']].values
    fig_data = np.array(fig_data, copy=True)

    # NaN class
    fig_data[np.isnan(fig_data)] = 100.5

    # only mask non-committed if explicitly requested
    if fig_committed_area is True:
        fig_data[fig_dframe['committed_area'].values == 0] = 100.9

    fig_size = get_field(fig_kwargs, field_name='fig_size')
    if fig_size is None:
        fig_size = get_field(fig_kwargs, field_name='figsize')

    # title size
    fig_title_size = get_field(
        fig_kwargs, field_name='title_fontsize', field_mandatory=False)
    if fig_title_size is None:
        fig_title_size = 8

    # set colorbar tick fontsize
    fig_cbar_tick_size = get_field(
        fig_kwargs, field_name='cb_tick_fontsize', field_mandatory=False)
    if fig_cbar_tick_size is None:
        fig_cbar_tick_size = 4

    fig_map_pos = get_field(fig_kwargs, field_name='map_pos')

    m = CartoMap(figsize=fig_size, map_pos=fig_map_pos)

    fig_data_extent = get_field(fig_kwargs, field_name='data_extent')
    fig_grid_sampling = get_field(fig_kwargs, field_name='grid_sampling')
    fig_max_distance = get_field(fig_kwargs, field_name='max_dist')

    img, reg_grid = oversample(
        fig_lons, fig_lats, fig_data, fig_data_extent,
        grid_sampling=fig_grid_sampling,
        max_dist=fig_max_distance
    )

    # masked / no-data pixels
    img[img.mask] = 101.5

    fig_cmap = get_field(fig_kwargs, field_name='cmap')

    if fig_cmap is None:
        cmap_base = plt.get_cmap('tab20c')
        fig_cmap = ListedColormap(
            cmap_base([0, 1, 2, 4, 9, 10, 11, 13]),
            N=8
        )

    # bins for class values:
    # committed:
    #   0.5 -> PR >= 0.8
    #   1.5 -> 0.65 <= PR < 0.8
    #   2.5 -> 0.5 <= PR < 0.65
    #   3.5 -> PR < 0.5
    #
    # non-committed:
    #   4.5 -> PR >= 0.8
    #   5.5 -> 0.65 <= PR < 0.8
    #   6.5 -> 0.5 <= PR < 0.65
    #   7.5 -> PR < 0.5
    figure_bins = np.array([
        0, 1, 2, 3, 4, 5, 6, 7, 8,
        100.8, 101, np.inf
    ])

    x = np.digitize(img, figure_bins)

    color_arr = np.vstack([
        fig_cmap.colors,
        np.array([0.4, 0.4, 0.4, 1]),  # NaN
        np.array([0.8, 0.8, 0.8, 1]),  # non-committed mask
        np.array([0, 0, 0, 0])         # transparent / masked
    ])

    m.plot_img(color_arr[x - 1], fig_data_extent, **plot_kwargs)

    # set title fontsize
    if m.ax.get_title():
        m.ax.set_title(m.ax.get_title(), fontsize=fig_title_size)

    # last axis is usually the colorbar axis
    if len(m.fig.axes) > 1:
        cbar_ax = m.fig.axes[-1]
        cbar_ax.tick_params(labelsize=fig_cbar_tick_size)

    m.ax.add_feature(coast, lw=0.3, zorder=0.5)

    patch_result_nan = mpatches.Patch(
        color=(0.4, 0.4, 0.4),
        label='NaN'
    )

    handles = [patch_result_nan]

    if fig_committed_area is True:
        patch_not_committed = mpatches.Patch(
            color=(0.8, 0.8, 0.8),
            label='Non-committed area'
        )
        handles.append(patch_not_committed)

    m.ax.legend(
        handles=handles,
        loc='lower center',
        borderaxespad=0.,
        bbox_to_anchor=(0, -.03, 1., .102),
        ncol=4,
        frameon=False,
        prop={'size': 4}
    )

    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')
    if fig_dpi is None:
        fig_dpi = 300

    m.export_fig(fig_file_name, dpi=fig_dpi)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results classes snr
def plot_results_classes_snr(fig_file_name, fig_dframe, fig_kwargs, fig_committed_area=False):

    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    plot_kwargs = configure_plot_kwargs(fig_kwargs)

    fig_lons = fig_dframe['lon'].values
    fig_lats = fig_dframe['lat'].values

    fig_data = fig_dframe[fig_kwargs['parameter']].values
    fig_data = np.array(fig_data, copy=True)

    # NaN class
    fig_data[np.isnan(fig_data)] = 100.5

    # only mask non-committed if explicitly requested
    if fig_committed_area is True:
        fig_data[fig_dframe['committed_area'].values == 0] = 100.9

    fig_size = get_field(fig_kwargs, field_name='fig_size')
    if fig_size is None:
        fig_size = get_field(fig_kwargs, field_name='figsize')

    # title size
    fig_title_size = get_field(
        fig_kwargs, field_name='title_fontsize', field_mandatory=False)
    if fig_title_size is None:
        fig_title_size = 8

    # set colorbar tick fontsize
    fig_cbar_tick_size = get_field(
        fig_kwargs, field_name='cb_tick_fontsize', field_mandatory=False)
    if fig_cbar_tick_size is None:
        fig_cbar_tick_size = 4

    fig_map_pos = get_field(fig_kwargs, field_name='map_pos')

    m = CartoMap(figsize=fig_size, map_pos=fig_map_pos)

    fig_data_extent = get_field(fig_kwargs, field_name='data_extent')
    fig_grid_sampling = get_field(fig_kwargs, field_name='grid_sampling')
    fig_max_distance = get_field(fig_kwargs, field_name='max_dist')

    img, reg_grid = oversample(
        fig_lons, fig_lats, fig_data, fig_data_extent,
        grid_sampling=fig_grid_sampling,
        max_dist=fig_max_distance
    )

    img[img.mask] = 101.5

    fig_cmap = get_field(fig_kwargs, field_name='cmap')

    # bins for values:
    # 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5
    figure_bins = np.array([
        0, 1, 2, 3, 4, 5, 6, 7, 8,
        100.8, 101, np.inf
    ])

    x = np.digitize(img, figure_bins)

    color_arr = np.vstack([
        fig_cmap.colors,
        np.array([0.4, 0.4, 0.4, 1]),  # NaN
        np.array([0.8, 0.8, 0.8, 1]),  # non-committed mask
        np.array([0, 0, 0, 0])         # transparent
    ])

    m.plot_img(color_arr[x - 1], fig_data_extent, **plot_kwargs)

    # set title fontsize
    if m.ax.get_title():
        m.ax.set_title(m.ax.get_title(), fontsize=fig_title_size)

    # last axis is usually the colorbar axis
    if len(m.fig.axes) > 1:
        cbar_ax = m.fig.axes[-1]
        cbar_ax.tick_params(labelsize=fig_cbar_tick_size)

    m.ax.add_feature(coast, lw=0.3, zorder=0.5)

    patch_result_nan = mpatches.Patch(
        color=(0.4, 0.4, 0.4),
        label='NaN'
    )

    handles = [patch_result_nan]

    if fig_committed_area is True:
        patch_not_committed = mpatches.Patch(
            color=(0.8, 0.8, 0.8),
            label='Non-committed area'
        )
        handles.append(patch_not_committed)

    m.ax.legend(
        handles=handles,
        loc='lower center',
        borderaxespad=0.,
        bbox_to_anchor=(0, -.03, 1., .102),
        ncol=4,
        frameon=False,
        prop={'size': 4}
    )

    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')
    if fig_dpi is None:
        fig_dpi = 300

    m.export_fig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------
