"""
Library Features:

Name:          lib_figure_fx
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

    plot_kwargs = _configure_plot_kwargs(fig_kwargs)

    # get figure lons, lats and data
    fig_lons = fig_dframe['lon'].values
    fig_lats = fig_dframe['lat'].values
    fig_data = fig_dframe[fig_kwargs['parameter']].values

    # clip data (if activated)
    fig_clip = _get_field(fig_kwargs, field_name='clip')
    if fig_clip is not None:
        fig_data = np.clip(fig_data, fig_clip[0], fig_clip[1])
    # organize nan(s)
    idx_nan = np.isnan(fig_data)
    fig_data[idx_nan] = 100.5

    if fig_committed_area:
        fig_data[fig_dframe['committed_area'].values == 0] = 100.9

    # define map
    fig_size = _get_field(fig_kwargs, field_name='fig_size')
    fig_map_pos = _get_field(fig_kwargs, field_name='map_pos')
    m = CartoMap(figsize=fig_size, map_pos=fig_map_pos)

    # method to oversample data
    fig_data_extent = _get_field(fig_kwargs, field_name='data_extent')
    fig_grid_sampling = _get_field(fig_kwargs, field_name='grid_sampling')
    fig_max_distance = _get_field(fig_kwargs, field_name='max_dist')

    img, reg_grid = oversample(
        fig_lons, fig_lats, fig_data, fig_data_extent,
        grid_sampling=fig_grid_sampling,
        max_dist=fig_max_distance)

    # masked should be not shown
    img[img.mask] = 101.5

    # compute figure bin(s)
    fig_vmin = _get_field(fig_kwargs, field_name='vmin')
    fig_vmax = _get_field(fig_kwargs, field_name='vmax')
    fig_cmap = _get_field(fig_kwargs, field_name='cmap')

    figure_bins = np.linspace(fig_vmin, fig_vmax, len(fig_cmap.colors) + 1)
    figure_bins = np.concatenate((figure_bins, np.array([100.8, 101, np.inf])))

    x = np.digitize(img, figure_bins)
    y = [*fig_cmap.colors, *[np.array([0.4, 0.4, 0.4, 1]),
                             # np.array([0.9098, 0.8253, 0.5176, 1]),
                             # np.array([0.866, 0.235, 0.235, 1]),
                             np.array([0.8, 0.8, 0.8, 1]),
                             np.array([0, 0, 0, 0])]]
    color_arr = np.vstack(y)

    fig_data_extent = _get_field(fig_kwargs, field_name='data_extent')

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
        handles = [patch_result_nan, patch_not_committed]
    else:
        # handles = [patch_land, patch_result_nan]
        handles = [patch_result_nan]

    # add legend to map
    m.ax.legend(handles=handles, loc='lower center',
                borderaxespad=0., bbox_to_anchor=(0, -.03, 1., .102),
                ncol=4, frameon=False, prop={'size': 4})
    # save figure
    fig_dpi = _get_field(fig_kwargs, field_name='fig_dpi')
    m.export_fig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results pie
def plot_results_pie(fig_file_name, fig_perc_comm, fig_perc_global, fig_kwargs):

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    fig_size = _get_field(fig_kwargs, field_name='fig_size')
    fig, ax = plt.subplots(1, 2, figsize=fig_size)

    size = 0.3
    cmap = _get_field(fig_kwargs, field_name='cmap')
    if isinstance(cmap, list):
        comm_colors = cmap[0]
        global_colors = cmap[0]
    else:
        cmap = plt.get_cmap('tab20c')
        comm_colors = cmap([4, 2, 1, 0, 16])
        global_colors = cmap([13, 11, 10, 9, 16])

    wedges, texts = ax[0].pie(
        fig_perc_comm, radius=0.6, colors=comm_colors,
        wedgeprops=dict(width=size, edgecolor='w'))

    ax[0].set(aspect="equal")

    labels = [
        '{:3.1f}%: -1 $\leq$ R < 0.5'.format(fig_perc_comm[0]),
        '{:3.1f}%: 0.5 $\leq$ R < 0.65'.format(fig_perc_comm[1]),
        '{:3.1f}%: 0.65 $\leq$ R < 0.8'.format(fig_perc_comm[2]),
        '{:3.1f}%: 0.8 $\leq$ R $\leq$ 1'.format(fig_perc_comm[3]),
        '{:3.1f}%: NaN'.format(fig_perc_comm[4])]

    ax[0].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                 [labels[3], labels[2], labels[1], labels[0], labels[4]],
                 title="Committed area", loc="center",
                 fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    wedges, texts = ax[1].pie(
        fig_perc_global, radius=0.6, colors=global_colors,
        wedgeprops=dict(width=size, edgecolor='w'))

    ax[1].set(aspect="equal")

    if fig_perc_global.__len__() > 0:
        labels = [
            '{:3.1f}%: -1 $\leq$ R < 0.5'.format(fig_perc_global[0]),
            '{:3.1f}%: 0.5 $\leq$ R < 0.65'.format(fig_perc_global[1]),
            '{:3.1f}%: 0.65 $\leq$ R < 0.8'.format(fig_perc_global[2]),
            '{:3.1f}%: 0.8 $\leq$ R $\leq$ 1'.format(fig_perc_global[3]),
            '{:3.1f}%: NaN'.format(fig_perc_global[4])]

        ax[1].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                     [labels[3], labels[2], labels[1], labels[0], labels[4]],
                     title="Global", loc="center",
                     fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    # plt.subplots_adjust(hspace=0.5, left=0.214, bottom=0, right=0.98, top=0.78, wspace=0.512)

    plt.subplots_adjust(hspace=0.5, left=0.28, bottom=0, right=1, top=0.78, wspace=0.614)

    fig_dpi = _get_field(fig_kwargs, field_name='fig_dpi')
    fig.savefig(fig_file_name, dpi=fig_dpi)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure plot kwargs
def _configure_plot_kwargs(figure_settings, plot_keys=None):

    if plot_keys is None:
        plot_keys = ['cb_label', 'show_cb', 'title', 'vmin', 'vmax', 'cmap',
                     'title_fontsize', 'title_fontweight', 'cb_fontsize', 'cb_fontweight']

    plot_settings = {}
    for figure_key, figure_value in figure_settings.items():
        if figure_key in plot_keys:
            plot_settings[figure_key] = figure_value
    return plot_settings
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get field
def _get_field(field_obj, field_name):
    if field_name in list(field_obj.keys()):
        return field_obj[field_name]
    else:
        logging.error(' ===> Field name "' + field_name + '" is not available in the field obj')
        raise KeyError('Check the figure settings and add the field needed by the algorithm')
# ----------------------------------------------------------------------------------------------------------------------
