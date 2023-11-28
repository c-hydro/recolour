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
# method to plot committed area
def plot_committed_area(fig_file_name, fig_dframe, fig_kwargs):

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    plot_kwargs = _configure_plot_kwargs(fig_kwargs)

    # get figure lons, lats and data
    fig_lons = fig_dframe['lon'].values
    fig_lats = fig_dframe['lat'].values

    fig_land = fig_dframe['land_flag'] == 1
    fig_data = fig_dframe[fig_kwargs['parameter']][fig_land].values - 0.5
    fig_data = fig_data[:, 0]

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

    figure_bins = np.array([-0.5, 0.2, 0.5, 1.5])
    figure_bins = np.concatenate((figure_bins, np.array([100, 101, np.inf])))

    x = np.digitize(img, figure_bins)
    y = [*fig_cmap.colors, *[np.array([0.4, 0.4, 0.4, 1]),
                             np.array([0.9098, 0.8253, 0.5176, 1]),
                             # np.array([0.866, 0.235, 0.235, 1]),
                             # np.array([0.8, 0.8, 0.8, 1]),
                             np.array([0, 0, 0, 0])]]
    color_arr = np.vstack(y)

    fig_data_extent = _get_field(fig_kwargs, field_name='data_extent')

    m.plot_img(color_arr[x - 1], fig_data_extent, **plot_kwargs)
    m.ax.add_feature(coast, lw=0.3, zorder=0.5)

    patch_committed = mpatches.Patch(
        color=fig_cmap.colors[2], label='Committed area')

    patch_global = mpatches.Patch(
        color=fig_cmap.colors[0], label='Global')

    handles = [patch_committed, patch_global]

    # add legend to map
    m.ax.legend(handles=handles, loc='lower center',
                borderaxespad=0., bbox_to_anchor=(0, -.12, 1., .102),
                ncol=4, frameon=False)

    # save figure
    fig_dpi = _get_field(fig_kwargs, field_name='fig_dpi')
    m.export_fig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------


'''
# ----------------------------------------------------------------------------------------------------------------------
# method to plot ssm mean
def plot_ssm_mean(filename):

    with netCDF4.Dataset(filename) as nc:
        gpi = nc.variables['location_id'][:]
        data = nc.variables['mean_sm'][:].astype(np.float32)

    data[data >= 100] = np.nan
    data = np.ma.array(data, mask=np.isnan(data))

    df = pd.DataFrame({'mean': data}, index=gpi)

    df_grid = read_committed()
    df = df_grid.join(df)
    df = df[df['land_flag'] == 1]

    figsize = get_figsize(455.24411)
    map_pos = [0.05, 0.16, 0.9, 0.76]
    cb_pos = [0.3, 0.14, 0.4, 0.025]
    data_extent = [-180, 180, -60, 85]
    grid_sampling = 0.05
    max_dist = max_distance

    title = 'Metop ASCAT SSM CDR v7 12.5 km (H119) - Mean 2007-2020'
    cb_label = 'Saturation (%)'
    vmin = 0
    vmax = 100
    cmap = cmcrameri.cm.roma
    cmap = ListedColormap(cmap(np.linspace(0, 255, 30, dtype=int)), N=30)

    plot_kwargs = {
        'cb_label': cb_label, 'show_cb': True, 'title': title,
        'vmin': vmin, 'vmax': vmax, 'cmap': cmap}

    settings_kwargs = {'parameter': 'mean', 'committed_area': False,
                       'figsize': figsize, 'map_pos': map_pos,
                       'cb_pos': cb_pos, 'data_extent': data_extent,
                       'grid_sampling': grid_sampling, 'max_dist': max_dist,
                       'cmap': cmap}

    settings_kwargs['fig_filename'] = '/data-write/RADAR/warp/h119_ssm_mean.png'
    # plot_results(df, settings_kwargs, plot_kwargs)

    settings_kwargs['committed_area'] = True
    settings_kwargs['fig_filename'] = '/data-write/RADAR/warp/h119_ssm_mean_committed_area.png'
    plot_results(df, settings_kwargs, plot_kwargs)
# ----------------------------------------------------------------------------------------------------------------------
'''


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
