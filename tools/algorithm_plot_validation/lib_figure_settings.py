"""
Library Features:

Name:          lib_figure_settings
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import cmcrameri
import numpy as np

from copy import deepcopy

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get figure size
def get_figure_size(fig_width_pt=252.0):
    """
    Parameters
    ----------
    fig_width_pt : float
        Get this from LaTeX using \showthe\columnwidth
        Unit: pt

    Returns
    -------
    figsize : tuple
        Figure size.
        Unit: inches
    """
    inches_per_pt = 1.0 / 72.27
    # Aesthetic ratio
    golden_mean = (np.sqrt(5) - 1.0) / 2.0
    # width in inches
    fig_width = fig_width_pt * inches_per_pt
    # height in inches
    fig_height = fig_width * golden_mean

    return fig_width, fig_height
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get figure cmap
def get_figure_cmap(cmap_type='vik_r', cmap_n=None):
    if cmap_type == 'vik_r':
        if cmap_n is None:
            cmap_n = 32
        cmap_def = cmcrameri.cm.vik_r
        cmap_colors = ListedColormap(cmap_def(np.linspace(0, 255, cmap_n, dtype=int)), N=cmap_n)
    elif cmap_type == 'bam':
        if cmap_n is None:
            cmap_n = 24
        cmap_def = cmcrameri.cm.bam
        cmap_colors = ListedColormap(cmap_def(np.linspace(0, 255, cmap_n, dtype=int)), N=cmap_n)
    elif cmap_type == 'tab20c_committed_area':
        if cmap_n is None:
            cmap_n = 3
        cmap_def = plt.get_cmap('tab20')
        cmap_colors = ListedColormap(cmap_def([2, 1, 0]), N=cmap_n)
    elif cmap_type == 'tab20c_map':
        if cmap_n is None:
            cmap_n = 8
        cmap_def = plt.get_cmap('tab20c')
        cmap_colors = ListedColormap(cmap_def([0, 1, 2, 4, 9, 10, 11, 13]), N=cmap_n)
    elif cmap_type == 'tab20c_pie':
        cmap = plt.get_cmap('tab20c')
        cmap_colors_comm = cmap([4, 2, 1, 0, 16])
        cmap_colors_global = cmap([13, 11, 10, 9, 16])

        # test
        # cmap = ListedColormap(cmap([0, 1, 2, 4, 9, 10, 11, 13, 16]), N=9)
        # comm_colors = cmap([3, 2, 1, 0, 8])
        # global_colors = cmap([7, 6, 5, 4, 8])

        cmap_colors = [cmap_colors_comm, cmap_colors_global]
    elif cmap_type is None:
        logging.warning(' ===> CMap type is not defined and CMap colors will be set to NoneType')
        cmap_colors = None
    else:
        logging.error(' ===> CMap type "' + cmap_type + '" is not expected by the algorithm')
        raise NotImplemented('Case not implemented yet')
    return cmap_colors
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize figure extent
def organize_figure_extent(fig_dframe, fig_ext_default=None, fig_ext_domain=True):

    if fig_ext_default:
        fig_ext_default = [-180, 180, -60, 85]

    if fig_ext_domain:
        fig_geo_x_arr, fig_geo_y_arr = fig_dframe['lon'].values, fig_dframe['lat'].values

        fig_geo_x_min, fig_geo_x_max = np.min(fig_geo_x_arr), np.max(fig_geo_x_arr)
        fig_geo_y_min, fig_geo_y_max = np.min(fig_geo_y_arr), np.max(fig_geo_y_arr)

        domain_limits = [fig_geo_x_min, fig_geo_x_max, fig_geo_y_min, fig_geo_y_max]
        domain_offset = [-1, +1, -1, +1]
        fig_extent_geo = [x+y for x,y in zip(domain_offset, domain_limits)]
    else:
        fig_extent_geo = fig_ext_default.copy()

    return fig_extent_geo
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize figure settings
def organize_figure_settings(
        figure_settings, figure_type='data',
        figure_filename_tmpl='pearson_r_reference_vs_y.png',
        figure_parameter='xy_pr', figure_committed_area=False,
        figure_parameter_data='data', figure_parameter_p_r='p_r', figure_parameter_r='r',
        figure_title_label='title', figure_x_label='x_label', figure_y_label='y_label',
        figure_colorbar_label='colorbar', figure_colorbar_show=True,
        figure_colorbar_extent=None, figure_colorbar_ticks=None,
        figure_cmap_type='vik_r', figure_cmap_n=32,
        figure_vmin=0, figure_vmax=1, figure_data_extent=None,
        figure_title_fontsize=14, figure_title_fontweight=None,
        figure_cbar_fontsize=6, figure_cbar_fontweight=None,
        figure_lim_min=-10, figure_lim_max=10,
        figure_lim_thr=0, figure_lim_target=3, figure_lim_optimal=6,
        figure_palette_type='Set2'):

    # set figure size ptn
    figure_size_pt = 252.0
    if 'figure_size' in list(figure_settings.keys()):
        figure_size_pt = figure_settings['figure_size']
    figure_size = get_figure_size(figure_size_pt)
    # set figure dpi
    figure_dpi = 500
    if 'figure_dpi' in list(figure_settings.keys()):
        figure_dpi = figure_settings['figure_dpi']
    # set figure map position
    figure_map_pos = [0.05, 0.16, 0.9, 0.76]
    if 'map_pos' in list(figure_settings.keys()):
        figure_map_pos = figure_settings['map_pos']
    # set figure colorbar position
    figure_cb_pos = [0.3, 0.14, 0.4, 0.025]
    if 'cb_pos' in list(figure_settings.keys()):
        figure_cb_pos = figure_settings['cb_pos']
    # set figure grid_sampling
    figure_grid_sampling = 0.05
    if 'grid_sampling' in list(figure_settings.keys()):
        figure_grid_sampling = figure_settings['grid_sampling']
    # set figure max_distance point(s)
    figure_max_distance = 12000
    if 'max_distance' in list(figure_settings.keys()):
        figure_max_distance = figure_settings['max_distance']
    # set figure clip
    figure_clip = None
    if 'clip' in list(figure_settings.keys()):
        figure_clip = figure_settings['clip']
        if (figure_clip is not None) and (not isinstance(figure_clip, list)):
            logging.error(' ===> Field "clip" must be defined by NoneType or list with two elements')
            raise TypeError('Defined "clip" in the setting file using a list or a NoneType object')

    # set figure data_extent
    if figure_data_extent is None:
        if 'data_extent_default' in list(figure_settings.keys()):
            figure_data_extent = figure_settings['data_extent_default']
        else:
            figure_data_extent = [-180, 180, -60, 85]

    # get figure cmap object
    figure_cmap_obj = get_figure_cmap(cmap_type=figure_cmap_type, cmap_n=figure_cmap_n)
    # organize figure committed area conditions
    figure_committed_area = organize_figure_committed_area(figure_committed_area)

    # set figure committed area and filename(s)
    figure_file_list = organize_figure_filename(
        file_name_tmpl=figure_filename_tmpl, file_filter_committed_area=figure_committed_area)

    # define figure file and committed area flag
    if figure_type == 'data':
        figure_file_obj = deepcopy(figure_file_list)
        figure_committed_area_obj = deepcopy(figure_committed_area)
    elif figure_type == 'stats_pearson_pie' or figure_type == 'stats_snr_pie':
        figure_file_obj = figure_file_list[0]
        figure_committed_area_obj = figure_committed_area[0]
    elif figure_type == 'stats_pearson_box' or figure_type == 'stats_snr_box':
        figure_file_obj = figure_file_list[0]
        figure_committed_area_obj = deepcopy(figure_committed_area)
    elif figure_type == 'committed_area':
        figure_file_obj = deepcopy(figure_file_list)[0]
        figure_committed_area_obj = True
    else:
        logging.error(' ===> Figure type "' + figure_type + '" is not supported')
        raise NotImplemented('Case not implemented yet')

    # adapt the ticks to expected python idxs
    if figure_colorbar_ticks is not None:
        if isinstance(figure_colorbar_ticks, list):
            figure_colorbar_ticks[1] = figure_colorbar_ticks[1] + 1

    # organize settings kwargs
    settings_kwargs = {
        'fig_filename': figure_file_obj, 'fig_committed_area': figure_committed_area_obj,
        'fig_size': figure_size, 'fig_dpi': figure_dpi,
        'title': figure_title_label,
        'x_label': figure_x_label, 'y_label': figure_y_label,
        'title_fontsize': figure_title_fontsize, 'title_fontweight': figure_title_fontweight,
        'cb_pos': figure_cb_pos,
        'cb_label': figure_colorbar_label, 'show_cb': figure_colorbar_show,
        'cb_extend': figure_colorbar_extent, 'cb_ticks': figure_colorbar_ticks,
        'cb_fontsize': figure_cbar_fontsize, 'cb_fontweight': figure_cbar_fontweight,
        'parameter': figure_parameter,
        'parameter_data': figure_parameter_data,
        'parameter_p_r': figure_parameter_p_r, 'parameter_r': figure_parameter_r,
        'committed_area': figure_committed_area,
        'map_pos': figure_map_pos,
        'data_extent': figure_data_extent,
        'grid_sampling': figure_grid_sampling, 'max_dist': figure_max_distance,
        'cmap': figure_cmap_obj,
        "vmin": figure_vmin, "vmax": figure_vmax,
        "clip": figure_clip,
        'lim_min': figure_lim_min, 'lim_max': figure_lim_max,
        'lim_thr': figure_lim_thr, 'lim_target': figure_lim_target, 'lim_optimal': figure_lim_optimal,
        'palette_type': figure_palette_type}

    return figure_file_obj, figure_committed_area_obj, settings_kwargs
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize committed area values
def organize_figure_committed_area(comm_area_values, comm_area_expected=None):

    if comm_area_expected is None:
        comm_area_expected = [True, False, None]

    if not isinstance(comm_area_values, list):
        comm_area_values = [comm_area_values]

    comm_area_verified = None
    for comm_area_step in comm_area_values:
        if comm_area_step in comm_area_expected:
            if comm_area_verified is None:
                comm_area_verified = []
            comm_area_verified.append(comm_area_step)
        else:
            logging.warning(' ===> Committed area type "' + str(comm_area_step) + '" is not supported by the algorithm')

    if comm_area_verified is None:
        logging.warning(' ===> Committed area type is defined only by NoneType')
        comm_area_verified = [comm_area_verified]

    return comm_area_verified
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to set figure type (committed_area or global)
def organize_figure_filename(file_name_tmpl='analysis_{data_type}.png', file_filter_committed_area=False):

    if not isinstance(file_filter_committed_area, list):
        file_filter_committed_area = [file_filter_committed_area]

    # iterate over committed area flags
    file_name_list = []
    for file_filter_carea in file_filter_committed_area:
        if file_filter_carea:
            file_name_filter = file_name_tmpl.format(data_type='committed_area')
        else:
            file_name_filter = file_name_tmpl.format(data_type='global')
        file_name_list.append(file_name_filter)

    # remove list duplicates
    file_name_list = list(set(file_name_list))

    # check filename and filter list
    if file_name_list.__len__() != file_filter_committed_area.__len__():
        logging.warning(' ===> Filename(s) defined by committed area flag are less then the expected string(s)')
        logging.warning(' ===> Equal filename(s) will be overwritten by the plotting functions')

    return file_name_list
# ----------------------------------------------------------------------------------------------------------------------
