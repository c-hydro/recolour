"""
Library Features:

Name:          lib_utils_model
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""


# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import json
import cartopy
import numpy as np
import pandas as pd

import geopandas as gpd
import contextily as ctx

from shapely.geometry import box

from lib_utils_generic import invert_dict
from lib_info_args import logger_name

import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.io.img_tiles as cimgt

from copy import deepcopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# logging
log_stream = logging.getLogger(logger_name)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter model data
def filter_model_data(
        dframe_data, dframe_fields=None,
        var_tag_rain='rain', var_tag_airt='air_temperature', var_tag_sm='soil_moisture',
        var_no_data=-9999.0,
        interp_limit_airt=2, interp_limit_sm=2,
        interp_direction_airt='both', interp_direction_sm='both',
        rain_min=0.0, rain_max=None, airt_min=-30.0, airt_max=60.0, sm_min=None, sm_max=None):

    # sort data by index
    dframe_data = dframe_data.sort_index().copy()

    # rain no data
    dframe_data.loc[dframe_data[var_tag_rain] == var_no_data, var_tag_rain] = np.nan

    # rain physical limits
    rain_mask = dframe_data[var_tag_rain] < rain_min
    if rain_max is not None:
        rain_mask |= dframe_data[var_tag_rain] > rain_max
    dframe_data.loc[rain_mask, var_tag_rain] = np.nan

    # remove rows where all values are nan
    dframe_data = dframe_data.dropna(how='all')

    # air temperature no data and limits
    mask_airt = (
        (dframe_data[var_tag_airt] == var_no_data) |
        (dframe_data[var_tag_airt] < airt_min) |
        (dframe_data[var_tag_airt] > airt_max)
    )
    dframe_data.loc[mask_airt, var_tag_airt] = np.nan

    # soil moisture: detect range if not provided
    sm_values = dframe_data[var_tag_sm].replace(var_no_data, np.nan).dropna()

    if sm_min is None:
        sm_min = 0.0
    if sm_max is None:
        if not sm_values.empty and sm_values.quantile(0.95) <= 1.5:
            sm_max = 1.0
        else:
            sm_max = 100.0

    # soil moisture no data and limits
    mask_sm = (
        (dframe_data[var_tag_sm] == var_no_data) |
        (dframe_data[var_tag_sm] < sm_min) |
        (dframe_data[var_tag_sm] > sm_max)
    )
    dframe_data.loc[mask_sm, var_tag_sm] = np.nan

    # fill nans with interpolation
    dframe_data[var_tag_airt] = dframe_data[var_tag_airt].interpolate(
        limit=interp_limit_airt,
        limit_direction=interp_direction_airt
    )

    dframe_data[var_tag_sm] = dframe_data[var_tag_sm].interpolate(
        limit=interp_limit_sm,
        limit_direction=interp_direction_sm
    )

    # dframe fields default
    if dframe_fields is None:
        dframe_fields = {}

    # organize file fields
    tmp_fields = invert_dict(dframe_fields)
    dframe_data = dframe_data.rename(columns=tmp_fields)

    return dframe_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter model results
def filter_model_results(
        dframe_data,
        var_tag_rain='rain',
        var_tag_airt='air_temperature',
        var_tag_theta_mod='theta_simulated',
        var_tag_theta_obs='theta_observed',
        var_no_data=-9999.0,
        rain_min=0.0, rain_max=None,
        airt_min=-30.0, airt_max=60.0,
        theta_mod_min=0.0, theta_mod_max=1.0,
        theta_obs_min=0.0, theta_obs_max=100.0):

    # sort data by index
    dframe_data = dframe_data.sort_index().copy()

    # rain no data and limits
    mask_rain = (
        (dframe_data[var_tag_rain] == var_no_data) |
        (dframe_data[var_tag_rain] < rain_min)
    )

    if rain_max is not None:
        mask_rain |= (dframe_data[var_tag_rain] > rain_max)

    dframe_data.loc[mask_rain, var_tag_rain] = np.nan

    # air temperature no data and limits
    mask_airt = (
        (dframe_data[var_tag_airt] == var_no_data) |
        (dframe_data[var_tag_airt] < airt_min) |
        (dframe_data[var_tag_airt] > airt_max)
    )

    dframe_data.loc[mask_airt, var_tag_airt] = np.nan

    # theta simulated no data and limits
    mask_theta_mod = (
        (dframe_data[var_tag_theta_mod] == var_no_data) |
        (dframe_data[var_tag_theta_mod] < theta_mod_min) |
        (dframe_data[var_tag_theta_mod] > theta_mod_max)
    )

    dframe_data.loc[mask_theta_mod, var_tag_theta_mod] = np.nan

    # theta observed no data and limits
    mask_theta_obs = (
        (dframe_data[var_tag_theta_obs] == var_no_data) |
        (dframe_data[var_tag_theta_obs] < theta_obs_min) |
        (dframe_data[var_tag_theta_obs] > theta_obs_max)
    )

    dframe_data.loc[mask_theta_obs, var_tag_theta_obs] = np.nan

    dframe_data = dframe_data.sort_index()

    return dframe_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize model results
def organize_model_results(dframe_common,
                           values_results, values_time, dframe_fields=None,
                           var_tag_time='time'):

    dict_results = {'values_model': values_results}
    dframe_results = pd.DataFrame(dict_results, index=values_time)
    dframe_results.index.name = var_tag_time

    dframe_common = dframe_common.join(dframe_results)

    # dframe fields default
    if dframe_fields is None:
        dframe_fields = {}
    # organize file fields
    dframe_common = dframe_common.rename(columns=dframe_fields)

    return dframe_common
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize model metrics
def organize_model_metrics(data_metrics, data_registry, data_time=None, data_fields=None):

    if data_time is None:
        data_time = {}

    # organize model metrics and registry
    data_common = {**data_metrics, **data_registry, **data_time}

    # define data fields
    if data_fields is None:
        data_list = list(data_common.keys())
    else:
        data_list = list(data_fields.keys())

    # iterate over data fields
    data_filter = {}
    for data_key, data_value in data_common.items():
        if data_key in data_list:
            data_filter[data_key] = data_value

    dframe_filter = pd.DataFrame(data=data_filter, index=['info'])

    return dframe_filter
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize model data
def organize_model_data(dframe_data, var_tag_rain='values_k1', var_tag_airt='values_k2', var_tag_sm='values_k3'):

    # get rain and time data
    if var_tag_rain in list(dframe_data.columns):
        dframe_data_rain = dframe_data[var_tag_rain]
        values_rain = dframe_data_rain.values
        values_time = dframe_data_rain.index
    else:
        log_stream.error(' ===> Variable "' + var_tag_rain + '" not included in the dataframe obj')
        raise IOError('Check your dataframe obj')
    # get airt data
    if var_tag_airt in list(dframe_data.columns):
        dframe_data_airt = dframe_data[var_tag_airt]
        values_airt = dframe_data_airt.values
    else:
        log_stream.error(' ===> Variable "' + var_tag_airt + '" not included in the dataframe obj')
        raise IOError('Check your dataframe obj')
    # get sm data
    if var_tag_sm in list(dframe_data.columns):
        dframe_data_sm = dframe_data[var_tag_sm]
        values_sm = dframe_data_sm.values
    else:
        log_stream.error(' ===> Variable "' + var_tag_sm + '" not included in the dataframe obj')
        raise IOError('Check your dataframe obj')

    if len(values_rain) != len(values_airt) != len(values_sm):
        log_stream.error(' ===> Variables have different length')
        raise IOError('Check your dataframe obj')

    values_data = np.dstack((values_rain, values_airt,values_sm))[0, :, :]

    return values_data, values_time
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize model parameters
def organize_model_parameters(params_dict, parameters_list=None, parameters_mandatory=True):

    if parameters_list is None:
        parameters_list = ['w_p', 'w_max', 'alpha', 'm2', 'ks', 'kc', 'theta_min', 'theta_max']

    parameters_values = []
    for parameters_name in parameters_list:

        if parameters_name in list(params_dict.keys()):
            parameters_values.append(params_dict[parameters_name])
        else:
            if parameters_mandatory:
                log_stream.error(' ===> Variable "' + parameters_name + '" not included in the parameters dictionary')
                raise IOError('Check your parameters dictionary')
            else:
                log_stream.warning(' ===> Variable "' + parameters_name + '" not included in the parameters dictionary')

    # convert parameters values in array
    values_params = np.array(parameters_values)

    return values_params
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to configure time-series axes
def configure_time_axes(time_data, time_format='%Y-%m-%d'):

    tick_time_period = list(time_data)
    tick_time_idx = time_data
    tick_time_labels = [tick_label.strftime(time_format) for tick_label in time_data]

    return tick_time_period, tick_time_idx, tick_time_labels
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot model results using time-series
def plot_model_results_ts(file_name, dframe_results, dframe_metrics,
                       fig_spacing_x=None,
                       fig_fields=None,
                       fig_dpi=150, fig_show=True, 
                       y_min_sm=0.1, y_max_sm=0.6,
                       no_data_rain=-9999.0, no_data_air_t=-9999.0,
                       no_data_theta_obs=-9999.0, no_data_theta_sim=-9999.0,
                       **kwargs):

    # get figure fields
    description_string = 'No description available'
    if fig_fields is not None:
        if 'description' in list(fig_fields.keys()):
            description_string = fig_fields['description']

    # sort data by index (time)
    dframe_results = dframe_results.sort_index()

    # compute expected time period and data frame
    time_index_start, time_index_end = dframe_results.index[0], dframe_results.index[-1]
    time_index_resolution = dframe_results.index.resolution

    if time_index_resolution == 'hour':
        time_index_frequency = 'h'
    else:
        log_stream.error(' ===> Time resolution not expected in the dataframe obj')
        raise NotImplemented('Case not implemented yet')
    time_index_range = pd.date_range(start=time_index_start, end=time_index_end, freq=time_index_frequency)

    dframe_expected = pd.DataFrame(index=time_index_range)
    dframe_expected = dframe_expected.join(dframe_results)

    # get ts values
    values_time = dframe_expected.index
    values_rain = deepcopy(dframe_expected['rain'].values)
    values_air_t = deepcopy(dframe_expected['air_temperature'].values)
    values_theta_obs = deepcopy(dframe_expected['theta_observed'].values)
    values_theta_sim = deepcopy(dframe_expected['theta_simulated'].values)

    # nullify no data values
    values_rain[values_rain == no_data_rain] = np.nan
    values_air_t[values_air_t == no_data_air_t] = np.nan
    values_theta_obs[values_theta_obs == no_data_theta_obs] = np.nan
    values_theta_sim[values_theta_sim == no_data_theta_sim] = np.nan

    # get metrics values
    metrics_ns = dframe_metrics['ns'].values[0]
    metrics_ns_ln_q = dframe_metrics['ns_ln_q'].values[0]
    metrics_ns_rad_q = dframe_metrics['ns_rad_q'].values[0]
    metrics_kge = dframe_metrics['kge'].values[0]
    metrics_rmse = dframe_metrics['rmse'].values[0]
    metrics_rq = dframe_metrics['rq'].values[0]
    # get registry values
    registry_name = dframe_metrics['name'].values[0].strip()
    registry_catchment = dframe_metrics['catchment'].values[0].strip()
    registry_time = dframe_metrics['time'].values[0]

    # filter ts values (remove no data values)
    values_rain[values_rain == -9999] = np.nan
    values_theta_obs[values_theta_obs == -9999] = np.nan
    values_theta_sim[values_theta_sim == -9999] = np.nan

    #y_min_sm = -0.05
    #y_max_sm = 1.05
    y_min_rain = 0
    y_max_rain = np.nanmax([np.nanmax(values_rain[np.isfinite(values_rain)]), 25])
    y_min_air_t = np.nanmin([np.nanmin(values_air_t[np.isfinite(values_air_t)]), -25])
    y_max_air_t = np.nanmax([np.nanmax(values_air_t[np.isfinite(values_air_t)]), 50])

    # select time start and time end
    time_start_string, time_start_stamp = values_time[0].strftime('%Y-%m-%d %H:%M'), values_time[0]
    time_end_string, time_end_stamp = values_time[-1].strftime('%Y-%m-%d %H:%M'), values_time[-1]

    time_period_tick, time_idx_tick, time_labels_tick = configure_time_axes(values_time)

    if fig_spacing_x is not None:

        spacing_type = 'automatic'
        if 'type' in list(fig_spacing_x.keys()):
            spacing_type = fig_spacing_x['type']
        spacing_offset = 0
        if 'offset' in list(fig_spacing_x.keys()):
            spacing_offset = fig_spacing_x['offset']

        if spacing_type == 'days':
            time_period_tick_start = time_period_tick[0]
            time_period_tick_end = time_period_tick[-1] + pd.DateOffset(days=spacing_offset)
        elif spacing_type == 'months':
            time_period_tick_start = time_period_tick[0]
            time_period_tick_end = time_period_tick[-1] + pd.DateOffset(months=spacing_offset)
        else:
            time_period_tick_start = time_period_tick[0]
            time_period_tick_end = time_period_tick[-1]
    else:
        time_period_tick_start = time_period_tick[0]
        time_period_tick_end = time_period_tick[-1]

    # plot figure
    fig = plt.figure(figsize=(10, 7))
    fig.autofmt_xdate()

    # title
    s = ('Point -- Name: "' + registry_name + '" Catchment: "' + registry_catchment +
         '" Description: "' + description_string + '"\n'
         ' Time Ref: "' + registry_time + ' Time Period Start: "' +
         time_start_string + '" Time Period End: "' + time_end_string +
         '" UTC \n'
         f'NS: "{metrics_ns:.3f}" NS(lnSD): "{metrics_ns_ln_q:.3f}" NS(radSD): "{metrics_ns_rad_q:.3f}" '
         f'RQ: "{metrics_rq:.3f}" RMSE: "{metrics_rmse:.3f}" KGE: "{metrics_kge:.3f}"')

    # upper panel (soil moisture)
    ax1 = plt.axes([0.1, 0.5, 0.8, 0.40])
    ax1.set_title(s, fontsize=10, fontweight='bold')
    ax1.plot(values_time, values_theta_obs, 'g', linewidth=2, label=r'$\theta_{obs}$')
    ax1.plot(values_time, values_theta_sim, 'r', linewidth=1, label=r'$\theta_{sim}$')
    ax1.legend()
    ax1.set_ylabel('Relative Soil Moisture [-]')
    ax1.set_xlim([time_period_tick_start, time_period_tick_end])
    ax1.set_ylim([y_min_sm, y_max_sm])
    ax1.tick_params(labelbottom=False)
    ax1.grid(visible=True)
    #ax1.grid(b=True)

    # lower panel (rain)
    ax2 = plt.axes([0.1, 0.1, 0.8, 0.40])
    ax2.plot(values_time, values_rain, color=[.5, .5, .5], linewidth=1, label='rain')
    ax2.set_ylabel('Rain (mm/h)')
    ax2.set_ylim([y_min_rain, y_max_rain])
    ax2.set_xlim(time_period_tick_start, time_period_tick_end)
    ax2.tick_params(axis='x', labelrotation=45, labelsize=6)
    ax2.grid(visible=True)
    #ax2.grid(b=True)

    ax3 = ax2.twinx()
    ax3.plot(values_time, values_air_t, color='r', linewidth=0.2, label='air temperature')
    ax3.set_ylabel('Air Temperature (C)')
    ax3.set_ylim(y_min_air_t, y_max_air_t)
    ax3.set_xlim(time_period_tick_start, time_period_tick_end)

    # save figure
    plt.savefig(file_name, format='png', dpi=fig_dpi)

    if fig_show:
        plt.show()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to plot model results using maps
def plot_model_results_maps(
        file_name, dframe_results,
        fig_fields=None,
        fig_dpi=150, fig_show=True,
        y_min_sm=0.0, y_max_sm=1.0,
        no_data_theta_sim=-9999.0,
        lon_name='longitude', lat_name='latitude',
        theta_name='theta_simulated',
        time_name='time',
        buffer_ratio=0.10,
        fig_background='osm',
        fig_color_map_type=None,
        fig_color_map_label='soil moisture',
        fig_color_tick_loc=None,
        fig_color_tick_label=None,
        fig_zoom=8,
        **kwargs):

    try:
        import lib_img_tiles as cimgt_custom
    except ImportError:
        cimgt_custom = None

    # -------------------------------------------------------------------------
    # Create destination folder
    folder_name = os.path.dirname(file_name)
    if folder_name:
        os.makedirs(folder_name, exist_ok=True)

    # -------------------------------------------------------------------------
    # Figure metadata
    title_string = 'Soil Moisture Map'
    subtitle_string = None
    time_reference = None

    if fig_fields is not None:
        title_string = fig_fields.get('title', title_string)
        subtitle_string = fig_fields.get('subtitle', subtitle_string)
        time_reference = fig_fields.get('time_reference', None)

    # Data time
    time_data = None

    if time_name in dframe_results.columns:
        time_values = pd.to_datetime(dframe_results[time_name].dropna().unique())
        if len(time_values) > 0:
            time_data = time_values[0]

    elif isinstance(dframe_results.index, pd.DatetimeIndex):
        time_data = dframe_results.index[0]

    time_reference = (
        pd.to_datetime(time_reference).strftime('%Y-%m-%d %H:%M')
        if time_reference is not None else 'NA'
    )

    time_data = (
        pd.to_datetime(time_data).strftime('%Y-%m-%d %H:%M')
        if time_data is not None else 'NA'
    )

    title_lines = [
        title_string,
        f'Time Ref: {time_reference} UTC | Time Data: {time_data} UTC'
    ]

    if subtitle_string is not None:
        title_lines.append(subtitle_string)

    full_title = '\n'.join(title_lines)

    # -------------------------------------------------------------------------
    # Clean dataframe
    dframe_plot = dframe_results.copy()

    dframe_plot[theta_name] = dframe_plot[theta_name].replace(
        no_data_theta_sim, np.nan
    )

    dframe_plot = dframe_plot.dropna(
        subset=[lon_name, lat_name, theta_name]
    )

    if dframe_plot.empty:
        raise RuntimeError('No valid points available to plot')

    # -------------------------------------------------------------------------
    # Create 2D grid from dataframe
    dframe_grid = dframe_plot.pivot_table(
        index=lat_name,
        columns=lon_name,
        values=theta_name,
        aggfunc='mean'
    )

    lat_values = dframe_grid.index.values.astype(np.float32)
    lon_values = dframe_grid.columns.values.astype(np.float32)
    map_data = dframe_grid.values.astype(np.float32)

    map_lons_2d, map_lats_2d = np.meshgrid(lon_values, lat_values)

    # Mask values outside limits
    if y_min_sm is not None:
        map_data[map_data < y_min_sm] = np.nan
    if y_max_sm is not None:
        map_data[map_data > y_max_sm] = np.nan

    # -------------------------------------------------------------------------
    # Bounds with buffer
    lon_min = np.nanmin(lon_values)
    lon_max = np.nanmax(lon_values)
    lat_min = np.nanmin(lat_values)
    lat_max = np.nanmax(lat_values)

    dx = np.nanmedian(np.diff(np.sort(lon_values))) if lon_values.size > 1 else 0.01
    dy = np.nanmedian(np.diff(np.sort(lat_values))) if lat_values.size > 1 else 0.01

    lon_buffer = max((lon_max - lon_min) * buffer_ratio, dx)
    lat_buffer = max((lat_max - lat_min) * buffer_ratio, dy)

    lon_min -= lon_buffer
    lon_max += lon_buffer
    lat_min -= lat_buffer
    lat_max += lat_buffer

    # -------------------------------------------------------------------------
    # Colormap
    if fig_color_map_type is None:
        sm_colors = [
            [0.58431372549, 0.266666666667, 0.0, 1.0],
            [0.709803921569, 0.494117647059, 0.0, 1.0],
            [0.894117647059, 0.827450980392, 0.0, 1.0],
            [0.992156862745, 1.0, 0.63137254902, 1.0],
            [0.650980392157, 0.933333333333, 0.96862745098, 1.0],
            [0.309803921569, 0.662745098039, 0.976470588235, 1.0],
            [0.160784313725, 0.333333333333, 0.658823529412, 1.0],
            [0.109803921569, 0.219607843137, 0.549019607843, 1.0]
        ]
        fig_color_map_obj = colors.LinearSegmentedColormap.from_list(
            'sm.cmap', sm_colors
        )

    elif isinstance(fig_color_map_type, colors.Colormap):
        fig_color_map_obj = deepcopy(fig_color_map_type)

    elif isinstance(fig_color_map_type, str):
        fig_color_map_obj = plt.get_cmap(fig_color_map_type)

    else:
        raise NotImplementedError('Variable colormap case not implemented yet')

    if fig_color_tick_loc is None:
        fig_color_tick_loc = np.linspace(y_min_sm, y_max_sm, 6)

    if fig_color_tick_label is None:
        fig_color_tick_label = [f'{v:.2f}' for v in fig_color_tick_loc]

    # -------------------------------------------------------------------------
    # CRS
    plot_crs = cartopy.crs.Mercator()
    data_crs = cartopy.crs.PlateCarree()

    # Background
    if fig_background == 'osm':
        if cimgt_custom is not None:
            map_background = cimgt_custom.OSM()
        else:
            map_background = cimgt.OSM()

    elif fig_background == 'google':
        map_background = cimgt.GoogleTiles()

    else:
        raise NotImplementedError('Variable background case not implemented yet')

    # -------------------------------------------------------------------------
    # Figure
    fig = plt.figure(figsize=(12, 10))

    fig.suptitle(
        full_title,
        fontsize=12,
        color='black',
        fontweight='bold',
        y=0.97
    )

    ax = fig.add_axes([0.1, 0.10, 0.75, 0.75], projection=plot_crs)


    ax.set_extent(
        [lon_min, lon_max, lat_min, lat_max],
        crs=data_crs
    )

    gl = ax.gridlines(
        crs=data_crs,
        draw_labels=True,
        linewidth=2,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )

    # New API
    if hasattr(gl, 'top_labels'):
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = True
    # Old API
    else:
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlabels_bottom = True
        gl.ylabels_left = True

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color': 'gray', 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'color': 'gray', 'weight': 'bold'}

    ax.add_image(map_background, fig_zoom)

    sc = ax.pcolormesh(
        map_lons_2d,
        map_lats_2d,
        map_data,
        zorder=3,
        cmap=fig_color_map_obj,
        transform=data_crs,
        vmin=y_min_sm,
        vmax=y_max_sm
    )

    # Colorbar
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)

    cb = plt.colorbar(sc, cax=ax_cb, extend='both')
    cb.set_label(
        f'{fig_color_map_label} [-]',
        size=12,
        color='gray',
        weight='normal'
    )
    cb.ax.tick_params(labelsize=10, labelcolor='gray')

    if fig_color_tick_loc is not None:
        cb.set_ticks(fig_color_tick_loc)

    if fig_color_tick_label is not None:
        cb.set_ticklabels(fig_color_tick_label)

    fig.savefig(file_name, dpi=fig_dpi)

    if fig_show:
        plt.show()
    else:
        plt.close(fig)
# ----------------------------------------------------------------------------------------------------------------------
