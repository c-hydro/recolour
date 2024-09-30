"""
Library Features:

Name:          lib_utils_graph
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240109'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import os

import numpy as np
import pandas as pd

import seaborn as sns

import matplotlib.dates as mdates
import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib.font_manager import FontProperties

from copy import deepcopy

from lib_info_args import logger_name

# logging
warnings.filterwarnings('ignore')
logging.getLogger("matplotlib").setLevel(logging.WARNING)
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure time series data
def configure_time_series_info(ts_data_in, ts_metrics_in, fields=None, metrics=None):

    if fields is not None:

        ts_data_list, ts_metrics_list = [], []
        for ts_name_in, ts_name_out in fields.items():
            if ts_name_in != 'time':
                if ts_name_in in ts_data_in.columns:
                    ts_data_list.append(ts_name_in)
                else:
                    log_stream.warning(' ===> Field "' + ts_name_in +
                                       '" not in the data dataframe object. Field keep the standard name')
                if ts_name_in in ts_metrics_in.columns:
                    ts_metrics_list.append(ts_name_in)
                else:
                    log_stream.warning(' ===> Field "' + ts_name_in +
                                       '" not in the metrics dataframe object. Field keep the standard name')

    ts_metrics_list = sorted(ts_metrics_list)
    ts_data_list = sorted(ts_data_list)

    if set(ts_metrics_list) != set(ts_data_list):
        ts_common_list = set(ts_data_list).intersection(ts_metrics_list)
        ts_data_out = deepcopy(ts_data_in[ts_common_list])
        ts_metrics_out = deepcopy(ts_metrics_in[ts_common_list])
    else:
        ts_data_out = deepcopy(ts_data_in[ts_data_list])
        ts_metrics_out = deepcopy(ts_metrics_in[ts_data_list])

    if metrics is not None:
        ts_metrics_filtered = {}
        for ts_metrics_key, ts_metrics_fields in ts_metrics_out.items():
            tmp_metrics_fields = {k: v for k, v in ts_metrics_fields.items() if k in metrics}
            ts_metrics_filtered[ts_metrics_key] = tmp_metrics_fields
        ts_metrics_out = deepcopy(ts_metrics_filtered)

    return ts_data_out, ts_metrics_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to configure time-series axes
def configure_time_series_axes(dframe_data, time_format='%m-%d %H'):

    tick_time_period = list(dframe_data.index)
    tick_time_idx = dframe_data.index
    tick_time_labels = [tick_label.strftime(time_format) for tick_label in dframe_data.index]

    return tick_time_period, tick_time_idx, tick_time_labels
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure time-series lut
def configure_time_series_lut(ts_name_in, ts_fields=None):

    ts_name_out = None
    if ts_fields is not None:
        for ts_field_key_in, ts_field_key_out in ts_fields.items():
            if ts_name_in == ts_field_key_in:
                ts_name_out = ts_field_key_out
                break
    else:
        log_stream.error(' ===> Time-series fields must be defined in the configuration file')
        raise RuntimeError('Time-series fields not defined')

    if ts_name_out is None:
        log_stream.error(' ===> Time-series name "' + ts_name_in + '" is not defined in the time-series fields')
        raise RuntimeError('Time-series name not defined')

    return ts_name_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure time-series style
def configure_time_series_style(ts_name, obj_configuration=None):
        ts_configuration = None
        if obj_configuration is not None:
            if ts_name in list(obj_configuration.keys()):
                ts_configuration = obj_configuration[ts_name]
        return ts_configuration
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure time-series for heatmaps
def configure_time_series_heatmap(point_ts, point_columns_excluded=None, point_legend=None):

    if point_columns_excluded is None:
        point_columns_excluded = ['time']

    if not isinstance(point_columns_excluded, list):
        point_columns_excluded = [point_columns_excluded]

    point_label_list, point_name_list, point_arr = [], [], None
    for point_name, point_data in point_ts.items():

        if not point_name in point_columns_excluded:
            if point_arr is None:
                point_arr = point_data.values
            else:
                point_arr = np.vstack([point_arr, point_data.values])

            if point_legend is not None:
                point_label = configure_time_series_lut(point_name, point_legend)
            else:
                point_label = deepcopy(point_name)

            point_label_list.append(point_label)
            point_name_list.append(point_name)

    return point_arr, point_name_list, point_label_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure metrics for table
def configure_time_series_metrics(metrics_data):

    metrics_type_list, metrics_name_list, metrics_arr = None, [], None
    for metric_name, metric_data in metrics_data.items():

        if metrics_type_list is None:
            metrics_type_list = list(metric_data.keys())

        tmp_arr = np.array(list(metric_data.values()))

        if metrics_arr is None:
            metrics_arr = tmp_arr
        else:
            metrics_arr = np.vstack([metrics_arr, tmp_arr])

        metrics_name_list.append(metric_name)

    return metrics_arr, metrics_name_list, metrics_type_list
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to view time series
def view_time_series(file_name, point_ts_data,
                     point_metrics, point_registry, min_ts=0, max_ts=100,
                     file_fields=None, file_groups_name=None, file_groups_time=None, file_metrics=None,
                     fig_title='figure', fig_label_axis_x='time', fig_label_axis_y='value',
                     fig_legend=None, fig_style=None,
                     fig_spacing_x=None,
                     fig_cbar='coolwarm', fig_dpi=120):

    # ------------------------------------------------------------------------------------------------------------------
    # figure style
    fig_style_ts, fig_style_hm = None, None
    if 'time_series' in list(fig_style.keys()):
        fig_style_ts = fig_style['time_series']
        for fig_key, fig_value in fig_legend.items():
            if fig_key in fig_style_ts:
                fig_style_ts[fig_value] = fig_style_ts[fig_key]
                fig_style_ts.pop(fig_key)

    if 'heatmap' in list(fig_style.keys()):
        fig_style_hm = fig_style['heatmap']
    else:
        fig_style_hm = {
            "line_width": 0.2, "line_color": "black", "line_weight": "bold",
            "cmap": fig_cbar, "cbar_label": "variable [-]"}
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # figure registry
    registry_tag, registry_name, registry_code = point_registry['tag'], point_registry['name'], point_registry['code']
    registry_geo_x, registry_geo_y = point_registry['longitude'], point_registry['latitude']

    point_ts_data[point_ts_data < min_ts] = np.nan
    point_ts_data[point_ts_data > max_ts] = np.nan

    # figure fields
    if 'time' in point_ts_data.index.name:
        time_period = point_ts_data.index
        time_start, time_end = time_period[0], time_period[-1]
    else:
        log_stream.error(' ===> Time field not in the dataframe object. Time field must be included in dataframe')
        raise RuntimeError('Time field not in the dataframe object')

    if file_groups_time is not None:
        from pandas.tseries.offsets import DateOffset
        time_sub_period = file_groups_time['time_sub_period']
        time_sub_n = file_groups_time['time_sub_name']

        time_sub_period_start = deepcopy(time_start)
        time_sub_period_end = deepcopy(time_end)

        time_period_collections, time_period_id = {}, 1
        time_sub_period_tmp = deepcopy(time_sub_period_start)
        while time_sub_period_tmp < time_sub_period_end:
            time_sub_period_lower = deepcopy(time_sub_period_tmp)
            time_sub_period_tmp = time_sub_period_tmp + DateOffset(**time_sub_period)
            time_sub_period_upper = deepcopy(time_sub_period_tmp)

            time_period_collections[time_period_id] = [time_sub_period_lower, time_sub_period_upper]

            time_period_id += 1

        time_period_max = np.max(list(time_period_collections.keys()))

        if time_period_collections[time_period_max][1] > time_sub_period_end:
            time_period_collections[time_period_max][1] = time_sub_period_end

    else:
        time_sub_n = '{:02d}'
        time_period_collections, time_period_id = {}, 0
        time_period_collections[time_period_id] = [time_start, time_end]

    # iterate over period(s)
    for time_period_id, time_period_range in time_period_collections.items():

        # select time start and time end
        time_start_string, time_start_stamp = time_period_range[0].strftime('%Y%m%d%H%M'), time_period_range[0]
        time_end_string, time_end_stamp = time_period_range[1].strftime('%Y%m%d%H%M'), time_period_range[1]

        # get time-series data for selected period
        point_ts_string = time_sub_n.format(time_period_id)
        point_ts_period = point_ts_data.loc[time_start_stamp:time_end_stamp]

        # adjust time period information  (for heatmap)
        point_ts_period['day'] = pd.to_datetime(point_ts_period.index).date
        point_ts_period = point_ts_period.rename(columns=fig_legend)
        pivot_ts_period = point_ts_period.pivot_table(columns='day', aggfunc='mean')
        pivot_ts_period.index.name = None

        point_ts_period = point_ts_period.drop(['day'], axis=1)

        time_range = pd.date_range(start=time_start_stamp, end=time_end_stamp, freq='12H')

        day_space = 1
        if fig_spacing_x is not None:
            if fig_spacing_x['type'] == 'days':
                if 'period' in fig_spacing_x:
                    day_space = fig_spacing_x['period']
        day_list_generic = []
        for time_step in time_range:
            if time_step.day not in day_list_generic:
                day_list_generic.append(time_step.day)
        day_list_ticks = day_list_generic[0::day_space]
        if time_range[-1].day != day_list_ticks[-1]:
            day_list_ticks.append(time_range[-1].day)

        # define file name
        if ('time_start' not in file_name) and ('time_end' not in file_name):
            file_name_def = file_name.format(
                group_name=file_groups_name, time_sub_name=point_ts_string)
        elif ('time_start' in file_name) and ('time_end' in file_name):
            file_name_def = file_name.format(
                group_name=file_groups_name, time_sub_name=point_ts_string,
                time_start=time_start_string, time_end=time_end_string)
        elif ('time_start' in file_name) and ('time_end' not in file_name):
            file_name_def = file_name.format(
                group_name=file_groups_name, time_sub_name=point_ts_string, time_start=time_start_string)
        elif ('time_start' not in file_name) and ('time_end' in file_name):
            file_name_def = file_name.format(
                group_name=file_groups_name, time_sub_name=point_ts_string, time_end=time_end_string)
        else:
            log_stream.error(' ===> File name not correctly defined')
            raise RuntimeError('File name not correctly defined')

        # configure time-series axes
        [tick_time_period, tick_time_idx, tick_time_labels] = configure_time_series_axes(point_ts_period)
        # configure time-series data for heatmaps
        point_arr, point_name, point_label = configure_time_series_heatmap(
            point_ts_period, point_columns_excluded=['time', 'date'], point_legend=None)
        # configure time-series data for heatmaps
        metrics_arr, metrics_ts, metrics_type = configure_time_series_metrics(point_metrics)

        point_arr[point_arr < min_ts] = np.nan
        point_arr[point_arr > max_ts] = np.nan

        metrics_arr[metrics_arr == -9999] = np.nan
        metrics_arr[metrics_arr == -9998] = np.nan

        # open figure
        fig = plt.figure(figsize=(18, 10))
        fig.autofmt_xdate()

        # super title
        fig_title = fig_title.format(
            point_name=registry_name, point_tag=registry_tag, point_code=registry_code,
            time_start=time_start_string, time_end=time_end_string, group_name=file_groups_name)
        plt.suptitle(fig_title, size=12, color='black', weight='bold')
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # subplot 1 (heatmap)
        ax1 = plt.subplot(3, 9, (1, 7))
        ax1.set_xticklabels([])

        # create heatmap
        ax1 = sns.heatmap(data=pivot_ts_period, vmin=min_ts, vmax=max_ts,
                          annot=True, fmt='.1f', center=None,
                          linewidths=fig_style_hm['line_width'],
                          linecolor=fig_style_hm['line_color'],
                          cbar=True, square=True,
                          annot_kws={
                              'size': fig_style_hm['text_size'],
                              'color': fig_style_hm['text_color'],
                              "weight": fig_style_hm['text_weight']
                          },
                          cbar_kws={
                              "shrink": 0.5, 'aspect': 10,
                              'label': fig_style_hm['cbar_label']
                          },
                          ax=ax1, cmap=fig_style_hm['cmap'])

        # rotate the tick labels and set their alignment
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, fontsize=6) # ha="right", rotation_mode="anchor"
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # subplot 2 (metrics)
        ax2 = plt.subplot(3, 9, (8, 9))

        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.axis("off")

        metrics_text, metrics_rows, metrics_cols = [], [], metrics_type
        for point_id, metrics_name in enumerate(metrics_ts):
            metrics_row = metrics_arr[point_id, :]
            metrics_text.append(['%1.2f' % x for x in metrics_row])

            metrics_row = configure_time_series_lut(metrics_name, fig_legend)
            metrics_rows.append(metrics_row)

        # add a table in the center
        table = plt.table(cellText=metrics_text,
                          rowLabels=metrics_rows, colLabels=metrics_cols,
                          loc='center')

        for (row, col), cell in table.get_celld().items():
            if (row == 0) or (col == -1):
                cell.set_text_props(fontproperties=FontProperties(weight='bold'))

        # Adjust layout to make room for the table:
        plt.subplots_adjust(left=0.2, bottom=0.2)
        # ------------------------------------------------------------------------------------------------------------------

        # ------------------------------------------------------------------------------------------------------------------
        # subplot 3 (series)
        ax3 = plt.subplot(3, 9, (10, 27))

        # iterate over series
        px_obj, label_obj = [], []
        for point_id, (point_key, point_data) in enumerate(point_ts_period.items()):

            if point_key not in ['time', 'day']:

                # get style and label
                point_style = configure_time_series_style(point_key, fig_style_ts)
                # point_label = configure_time_series_lut(point_key, fig_legend)

                if point_style is None:
                    point_px = ax3.plot(point_data.index, point_data.values,
                                        label=point_key, lw=1)
                else:
                    point_px = ax3.plot(point_data.index, point_data.values,
                                        label=point_key, **point_style)

                px_obj.append(point_px[0])
                label_obj.append(point_key)

        px_obj, label_obj = tuple(px_obj), tuple(label_obj)

        ax3.set_xlabel(fig_label_axis_x, color='#000000')
        ax3.set_xlim(tick_time_period[0], tick_time_period[-1])
        # ax3.set_ylabel(fig_label_axis_y, color='#000000')
        ax3.set_ylim(min_ts, max_ts)
        ax3.grid(b=True)

        if fig_spacing_x is None:
            ax3.set_xticks(tick_time_idx)
            ax3.set_xticklabels(tick_time_labels, rotation=90, fontsize=6)
        elif fig_spacing_x['type'] == 'days':
            ax3.xaxis.set_major_locator(mdates.DayLocator(day_list_ticks))
            ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d %H"))
            ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0, fontsize=6)
        elif fig_spacing_x['type'] == 'automatic':
            pass
        else:
            log_stream.error(' ===> Spacing type not correctly defined')
            raise RuntimeError('Spacing type not correctly defined')

        ax4 = ax3.twinx()
        ax4.set_ylabel(fig_label_axis_y, rotation=-90, va="bottom", color='#000000')
        ax4.set_ylim(min_ts, max_ts)

        # add legend
        legend = ax3.legend(px_obj, label_obj, frameon=False, ncol=2, loc=0)
        ax3.add_artist(legend)

        fig.tight_layout()

        # save figure
        file_folder_def, _ = os.path.split(file_name_def)
        os.makedirs(file_folder_def, exist_ok=True)
        fig.savefig(file_name_def, dpi=fig_dpi)

        # close figure
        # plt.show()
        plt.close()

# ----------------------------------------------------------------------------------------------------------------------
