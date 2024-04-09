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
def configure_time_series_heatmap(point_ts, point_legend=None):

    point_label_list, point_name_list, point_arr = [], [], None
    for point_name, point_data in point_ts.items():
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
def view_time_series(file_name, point_ts,
                     point_metrics, point_registry, min_ts=0, max_ts=100,
                     file_fields=None, file_groups=None, file_metrics=None,
                     fig_title='figure', fig_label_axis_x='time', fig_label_axis_y='value',
                     fig_legend=None, fig_style=None,
                     fig_cbar='coolwarm', fig_cbar_kw={}, fig_dpi=120):

    # ------------------------------------------------------------------------------------------------------------------
    # figure registry
    registry_tag, registry_name, registry_code = point_registry['tag'], point_registry['name'], point_registry['code']
    registry_geo_x, registry_geo_y = point_registry['longitude'], point_registry['latitude']

    point_ts[point_ts < min_ts] = np.nan
    point_ts[point_ts > max_ts] = np.nan

    # figure fields
    if 'time' in point_ts.index.name:
        time_period = point_ts.index
        time_start, time_end = time_period[0], time_period[-1]
    else:
        log_stream.error(' ===> Time field not in the dataframe object. Time field must be included in dataframe')
        raise RuntimeError('Time field not in the dataframe object')

    # configure time-series axes
    [tick_time_period, tick_time_idx, tick_time_labels] = configure_time_series_axes(point_ts)
    # configure time-series data for heatmaps
    point_arr, point_name, point_label = configure_time_series_heatmap(point_ts, fig_legend)
    # configure time-series data for heatmaps
    metrics_arr, metrics_ts, metrics_type = configure_time_series_metrics(point_metrics)

    point_arr[point_arr < min_ts] = np.nan
    point_arr[point_arr > max_ts] = np.nan

    metrics_arr[metrics_arr == -9999] = np.nan
    metrics_arr[metrics_arr == -9998] = np.nan

    # open figure
    fig = plt.figure(figsize=(17, 10))
    fig.autofmt_xdate()

    # super title
    fig_title = fig_title.format(point_name=registry_name, point_tag=registry_tag, point_code=registry_code)
    plt.suptitle(fig_title, size=12, color='black', weight='bold')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # subplot 1 (heatmap)
    ax1 = plt.subplot(3, 9, (1, 7))
    ax1.set_xticklabels([])

    # create heatmap
    image_norm = mpl.colors.Normalize(vmin=min_ts, vmax=max_ts)
    image_renderer = ax1.imshow(point_arr, cmap=fig_cbar, norm=image_norm)

    # create colorbar
    cbar = ax1.figure.colorbar(image_renderer, ax=ax1, shrink=0.6, **fig_cbar_kw)
    cbar.ax.set_ylabel(fig_label_axis_y, rotation=-90, va="bottom", color='#000000')

    # show all ticks and label them with the respective list entries
    ax1.set_xticks(np.arange(len(tick_time_labels)))
    ax1.set_xticklabels(tick_time_labels)
    ax1.set_yticks(np.arange(len(point_label)))
    ax1.set_yticklabels(point_label)

    # rotate the tick labels and set their alignment.
    plt.setp(ax1.get_xticklabels(), rotation=90, fontsize=6, ha="right", rotation_mode="anchor")
    # create grid
    ax1.set_xticks(np.arange(point_arr.shape[1] + 1) - .5, minor=True)
    ax1.set_yticks(np.arange(point_arr.shape[0] + 1) - .5, minor=True)
    ax1.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax1.tick_params(which="minor", bottom=False, left=False)

    # create annotations
    for i in range(len(point_label)):
        for j in range(len(tick_time_period)):
            text = ax1.text(
                j, i, point_arr[i, j].round(1),
                ha="center", va="center", color="w", fontsize=6, fontweight='bold')
    # set title
    # ax1.set_title(fig_title, size=12, color='black', weight='bold')
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
                          rowLabels=metrics_rows,
                          colLabels=metrics_cols,
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

    # iterate over serie(s)
    px_obj, label_obj = [], []
    for point_id, (point_key, point_data) in enumerate(point_ts.items()):

        # get style and label
        point_style = configure_time_series_style(point_key, fig_style)
        point_label = configure_time_series_lut(point_key, fig_legend)

        if point_style is None:
            point_px = ax3.plot(point_data.index, point_data.values,
                                label=point_label, lw=1)
        else:
            point_px = ax3.plot(point_data.index, point_data.values,
                                label=point_label, **point_style)

        px_obj.append(point_px[0])
        label_obj.append(point_label)

    px_obj, label_obj = tuple(px_obj), tuple(label_obj)

    ax3.set_xlabel(fig_label_axis_x, color='#000000')
    ax3.set_xlim(tick_time_period[0], tick_time_period[-1])
    # ax3.set_ylabel(fig_label_axis_y, color='#000000')
    ax3.set_ylim(min_ts, max_ts)
    ax3.grid(b=True)

    ax3.set_xticks(tick_time_idx)
    ax3.set_xticklabels(tick_time_labels, rotation=90, fontsize=6)

    ax4 = ax3.twinx()
    ax4.set_ylabel(fig_label_axis_y, rotation=-90, va="bottom", color='#000000')
    ax4.set_ylim(min_ts, max_ts)

    # add legend
    legend = ax3.legend(px_obj, label_obj, frameon=False, ncol=2, loc=0)
    ax3.add_artist(legend)

    fig.tight_layout()

    # save figure
    file_path, file_folder = os.path.split(file_name)
    os.makedirs(file_path, exist_ok=True)
    fig.savefig(file_name, dpi=fig_dpi)

    # close figure
    # plt.show()
    plt.close()

# ----------------------------------------------------------------------------------------------------------------------
