"""
Library Features:

Name:          lib_notebook_figure
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
from copy import deepcopy
from datetime import date

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot time series by year
def plot_ts_by_year(ts_dframe, tag_ref='sm_obs', tag_kn=None, tag_year=None,
                    figure_file_name=None,
                    figure_x_label='time', figure_y_label='sm [%]',
                    figure_title='Soil Moisture Time Series', figure_point='point',
                    figure_marker=None,
                    figure_legend=None, figure_color=None,
                    figure_ylim=None, figure_size=(14, 10), figure_dpi=150):

    if tag_year is None:
        tag_year = []
    if figure_ylim is None:
        figure_ylim = [0, 100]
    if figure_color is None:
        figure_color = ['black', 'red', 'blue']
    if figure_legend is None:
        figure_legend = ['obs', 'k1', 'k2']
    if figure_marker is None:
        figure_marker = ['o', None, None]
    if tag_kn is None:
        tag_kn = ['sm_smap_exp_t02', 'sm_smap_exp_t02']

    year_n = len(tag_year)

    if year_n == 0:
        return

    # open figure
    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)
    fig.autofmt_xdate()

    # iterate over years
    for year_id, year_tag in enumerate(tag_year):

        ax = plt.subplot(year_n, 1, year_id + 1)

        figure_label = [tag_ref] + tag_kn

        # check if year is available in the dataset
        if year_tag in ts_dframe['year'].values:

            year_dframe = ts_dframe[ts_dframe['year'] == year_tag]
            select_dframe = year_dframe[figure_label]

            year_start = date(year_tag, 1, 1).strftime('%Y-%m-%d')
            year_end = date(year_tag, 12, 31).strftime('%Y-%m-%d')
            year_range = pd.date_range(start=year_start, end=year_end, freq='D')

            null_data = np.zeros(shape=(len(year_range), len(figure_label)))
            null_data[:, :] = np.nan
            null_dict = {}
            for ts_id, ts_label in enumerate(figure_label):
                null_dict[ts_label] = null_data[:, ts_id]
            base_dframe = pd.DataFrame(data=null_dict, index=year_range)

            base_dframe.update(select_dframe)
            vars_dframe = deepcopy(base_dframe)

            print()

        else:

            year_start = date(year_tag, 1, 1).strftime('%Y-%m-%d')
            year_end = date(year_tag, 12, 31).strftime('%Y-%m-%d')
            year_range = pd.date_range(start=year_start, end=year_end, freq='D')

            null_data = np.zeros(shape=(len(year_range), len(figure_label)))
            null_data[:, :] = np.nan
            null_dict = {}
            for ts_id, ts_label in enumerate(figure_label):
                null_dict[ts_label] = null_data[:, ts_id]

            vars_dframe = pd.DataFrame(data=null_dict, index=year_range)

            print(' ===> Year ' + str(year_tag) + ' not available in the dataset')

        for ts_label, ts_marker, ts_color in zip(figure_label, figure_marker, figure_color):
            if ts_marker is not None:
                vars_dframe[ts_label].plot(figsize=figure_size, ax=ax, marker=ts_marker, ms=1, color=ts_color)
            else:
                vars_dframe[ts_label].plot(figsize=figure_size, ax=ax, color=ts_color)



        ax.set_xlabel(figure_x_label)
        ax.set_ylabel(figure_y_label)
        ax.grid(True)

        ax.legend(figure_legend)

        plt.ylim(figure_ylim)

        if year_id == 0:
            ts_title_label = figure_title + ' -- ' + figure_point
            plt.title(ts_title_label)

    fig.tight_layout()
    plt.show()

    # save figure
    if figure_file_name is not None:
        file_path, file_folder = os.path.split(figure_file_name)
        os.makedirs(file_path, exist_ok=True)
        fig.savefig(figure_file_name, dpi=figure_dpi)

    plt.close()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot generic time series
def plot_ts_generic(ts_dframe, tag_ref='sm_obs', tag_kn=None,
                    figure_file_name=None,
                    figure_x_label='time', figure_y_label='sm [%]',
                    figure_title='Soil Moisture Time Series', figure_point='point',
                    figure_marker=None,
                    figure_legend=None, figure_color=None,
                    figure_ylim=None, figure_size=(10, 5), figure_dpi=150):

    if figure_ylim is None:
        figure_ylim = [0, 100]
    if figure_color is None:
        figure_color = ['black', 'red', 'blue']
    if figure_legend is None:
        figure_legend = ['obs', 'k1', 'k2']
    if figure_marker is None:
        figure_marker = ['o', None, None]
    if tag_kn is None:
        tag_kn = ['sm_smap_exp_t02', 'sm_smap_exp_t02']

    ts_title_label = figure_title + ' (' + figure_point + ')'

    fig, ax = plt.subplots()

    figure_label = [tag_ref] + tag_kn

    for ts_label, ts_marker, ts_color in zip(figure_label, figure_marker, figure_color):
        if ts_marker is not None:
            ts_dframe[ts_label].plot(figsize=figure_size, ax=ax, marker=ts_marker, color=ts_color)
        else:
            ts_dframe[ts_label].plot(figsize=figure_size, ax=ax, color=ts_color)

    ax.set_xlabel(figure_x_label)
    ax.set_ylabel(figure_y_label)
    ax.grid(True)

    ax.legend(figure_legend)

    plt.ylim(figure_ylim)
    plt.title(ts_title_label)
    plt.show()

    if figure_file_name is not None:
        fig.savefig(figure_file_name, dpi=figure_dpi)
# ----------------------------------------------------------------------------------------------------------------------
