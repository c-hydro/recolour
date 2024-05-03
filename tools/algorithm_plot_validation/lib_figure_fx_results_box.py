"""
Library Features:

Name:          lib_figure_fx_results_box
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240411'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import seaborn as sns

import matplotlib.pylab as plt

from lib_figure_fx_tools import get_field
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results box snr
def plot_results_box_snr(fig_file_name,
                         fig_df_obj, fig_perc_obj, fig_kwargs,
                         variable_data="x_snr", variable_p_r="xy_p_r", variable_r="xy_r",
                         variable_type='type'):

    # sub dataframe obj
    fig_df_sub = fig_df_obj.loc[:, [variable_data, variable_p_r, variable_r, variable_type]].dropna()

    # active filter(s)
    # fig_df_sub = fig_df_sub.loc[fig_df_sub[variable_r] > 0.05]
    # df_obj.loc[(df_obj[variable_p_r] > 0.05) | (df_obj[variable_r] < 0.5), variable_data] = np.nan

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    fig_size = get_field(fig_kwargs, field_name='fig_size')
    fig_x_label = get_field(fig_kwargs, field_name='x_label')
    fig_y_label = get_field(fig_kwargs, field_name='y_label')
    fig_palette_type = get_field(fig_kwargs, field_name='palette_type')
    fig_lim_min = get_field(fig_kwargs, field_name='lim_min')
    fig_lim_max = get_field(fig_kwargs, field_name='lim_max')
    fig_lim_thr = get_field(fig_kwargs, field_name='lim_thr')
    fig_lim_target = get_field(fig_kwargs, field_name='lim_target')
    fig_lim_optimal = get_field(fig_kwargs, field_name='lim_optimal')
    fig_title = get_field(fig_kwargs, field_name='title')
    fig_suptitle = get_field(fig_kwargs, field_name='suptitle', field_mandatory=False, field_default=None)
    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')

    fig, ax = plt.subplots(1, sharey=True, figsize=(3, 6))
    flierprops = dict(markerfacecolor='0.75', markersize=5, linestyle='none')

    sns.boxplot(
        data=fig_df_sub, y=variable_data, x=variable_type,
        palette=fig_palette_type,
        saturation=1., linewidth=1.5, width=0.2,
        ax=ax,
        flierprops=flierprops, showfliers=False, whis=[5, 95], )

    xticklabels = []
    for label in ax.get_xticklabels():
        xticklabels.append('{:}'.format(label.get_text()))

    if fig_x_label is not None:
        ax.set_xticklabels(fig_x_label, fontsize=13)
    else:
        ax.set_xticklabels(xticklabels, fontsize=13)
    ax.set_ylabel(fig_y_label, fontsize=13)
    ax.set_xlabel('')

    if fig_title is not None:
        ax.set_title(fig_title)
    if fig_suptitle is not None:
        fig.suptitle(fig_suptitle)

    ax.set_ylim([fig_lim_min, fig_lim_max])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    plt.axhline(y=fig_lim_optimal, linewidth=2, color='b', linestyle='-', label='optimal')
    plt.axhline(y=fig_lim_target, linewidth=2, color='g', linestyle='-', label='target')
    plt.axhline(y=fig_lim_thr, linewidth=2, color='r', linestyle='-', label='threshold')

    for idx, label in enumerate(fig_perc_obj):
        label_perc = fig_perc_obj[label]

        plt.text(idx + 0.25, fig_lim_optimal + 0.25, label_perc['optimal'], color='b', fontsize=8)
        plt.text(idx + 0.25, fig_lim_target + 0.25, label_perc['target'], color='g', fontsize=8)
        plt.text(idx + 0.25, fig_lim_thr + 0.25, label_perc['threshold'], color='r', fontsize=8)

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(
        handles, labels,
        ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1), prop={'size': 6})
    lgd.get_frame().set_linewidth(0.0)
    ax.grid('off')

    fig.savefig(fig_file_name, dpi=fig_dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results box pearson
def plot_results_box_pearson(fig_file_name,
                             fig_df_obj, fig_perc_obj, fig_kwargs,
                             variable_data="xy_r", variable_p_r="xy_p_r", variable_r="xy_r",
                             variable_type='type'):

    # sub dataframe obj
    if variable_data == variable_r:
        fig_df_sub = fig_df_obj.loc[:, [variable_data, variable_p_r, variable_type]].dropna()
    else:
        fig_df_sub = fig_df_obj.loc[:, [variable_data, variable_p_r, variable_r, variable_type]].dropna()

    # active filter(s)
    # fig_df_sub = fig_df_sub.loc[fig_df_sub[variable_p_r] > 0.05]
    # df_obj = df_obj.loc[(df_obj[variable_p_r] > 0.05) | (df_obj[variable_r] < 0.5), variable_data] = np.nan

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    fig_size = get_field(fig_kwargs, field_name='fig_size')
    fig_x_label = get_field(fig_kwargs, field_name='x_label')
    fig_y_label = get_field(fig_kwargs, field_name='y_label')
    fig_palette_type = get_field(fig_kwargs, field_name='palette_type')
    fig_lim_min = get_field(fig_kwargs, field_name='lim_min')
    fig_lim_max = get_field(fig_kwargs, field_name='lim_max')
    fig_lim_thr = get_field(fig_kwargs, field_name='lim_thr')
    fig_lim_target = get_field(fig_kwargs, field_name='lim_target')
    fig_lim_optimal = get_field(fig_kwargs, field_name='lim_optimal')
    fig_title = get_field(fig_kwargs, field_name='title')
    fig_suptitle = get_field(fig_kwargs, field_name='suptitle', field_mandatory=False, field_default=None)
    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')

    fig, ax = plt.subplots(1, sharey=True, figsize=(3, 6))
    flierprops = dict(markerfacecolor='0.75', markersize=5, linestyle='none')

    sns.boxplot(
        data=fig_df_sub, y=variable_data, x=variable_type,
        palette=fig_palette_type,
        saturation=1., linewidth=1.5, width=0.2,
        ax=ax,
        flierprops=flierprops, showfliers=False, whis=[5, 95], )

    xticklabels = []
    for label in ax.get_xticklabels():
        xticklabels.append('{:}'.format(label.get_text()))

    if fig_x_label is not None:
        ax.set_xticklabels(fig_x_label, fontsize=13)
    else:
        ax.set_xticklabels(xticklabels, fontsize=13)
    ax.set_ylabel(fig_y_label, fontsize=13)
    ax.set_xlabel('')

    if fig_title is not None:
        ax.set_title(fig_title)
    if fig_suptitle is not None:
        fig.suptitle(fig_suptitle)

    ax.set_ylim([fig_lim_min, fig_lim_max])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    plt.axhline(y=fig_lim_optimal, linewidth=2, color='b', linestyle='-', label='optimal')
    plt.axhline(y=fig_lim_target, linewidth=2, color='g', linestyle='-', label='target')
    plt.axhline(y=fig_lim_thr, linewidth=2, color='r', linestyle='-', label='threshold')

    for idx, label in enumerate(fig_perc_obj):
        label_perc = fig_perc_obj[label]

        plt.text(idx + 0.25, fig_lim_optimal + 0.02, label_perc['optimal'], color='b', fontsize=8)
        plt.text(idx + 0.25, fig_lim_target + 0.02, label_perc['target'], color='g', fontsize=8)
        plt.text(idx + 0.25, fig_lim_thr + 0.02, label_perc['threshold'], color='r', fontsize=8)

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(
        handles, labels,
        ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1), prop={'size': 6})

    lgd.get_frame().set_linewidth(0.0)
    ax.grid('off')

    plt.show()

    fig.savefig(fig_file_name, dpi=fig_dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')

# ----------------------------------------------------------------------------------------------------------------------
