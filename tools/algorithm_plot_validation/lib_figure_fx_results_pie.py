"""
Library Features:

Name:          lib_figure_fx_results_pie
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240411'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import matplotlib.pylab as plt

from lib_figure_fx_tools import get_field
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results pie pearson
def plot_results_pie_pearson(fig_file_name, fig_perc_comm, fig_perc_global, fig_kwargs):

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    # settings figure
    fig_committed_area = get_field(fig_kwargs, field_name='fig_committed_area')
    fig_size = get_field(fig_kwargs, field_name='fig_size')
    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')

    # open and plot figure
    if fig_committed_area:
        fig, ax = plt.subplots(1, 2, figsize=fig_size)
    else:
        fig, ax = plt.subplots(1, 1, figsize=fig_size)
        ax = [ax]   # just to pass the same structure

    # cmap and colors
    size = 0.3
    cmap = get_field(fig_kwargs, field_name='cmap')
    if isinstance(cmap, list):
        comm_colors = cmap[0]
        global_colors = cmap[1]
    else:
        cmap = plt.get_cmap('tab20c')
        comm_colors = cmap([4, 2, 1, 0, 16])
        global_colors = cmap([13, 11, 10, 9, 16])

    # global case
    wedges, texts = ax[0].pie(
        fig_perc_global, radius=0.6, colors=global_colors,
        wedgeprops=dict(width=size, edgecolor='w'))

    ax[0].set(aspect="equal")

    labels = [
        '{:3.1f}%: -1 $\leq$ R < 0.5'.format(fig_perc_global[0]),
        '{:3.1f}%: 0.5 $\leq$ R < 0.65'.format(fig_perc_global[1]),
        '{:3.1f}%: 0.65 $\leq$ R < 0.8'.format(fig_perc_global[2]),
        '{:3.1f}%: 0.8 $\leq$ R $\leq$ 1'.format(fig_perc_global[3]),
        '{:3.1f}%: NaN'.format(fig_perc_global[4])]

    ax[0].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                 [labels[3], labels[2], labels[1], labels[0], labels[4]],
                 title="Global", loc="center",
                 fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    # committed case
    if fig_committed_area:
        wedges, texts = ax[1].pie(
            fig_perc_comm, radius=0.6, colors=comm_colors,
            wedgeprops=dict(width=size, edgecolor='w'))

        ax[1].set(aspect="equal")

        labels = [
            '{:3.1f}%: -1 $\leq$ R < 0.5'.format(fig_perc_comm[0]),
            '{:3.1f}%: 0.5 $\leq$ R < 0.65'.format(fig_perc_comm[1]),
            '{:3.1f}%: 0.65 $\leq$ R < 0.8'.format(fig_perc_comm[2]),
            '{:3.1f}%: 0.8 $\leq$ R $\leq$ 1'.format(fig_perc_comm[3]),
            '{:3.1f}%: NaN'.format(fig_perc_comm[4])]

        ax[1].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                     [labels[3], labels[2], labels[1], labels[0], labels[4]],
                     title="Committed area", loc="center",
                     fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    # plt.subplots_adjust(hspace=0.5, left=0.214, bottom=0, right=0.98, top=0.78, wspace=0.512)
    plt.subplots_adjust(hspace=0.5, left=0.28, bottom=0, right=1, top=0.78, wspace=0.614)

    # save figure
    fig.savefig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot results pie snr
def plot_results_pie_snr(fig_file_name, fig_perc_comm, fig_perc_global, fig_kwargs):

    # define folder to save figure
    folder_name, file_name = os.path.split(fig_file_name)
    os.makedirs(folder_name, exist_ok=True)

    # settings figure
    fig_dpi = get_field(fig_kwargs, field_name='fig_dpi')
    fig_committed_area = get_field(fig_kwargs, field_name='fig_committed_area')
    fig_size = get_field(fig_kwargs, field_name='fig_size')

    # open and plot figure
    if fig_committed_area:
        fig, ax = plt.subplots(1, 2, figsize=fig_size)
    else:
        fig, ax = plt.subplots(1, 1, figsize=fig_size)
        ax = [ax]   # just to pass the same structure

    # cmap and colors
    size = 0.3
    cmap = get_field(fig_kwargs, field_name='cmap')
    if isinstance(cmap, list):
        comm_colors = cmap[0]
        global_colors = cmap[1]
    else:
        cmap = plt.get_cmap('tab20c')
        comm_colors = cmap([4, 2, 1, 0, 16])
        global_colors = cmap([13, 11, 10, 9, 16])

    # global case
    wedges, texts = ax[0].pie(
        fig_perc_global, radius=0.6, colors=global_colors,
        wedgeprops=dict(width=size, edgecolor='w'))

    ax[0].set(aspect="equal")

    labels = [
        '{:3.1f}%: SNR < 0 dB'.format(fig_perc_global[0]),
        '{:3.1f}%: 0 dB $\leq$ SNR < 3 dB'.format(fig_perc_global[1]),
        '{:3.1f}%: 3 dB $\leq$ SNR < 6 dB'.format(fig_perc_global[2]),
        '{:3.1f}%: SNR $\geq$ 6 dB'.format(fig_perc_global[3]),
        '{:3.1f}%: NaN'.format(fig_perc_global[4])]

    ax[0].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                 [labels[3], labels[2], labels[1], labels[0], labels[4]],
                 title="Global", loc="center",
                 fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    # committed case
    if fig_committed_area:
        wedges, texts = ax[1].pie(
            fig_perc_comm, radius=0.6, colors=comm_colors,
            wedgeprops=dict(width=size, edgecolor='w'))

        ax[1].set(aspect="equal")

        labels = [
            '{:3.1f}%: SNR < 0 dB'.format(fig_perc_comm[0]),
            '{:3.1f}%: 0 dB $\leq$ SNR < 3 dB'.format(fig_perc_comm[1]),
            '{:3.1f}%: 3 dB $\leq$ SNR < 6 dB'.format(fig_perc_comm[2]),
            '{:3.1f}%: SNR $\geq$ 6 dB'.format(fig_perc_comm[3]),
            '{:3.1f}%: NaN'.format(fig_perc_comm[4])]

        ax[1].legend([wedges[3], wedges[2], wedges[1], wedges[0], wedges[4]],
                     [labels[3], labels[2], labels[1], labels[0], labels[4]],
                     title="Committed area", loc="center",
                     fontsize=7, title_fontsize=8, bbox_to_anchor=(-0.3, 0.5))

    # plt.subplots_adjust(hspace=0.5, left=0.214, bottom=0, right=0.98, top=0.78, wspace=0.512)
    plt.subplots_adjust(hspace=0.5, left=0.28, bottom=0, right=1, top=0.78, wspace=0.614)

    # save figure
    fig.savefig(fig_file_name, dpi=fig_dpi)

# ----------------------------------------------------------------------------------------------------------------------
