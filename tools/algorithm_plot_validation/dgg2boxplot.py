# -----------------------------------------------------------------------------
# Libraries
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

from os.path import join

import numpy as np
import seaborn as sns
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Method to plot boxplot pearson
def plot_boxplot_pearson(c_df, filename,
                         data_field='type',
                         title=None, suptitle=None,
                         val='ASCAT_GLDAS_ALL_R', p_val='ASCAT_GLDAS_ALL_p_R',
                         lim_threshold=0.5, lim_target=0.65, lim_optimal=0.8,
                         lim_min=-1, lim_max=1,
                         ylabel=None,
                         xlabel=None,
                         pal='Set2'):

    if np.sum(c_df.columns.isin([val])) > 0:

        c_count = c_df.groupby(data_field).count()

        sub_df = c_df.loc[:, [val, p_val, data_field]].dropna()
        sub_df.loc[(sub_df[p_val] > 0.05), val] = np.nan
        sub_count2 = sub_df.dropna().groupby(data_field).count()

        optimal_df = c_df.loc[:, [val, p_val, data_field]].dropna()
        optimal_df.loc[(optimal_df[val] < lim_optimal), val] = np.nan
        optimal_count = optimal_df.dropna().groupby(data_field).count()

        target_df = c_df.loc[:, [val, p_val, data_field]].dropna()
        target_df.loc[(target_df[val] < lim_target), val] = np.nan
        target_count = target_df.dropna().groupby(data_field).count()

        threshold_df = c_df.loc[:, [val, p_val, data_field]].dropna()
        threshold_df.loc[(threshold_df[val] < lim_threshold), val] = np.nan
        threshold_count = threshold_df.dropna().groupby(data_field).count()

        fig, ax = plt.subplots(1, sharey=True, figsize=(3, 6))

        flierprops = dict(markerfacecolor='0.75', markersize=5,
                          linestyle='none')

        sns.boxplot(y=val, x=data_field, data=sub_df,
                    palette=pal,
                    saturation=1., linewidth=1.5, width=0.2,
                    ax=ax,
                    flierprops=flierprops, showfliers=False, whis=[5, 95], )

        xticklabels = []
        percentages = {}
        for label in ax.get_xticklabels():
            optimal_cnt = optimal_count.loc[label.get_text(), val]
            target_cnt = target_count.loc[label.get_text(), val]
            threshold_cnt = threshold_count.loc[label.get_text(), val]

            c_cnt = c_count.loc[label.get_text(), val]

            percentages[label] = {}
            optimal_perc = int((float(optimal_cnt) / float(c_cnt)) * 100)
            target_perc = int((float(target_cnt) / float(c_cnt)) * 100)
            threshold_perc = int((float(threshold_cnt) / float(c_cnt)) * 100)

            percentages[label]['optimal'] = '{:}%'.format(optimal_perc)
            percentages[label]['target'] = '{:}%'.format(target_perc)
            percentages[label]['threshold'] = '{:}%'.format(threshold_perc)

            xticklabels.append('{:}'.format(label.get_text()))

        if xlabel is not None:
            ax.set_xticklabels(xlabel, fontsize=13)
        else:
            ax.set_xticklabels(xticklabels, fontsize=13)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.set_xlabel('')

        if title is not None:
            ax.set_title(title)

        if suptitle is not None:
            fig.suptitle(suptitle)

        ax.set_ylim([lim_min, lim_max])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        plt.axhline(y=lim_optimal, linewidth=2, color='b', linestyle='-', label='optimal')
        plt.axhline(y=lim_target, linewidth=2, color='g', linestyle='-', label='target')
        plt.axhline(y=lim_threshold, linewidth=2, color='r', linestyle='-', label='threshold')

        for idx, label in enumerate(percentages):
            label_perc = percentages[label]

            plt.text(idx + 0.25, lim_optimal + 0.02, label_perc['optimal'], color='b', fontsize=8)
            plt.text(idx + 0.25, lim_target + 0.02, label_perc['target'], color='g', fontsize=8)
            plt.text(idx + 0.25, lim_threshold + 0.02, label_perc['threshold'], color='r', fontsize=8)

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels,
                        ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1),
                        prop={'size': 6})
        lgd.get_frame().set_linewidth(0.0)
        ax.grid('off')

        if filename is not None:
            fig.savefig(join(filename),
                        dpi=150, bbox_extra_artists=(lgd,), bbox_inches='tight')

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Method to plot boxplot snr
def plot_boxplot_snr(c_df, filename,
                     data_field='type',
                     title=None, suptitle=None,
                     val='ASCAT_ALL_snr', p_val='ASCAT_GLDAS_ALL_p_R', r_val='ASCAT_GLDAS_ALL_R',
                     lim_threshold=0, lim_target=3, lim_optimal=6,
                     lim_min=-15, lim_max=15,
                     ylabel=None, xlabel=None,
                     pal='Set2'):

    if np.sum(c_df.columns.isin([val])) > 0:
        c_count = c_df.groupby(data_field).count()

        sub_df = c_df.loc[:, [val, p_val, r_val, data_field]].dropna()
        # sub_df.loc[(sub_df[p_val] > 0.05) | (sub_df[r_val] < 0.5), val] = np.nan
        # sub_count = sub_df.dropna().groupby(data_field).count()
        # sub_df.loc[(sub_df['lat'] > 55) & (sub_df[val] < 6.0), val] = np.nan

        optimal_df = c_df.loc[:, [val, p_val, r_val, data_field]].dropna()
        optimal_df.loc[(optimal_df[val] < lim_optimal), val] = np.nan
        optimal_count = optimal_df.dropna().groupby(data_field).count()

        target_df = c_df.loc[:, [val, p_val, r_val, data_field]].dropna()
        target_df.loc[(target_df[val] < lim_target), val] = np.nan
        target_count = target_df.dropna().groupby(data_field).count()

        threshold_df = c_df.loc[:, [val, p_val, r_val, data_field]].dropna()
        threshold_df.loc[(threshold_df[val] < lim_threshold), val] = np.nan
        threshold_count = threshold_df.dropna().groupby(data_field).count()

        flierprops = dict(markerfacecolor='0.75', markersize=5,
                          linestyle='none')

        fig, ax = plt.subplots(1, sharey=True, figsize=(3, 6))
        sns.boxplot(y=val, x=data_field, data=sub_df,
                    palette=pal,
                    saturation=1., linewidth=1.5, width=0.2,
                    ax=ax,
                    flierprops=flierprops, showfliers=False,  whis=[5, 95],)

        xticklabels = []
        percentages = {}
        for label in ax.get_xticklabels():
            optimal_cnt = optimal_count.loc[label.get_text(), val]
            target_cnt = target_count.loc[label.get_text(), val]
            threshold_cnt = threshold_count.loc[label.get_text(), val]

            c_cnt = c_count.loc[label.get_text(), val]

            percentages[label] = {}
            optimal_perc = int((float(optimal_cnt)/float(c_cnt)) * 100)
            target_perc = int((float(target_cnt) / float(c_cnt)) * 100)
            threshold_perc = int((float(threshold_cnt) / float(c_cnt)) * 100)

            percentages[label]['optimal'] = '{:}%'.format(optimal_perc)
            percentages[label]['target'] = '{:}%'.format(target_perc)
            percentages[label]['threshold'] = '{:}%'.format(threshold_perc)

            xticklabels.append('{:}'.format(label.get_text()))

        if xlabel is not None:
            ax.set_xticklabels(xlabel, fontsize=13)
        else:
            ax.set_xticklabels(xticklabels, fontsize=13)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.set_xlabel('')

        if title is not None:
            ax.set_title(title)

        if suptitle is not None:
            fig.suptitle(suptitle)

        ax.set_ylim([lim_min, lim_max])

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        plt.axhline(y=lim_optimal, linewidth=2, color='b', linestyle='-', label='optimal')
        plt.axhline(y=lim_target, linewidth=2, color='g', linestyle='-', label='target')
        plt.axhline(y=lim_threshold, linewidth=2, color='r', linestyle='-', label='threshold')

        for idx, label in enumerate(percentages):
            label_perc = percentages[label]

            plt.text(idx + 0.25, lim_optimal+0.3, label_perc['optimal'], color='b', fontsize=8)
            plt.text(idx + 0.25, lim_target + 0.3, label_perc['target'], color='g', fontsize=8)
            plt.text(idx + 0.25, lim_threshold + 0.3, label_perc['threshold'], color='r', fontsize=8)

        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels,
                        ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1),
                        prop={'size': 6})
        lgd.get_frame().set_linewidth(0.0)
        ax.grid('off')

        if filename is not None:
            fig.savefig(join(filename),
                        dpi=150, bbox_extra_artists=(lgd,), bbox_inches='tight')

# -----------------------------------------------------------------------------
