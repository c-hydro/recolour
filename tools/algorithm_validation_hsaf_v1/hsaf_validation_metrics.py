# Copyright (c) 2017, Vienna University of Technology (TU Wien), Department
# of Geodesy and Geoinformation (GEO).
# All rights reserved.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY,
# DEPARTMENT OF GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
This module defines the metrics used for the validation.
"""
import copy
import numpy as np
import itertools

import pytesmo.metrics as metrics
import pytesmo.df_metrics as df_metrics


class BasicMetrics(object):

    """
    This class computes basic metrics. It also stores information about
    gpi, lat, lon and number of observations.
    """

    def __init__(self, result_path=None, other_name='k1'):

        self.result_path = result_path
        self.other_name = other_name

        metrics = {'R': np.float32([np.nan]),
                   'p_R': np.float32([np.nan]),
                   'rho': np.float32([np.nan]),
                   'p_rho': np.float32([np.nan]),
                   'n_obs': np.int32([0])}

        self.result_template = {'gpi': np.int32([-1]),
                                'lon': np.float32([np.nan]),
                                'lat': np.float32([np.nan])}

        self.seasons = ['ALL', 'DJF', 'MAM', 'JJA', 'SON']

        for season in self.seasons:
            for metric in metrics.keys():
                key = "{:}_{:}".format(season, metric)
                self.result_template[key] = metrics[metric].copy()

        self.month_to_season = np.array(['', 'DJF', 'DJF', 'MAM', 'MAM',
                                         'MAM', 'JJA', 'JJA', 'JJA', 'SON',
                                         'SON', 'SON', 'DJF'])

    def calc_metrics(self, data, gpi_info):
        """
        Method to calculate the desired statistics

        Parameters
        ----------
        data : pandas.DataFrame
            with 2 columns, the first column is the reference dataset
            named 'ref'
            the second column the dataset to compare against named 'other'
        gpi_info : tuple
            Grid point info (i.e. gpi, lon, lat)
        """
        dataset = copy.deepcopy(self.result_template)

        dataset['gpi'][0] = gpi_info[0]
        dataset['lon'][0] = gpi_info[1]
        dataset['lat'][0] = gpi_info[2]

        for season in self.seasons:

            if season != 'ALL':
                subset = self.month_to_season[data.index.month] == season
            else:
                subset = np.ones(len(data), dtype=bool)

            if subset.sum() < 10:
                continue

            x = data['ref'].values[subset]
            y = data[self.other_name].values[subset]
            R, p_R = metrics.pearsonr(x, y)
            rho, p_rho = metrics.spearmanr(x, y)

            dataset['{:}_n_obs'.format(season)][0] = subset.sum()
            dataset['{:}_R'.format(season)][0] = R
            dataset['{:}_p_R'.format(season)][0] = p_R
            dataset['{:}_rho'.format(season)][0] = rho
            dataset['{:}_p_rho'.format(season)][0] = p_rho

        return dataset


class HSAF_Metrics(object):

    """
    This class computes metrics as defined by the H-SAF consortium in
    order to prove the operational readiness of a product. It also stores
    information about gpi, lat, lon and number of observations.
    """

    def __init__(self, other_name1='k1', other_name2='k2',
                 dataset_names=None, seasonal_metrics=False):

        self.other_name1 = other_name1
        self.other_name2 = other_name2
        self.df_columns = ['ref', self.other_name1, self.other_name2]

        if dataset_names is None:
            self.ds_names = self.df_columns
        else:
            self.ds_names = dataset_names

        # create lut between df columns and dataset names
        self.ds_names_lut = {}
        for name, col in zip(self.ds_names, self.df_columns):
            self.ds_names_lut[col] = name

        self.tds_names = []
        for combi in itertools.combinations(self.df_columns, 2):
            self.tds_names.append("{:}_and_{:}".format(*combi))

        metrics_common = {'n_obs': np.int32([0])}

        metrics_sds = {'snr': np.float32([np.nan]),
                       'err_var': np.float32([np.nan]),
                       'beta': np.float32([np.nan])}

        metrics_tds = {'R': np.float32([np.nan]),
                       'p_R': np.float32([np.nan]),
                       'rho': np.float32([np.nan]),
                       'p_rho': np.float32([np.nan]),
                       'bias': np.float32([np.nan]),
                       'ubrmsd': np.float32([np.nan])}

        self.result_template = {'gpi': np.int32([-1]),
                                'lon': np.float32([np.nan]),
                                'lat': np.float32([np.nan])}

        if seasonal_metrics:
            self.seasons = ['ALL', 'DJF', 'MAM', 'JJA', 'SON']
        else:
            self.seasons = ['ALL']

        for season in self.seasons:
            # get template for common metric
            for metric in metrics_common.keys():
                key = "{:}_{:}".format(season, metric)
                self.result_template[key] = metrics_common[metric].copy()

            # get template for single-dataset metric
            for name in self.ds_names:
                for metric in metrics_sds.keys():
                    key = "{:}_{:}_{:}".format(name, season, metric)
                    self.result_template[key] = metrics_sds[metric].copy()

            # get template for two-dataset metric
            for tds_name in self.tds_names:
                split_tds_name = tds_name.split('_and_')
                tds_name_key = "{:}_{:}".format(self.ds_names_lut[
                    split_tds_name[0]],
                    self.ds_names_lut[
                    split_tds_name[1]])
                for metric in metrics_tds.keys():
                    key = "{:}_{:}_{:}".format(tds_name_key, season, metric)
                    self.result_template[key] = metrics_tds[metric].copy()

        self.month_to_season = np.array(['', 'DJF', 'DJF', 'MAM', 'MAM',
                                         'MAM', 'JJA', 'JJA', 'JJA', 'SON',
                                         'SON', 'SON', 'DJF'])

    def calc_metrics(self, data, gpi_info):
        """
        Method to calculate the desired statistics

        Parameters
        ----------
        data : pandas.DataFrame
            with 3 columns, the first column is the reference dataset
            named 'ref'
            the second and third column are the datasets to compare against
            named 'k1 and k2'
        gpi_info : tuple
            Grid point info (i.e. gpi, lon, lat)
        """
        dataset = copy.deepcopy(self.result_template)

        dataset['gpi'][0] = gpi_info[0]
        dataset['lon'][0] = gpi_info[1]
        dataset['lat'][0] = gpi_info[2]

        for season in self.seasons:

            if season != 'ALL':
                subset = self.month_to_season[data.index.month] == season
            else:
                subset = np.ones(len(data), dtype=bool)

            # number of observations
            n_obs = subset.sum()
            if n_obs < 10:
                continue

            dataset['{:}_n_obs'.format(season)][0] = n_obs

            # get single dataset metrics
            # calculate SNR
            x = data[self.df_columns[0]].values[subset]
            y = data[self.df_columns[1]].values[subset]
            z = data[self.df_columns[2]].values[subset]

            snr, err, beta = metrics.tcol_snr(x, y, z)

            for i, name in enumerate(self.ds_names):
                dataset['{:}_{:}_snr'.format(name, season)][0] = snr[i]
                dataset['{:}_{:}_err_var'.format(name, season)][0] = err[i]
                dataset['{:}_{:}_beta'.format(name, season)][0] = beta[i]

            # calculate Pearson correlation
            pearson_R, pearson_p = df_metrics.pearsonr(data)
            pearson_R = pearson_R._asdict()
            pearson_p = pearson_p._asdict()

            # calculate Spearman correlation
            spea_rho, spea_p = df_metrics.pearsonr(data)
            spea_rho = spea_rho._asdict()
            spea_p = spea_p._asdict()

            # calculate bias
            bias_nT = df_metrics.bias(data)
            bias_dict = bias_nT._asdict()

            # calculate ubRMSD
            ubRMSD_nT = df_metrics.ubrmsd(data)
            ubRMSD_dict = ubRMSD_nT._asdict()

            for tds_name in self.tds_names:
                R = pearson_R[tds_name]
                p_R = pearson_p[tds_name]
                rho = spea_rho[tds_name]
                p_rho = spea_p[tds_name]
                bias = bias_dict[tds_name]
                ubRMSD = ubRMSD_dict[tds_name]

                split_tds_name = tds_name.split('_and_')
                tds_name_key = "{:}_{:}".format(self.ds_names_lut[
                    split_tds_name[0]],
                    self.ds_names_lut[
                    split_tds_name[1]])

                dataset['{:}_{:}_R'.format(tds_name_key, season)][0] = R
                dataset['{:}_{:}_p_R'.format(tds_name_key, season)][0] = p_R
                dataset['{:}_{:}_rho'.format(tds_name_key, season)][0] = rho
                dataset['{:}_{:}_p_rho'.format(tds_name_key, season)][0] = \
                    p_rho
                dataset['{:}_{:}_bias'.format(tds_name_key, season)][0] = bias
                dataset['{:}_{:}_ubrmsd'.format(tds_name_key, season)][0] = \
                    ubRMSD

        return dataset
