"""
Library Features:

Name:          lib_utils_metrics
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import copy
import logging
import collections
import numpy as np
import itertools

import pandas as pd
import hydroeval as he
import pytesmo.df_metrics as df_metrics
import pytesmo.validation_framework.error_handling as eh
from numpy.lib.recfunctions import unstructured_to_structured
from pyswi.swi_ts.swi_ts import calc_swi_ts
from scipy.optimize import curve_fit, brute

import matplotlib.pyplot as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to initialize basic metrics
class BasicMetrics(object):
    """
    Compute basic metrics and store information about gpi, lat, lon, and the number of observations.

    Parameters
    ----------
    result_path : str, optional
        Path to store the results (default: None).
    other_name : str, optional
        A string representing the name of the dataset to compare against (default: 'k1').

    Attributes
    ----------
    result_path : str
        The path to store the results.
    other_name : str
        The name of the dataset to compare against.
    result_template : dict
        A template dictionary for storing computed metrics.
    seasons : list of str
        A list of seasons for which metrics are computed.
    month_to_season : numpy.ndarray
        A mapping from month index to season name.

    Methods
    -------
    calc_metrics(data, gpi_info)
        Calculate the desired statistics.

        Parameters
        ----------
        data : pandas.DataFrame
            A pandas DataFrame containing the reference and other datasets.
        gpi_info : tuple
            Grid point information (gpi, lon, lat).

        Returns
        -------
        dict
            A dictionary containing computed metrics for the provided data and grid point.
    """
    def __init__(self, result_path=None, other_name='k1',
                 dataset_names=None, seasonal_metrics=False, swi_option=False):

        self.result_path = result_path
        self.other_name = other_name
        self.df_columns = ['ref', self.other_name]
        self.swi_option = swi_option

        if self.swi_option:
            logging.warning(f' -------> SMAP missing, so SWI calculation not implemented.')

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


        # metrics' dictionaries
        metrics_common = {'n_obs': np.int32([0])}

        metrics_sds = {
            'nonan_n_obs': np.float32([np.nan]),
            'perc_nan': np.float32([np.nan]),
            'var': np.float32([np.nan]),
            'mean': np.float32([np.nan]),
            'max': np.float32([np.nan]),
            'min': np.float32([np.nan]),
        }

        metrics_tds = {'R': np.float32([np.nan]),
                       'p_R': np.float32([np.nan]),
                       'rho': np.float32([np.nan]),
                       'p_rho': np.float32([np.nan]),
                       'bias': np.float32([np.nan]),
                       'ubrmsd': np.float32([np.nan])}

        self.result_template = {'gpi': np.int32([-1]),
                                'lon': np.float32([np.nan]),
                                'lat': np.float32([np.nan]),
                                "status": np.int32([eh.UNCAUGHT])
                                }

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
                tds_name_key = "{:}_{:}".format(
                    self.ds_names_lut[
                        split_tds_name[0]],
                    self.ds_names_lut[
                        split_tds_name[1]])
                for metric in metrics_tds.keys():
                    key = "{:}_{:}_{:}".format(tds_name_key, season, metric)
                    self.result_template[key] = metrics_tds[metric].copy()

        self.month_to_season = np.array(['', 'DJF', 'DJF', 'MAM', 'MAM',
                                         'MAM', 'JJA', 'JJA', 'JJA', 'SON',
                                         'SON', 'SON', 'DJF'])
    # ------------------------------------------------------------------------

    def calc_metrics(self, data, gpi_info):
        """Method to calculate the desired statistics.

        Parameters
        ----------
        data : pandas.DataFrame
            A pandas DataFrame containing the reference and other datasets.
        gpi_info : tuple
            Grid point information (gpi, lon, lat).

        Returns
        -------
        dict
            A dictionary containing computed metrics for the provided data and grid point.
        """
        dataset = copy.deepcopy(self.result_template)

        dataset['gpi'][0] = gpi_info[0]
        dataset['lon'][0] = gpi_info[1]
        dataset['lat'][0] = gpi_info[2]

        # ensure data is properly sorted
        data = data[['ref', self.other_name]]

        for season in self.seasons:

            if season != 'ALL':
                subset = self.month_to_season[data.index.month] == season
            else:
                subset = np.ones(len(data), dtype=bool)

            data = data[subset]

            # calculate single-product metrics
            x = data.ref.values;
            y = data.k1.values
            data_vectors = [x, y]

            # check data usability
            skip = False
            if self.check_data(data_vectors, self.ds_names) == 0: skip = True
            if skip: continue

            nnan_obs = len(data) - data.isna().sum().values
            perc_nans= data.isna().sum().values / len(data)
            nanvars  = [np.nanvar(e)  for e in data_vectors]
            nanmeans = [np.nanmean(e) for e in data_vectors]
            nanmaxs  = [np.nanmax(e)  for e in data_vectors]
            nanmins  = [np.nanmin(e)  for e in data_vectors]


            for i, name in enumerate(self.ds_names):
                dataset['{:}_{:}_nonan_n_obs'.format(name, season)][0] = nnan_obs[i]
                dataset['{:}_{:}_perc_nan'.format(name, season)][0] = perc_nans[i]
                dataset['{:}_{:}_var'.format(name, season)][0] = nanvars[i]
                dataset['{:}_{:}_mean'.format(name, season)][0] = nanmeans[i]
                dataset['{:}_{:}_max'.format(name, season)][0] = nanmaxs[i]
                dataset['{:}_{:}_min'.format(name, season)][0] = nanmins[i]

            asdict = lambda x : x._asdict()
            # calculate pairs' metrics
            # Pearson's correlation coefficient (R) and
            # Spearman's (rho) as well as their respective
            # two-tailed p-values (p_R, p_rho)
            R, p_R = map(asdict, df_metrics.pearsonr(data))
            spea_rho, spea_p = map(asdict, df_metrics.spearmanr(data))
            bias_nT = asdict(df_metrics.bias(data))
            ubRMSD_nT = asdict(df_metrics.ubrmsd(data))

            for tds_name in self.tds_names:
                R = R[tds_name]
                p_R = p_R[tds_name]
                rho = spea_rho[tds_name]
                p_rho = spea_p[tds_name]
                bias = bias_nT[tds_name]
                ubRMSD = ubRMSD_nT[tds_name]

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

        dataset["status"][0] = eh.OK
        return dataset

    # ------------------------------------------------------------------------

    def check_data(self, data_vectors, names):
        """Performs checks to determine data usability.

        Calls same-name static method from ExtendedMetrics.

        Parameters
        ----------
        data_vectors : list of lists, ndarray
            List or array of vectors of data.
        names : list of strings
            List of the names of the columns in data_vectors, ordered.

        Returns
        -------
        exit_values : int
            Is 0 if data has issues (which are dumped to the log), 1 if all checks are passed

        """
        exit_value = ExtendedMetrics.check_data(data_vectors, names)
        return exit_value


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to initialize extended metrics
class ExtendedMetrics(object):
    """
    Compute metrics as defined by the H-SAF consortium to prove the operational readiness of a product.
    This class also stores information about gpi, lat, lon, and the number of observations.

    Parameters
    ----------
    other_name1 : str, optional
        A string representing the name of the first additional dataset (default: 'k1').
    other_name2 : str, optional
        A string representing the name of the second additional dataset (default: 'k2').
    dataset_names : list of str, optional
        A list of dataset names to use (default: None, uses ['ref', other_name1, other_name2']).
    seasonal_metrics : bool, optional
        A boolean indicating whether to compute seasonal metrics (default: False).

    Attributes
    ----------
    other_name1 : str
        The name of the first additional dataset.
    other_name2 : str
        The name of the second additional dataset.
    df_columns : list of str
        A list of column names for the datasets.
    ds_names : list of str
        A list of dataset names used for computation.
    ds_names_lut : dict
        A dictionary mapping column names to dataset names.
    tds_names : list of str
        A list of two-dataset combination names.
    result_template : dict
        A template dictionary for storing computed metrics.
    seasons : list of str
        A list of seasons for which metrics are computed.
    month_to_season : numpy.ndarray
        A mapping from month index to season name.

    Methods
    -------
    calc_metrics(data, gpi_info)
        Calculate the desired statistics.

        Parameters
        ----------
        data : pandas.DataFrame
            A pandas DataFrame containing the datasets to compare.
        gpi_info : tuple
            Grid point information (gpi, lon, lat).

        Returns
        -------
        dict
            A dictionary containing computed metrics for the provided data and grid point.
    """

    def __init__(self, other_name1='k1', other_name2='k2',
                 dataset_names=None, seasonal_metrics=False, swi_option=True):

        self.swi_option = swi_option

        self.other_name1 = other_name1
        self.other_name2 = other_name2
        self.df_columns = ['ref', self.other_name1, self.other_name2]

        self.ds_names = dataset_names if dataset_names is not None else self.df_columns

        # create lut between df columns and dataset names
        self.ds_names_lut = {}
        for name, col in zip(self.ds_names, self.df_columns):
            self.ds_names_lut[col] = name

        self.tds_names = []
        for combi in itertools.combinations(self.df_columns, 2):
            self.tds_names.append("{:}_and_{:}".format(*combi))

        self.tds_names_extend = []
        for combi in itertools.combinations(self.df_columns, 2):
            self.tds_names_extend.append("{:}_and_{:}".format(*combi))
            self.tds_names_extend.append("{:}_and_{:}".format(combi[1], combi[0]))

        # metrics' dictionaries
        metrics_common = {'n_obs': np.int32([0])}

        metrics_sds = {
            'nonan_n_obs': np.float32([np.nan]),
            'perc_nan': np.float32([np.nan]),
            'var': np.float32([np.nan]),
            'mean': np.float32([np.nan]),
            'max': np.float32([np.nan]),
            'min': np.float32([np.nan]),
            'snr': np.float32([np.nan]),
            'err_var': np.float32([np.nan]),
            'beta': np.float32([np.nan]),
            'w': np.float32([np.nan]),
        }

        metrics_tds = {
            'R': np.float32([np.nan]),
            'p_R': np.float32([np.nan]),
            'rho': np.float32([np.nan]),
            'p_rho': np.float32([np.nan]),
            'bias': np.float32([np.nan]),
            'ubrmsd': np.float32([np.nan])
        }

        metrics_tds_extend = {
            'w': np.float32([np.nan])
        }

        metrics_only = {
            'swi_t' : np.float32([np.nan]),
            'swi_nse' : np.float32([np.nan])
        }

        self.result_template = {'gpi': np.int32([-1]),
                                'lon': np.float32([np.nan]),
                                'lat': np.float32([np.nan]),
                                "status": np.int32([eh.UNCAUGHT])
        }

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

            # get template for two-dataset metric, extended names for weights
            for tds_name_e in self.tds_names_extend:
                split_tds_name = tds_name_e.split('_and_')
                tds_name_key = "{:}_{:}".format(self.ds_names_lut[
                    split_tds_name[0]],
                    self.ds_names_lut[
                    split_tds_name[1]])
                for metric in metrics_tds_extend.keys():
                    key = "{:}_{:}_{:}".format(tds_name_key, season, metric)
                    self.result_template[key] = metrics_tds_extend[metric].copy()

            # get template for SWI T value, only for SMAP metrics
            for metric in metrics_only.keys():
                key = "{:}_{:}_{:}".format('SMAP', season, metric)
                self.result_template[key] = metrics_only[metric].copy()

        self.month_to_season = np.array(['', 'DJF', 'DJF', 'MAM', 'MAM',
                                         'MAM', 'JJA', 'JJA', 'JJA', 'SON',
                                         'SON', 'SON', 'DJF'])

    def calc_metrics(self, data, gpi_info):
        """
        Calculate the desired statistics.

        Parameters
        ----------
        data : pandas.DataFrame
            A DataFrame containing the datasets to compare.
        gpi_info : tuple
            Grid point information (gpi, lon, lat).

        Returns
        -------
        dict
            A dictionary containing computed metrics for the provided data and grid point.

        Notes
        -----
        Calculated metrics rely on pytesmo.df_metrics, which contains wrappers for the metrics
        with the same names from scipy.stats. Wrappers are needed in order to compute all
        metrics pairwise if triplets of datasets are given.
        """

        # info start metrics
        logging.info(' ------> Compute metrics ... ')

        # get datasets
        dataset = copy.deepcopy(self.result_template)
        dataset['gpi'][0] = gpi_info[0]
        dataset['lon'][0] = gpi_info[1]
        dataset['lat'][0] = gpi_info[2]

        # define season-based subset
        for season in self.seasons:

            # info start season
            logging.info(' ------> Season "' + season + '" ... ')

            if season != 'ALL':
                subset = self.month_to_season[data.index.month] == season
            else:
                subset = np.ones(len(data), dtype=bool)

            # get single dataset metrics
            # calculate SNR
            """
            SNR is calculated via metrics.tcol_snr.
            DEPRECATED: Use the function `tcol_metrics` instead.
            NOTE: `tcol_metrics` is found in pytesmo.df_metrics and is a wrapper
            around the base method ``tcol_metrics`` from ``pytesmo.metrics.tcol``.

            Triple collocation based estimation of signal-to-noise ratio, absolute
            errors, and rescaling coefficients.

            Parameters
            ----------
            x: 1D numpy.ndarray
                first input dataset
            y: 1D numpy.ndarray
                second input dataset
            z: 1D numpy.ndarray
                third input dataset
            ref_ind: int
                Index of reference data set for estimating scaling
                coefficients. Default: 0 (x)

            Returns
            -------
            snr : namedtuple
                signal-to-noise (variance) ratio [dB] from the named columns.
            err_std_dev : namedtuple
                **SCALED** error standard deviation from the named columns
            beta : namedtuple
                Scaling coefficients (i_scaled = i * beta_i)

            Notes
            -----
            
            A complete documentation of this function can be found in metrics_docs.rst.

            References
            ----------
            .. [Gruber2015] Gruber, A., C. -H. Su, S. Zwieback, W. Crow, W. Dorigo,
            and W. Wagner. “Recent Advances in (Soil Moisture) Triple Collocation
            Analysis.” International Journal of Applied Earth Observation and
            Geoinformation, Advances in the Validation and Application of Remotely
            Sensed Soil Moisture - Part 1, 45 (March 1, 2016): 200–211.
            https://doi.org/10.1016/j.jag.2015.09.002.
            """

            # filter by season and order columns properly
            data = data[['ref', 'k1', 'k2']]
            data = data[subset]

            # number of observations
            n_obs = len(data)
            
            data_vectors = [data[col].values for col in data.columns]
            
            # check data usability
            skip = False
            if ExtendedMetrics.check_data(data_vectors, self.ds_names) == 0:
                skip = True
            if skip:
                # info end season
                dataset["status"][0] = eh.NO_GPI_DATA
                logging.info(' ------> Season "' + season + '" ... FAILED')
                continue

            # other version with nwise_apply (more elegant same results)
            # apply_cov = lambda x,y,z : np.cov(np.vstack((x,y,z)))
            # cov = df_metrics.nwise_apply(data, apply_cov, n=3, comm=True, must_include=[0])
            if self.swi_option:
                ## fit SWI(SMAP) against ref (HMC)
                # define dataframe with only ref and SMAP
                df_hmc_smap = data.drop(columns='k1', inplace=False)
                df_hmc_smap.dropna(inplace=True)

                # get julian dates of time series
                jd = df_hmc_smap.index.to_julian_date().values

                # prepare input data by rescaling them on their min/max
                # ref. Albergel et al., 2008, "From near-surface to root zone...

                # min-max normalization
                norm_min_max = lambda x: (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))

                sm_ts = {}
                sm_ts['sm_jd'] = jd
                sm_ts['sm'] = df_hmc_smap.k2.values.astype(np.float32)

                sm_ts_norm = {}
                sm_ts_norm['sm_jd'] = jd
                sm_ts_norm['sm'] = norm_min_max(df_hmc_smap.k2.values.astype(np.float32))

                sm_ts_ref = {}
                sm_ts_ref['sm_jd'] = jd
                sm_ts_ref['sm'] = df_hmc_smap.ref.values.astype(np.float32)

                sm_ts_ref_norm = {}
                sm_ts_ref_norm['sm_jd'] = jd
                sm_ts_ref_norm['sm'] = norm_min_max(df_hmc_smap.ref.values.astype(np.float32))

                # perform fit of t parameter against ref observations
                # curve fit cannot be used on integer values, and t is integer
                # use scipy.optimize.brute to minimize an error function
                t_values = (slice(1, 41, 1),) # create slice of t values from 1 to 40
                t_swi = brute(ExtendedMetrics.errfun, t_values, args=(sm_ts_norm, sm_ts_ref_norm['sm'],), finish=None)

                t_swi = 5

                # calculate ts of swi (on not-norm data) and merge the series with data df
                swi_ts = ExtendedMetrics.calc_swi(sm_ts, [int(t_swi)])
                swi_ts = swi_ts[f'swi_{int(t_swi)}']
                swi_ts_df = pd.Series(data=swi_ts, index=df_hmc_smap.index, name='k2_swi')
                data = data.merge(swi_ts_df, how='left', left_index=True, right_index=True)
                nse = he.evaluator(he.nse, swi_ts, sm_ts_ref['sm'])

                # out variables on nc file
                dataset['{:}_{:}_swi_t'.format('SMAP', season)][0] = t_swi
                dataset['{:}_{:}_swi_nse'.format('SMAP', season)][0] = nse

            # instead of dropping observations from all products when any is nan,
            # following metrics are computed product-wise
            # snr, err, beta instead drop all observations when one is nan
            x = data.ref.values
            y = data.k1.values
            if self.swi_option:
                z = data.k2_swi.values
            else:
                z = data.k2.values

            data_vectors = [x, y, z]

            nnan_obs  = len(data) - data.isna().sum().values
            perc_nans = data.isna().sum().values / len(data)
            nanvars   = [np.nanvar(e)  for e in data_vectors]
            nanmeans  = [np.nanmean(e) for e in data_vectors]
            nanmaxs   = [np.nanmax(e)  for e in data_vectors]
            nanmins   = [np.nanmin(e)  for e in data_vectors]

            # compute other metrics using whole dataframe
            # since the annidated function df_metrics.nwise_apply drops nan values
            snr, err, beta = df_metrics.tcol_metrics(data)
            w = ExtendedMetrics.calc_weights(snr[0])

            # calculate weights pair-wise
            w_pairs = ExtendedMetrics.calc_weights_pairs(snr[0])
            w_pairs_dict = w_pairs._asdict()

            # out variables on nc file
            dataset['{:}_n_obs'.format(season)][0] = n_obs
            for i, name in enumerate(self.ds_names):
                dataset['{:}_{:}_nonan_n_obs'.format(name, season)][0] = nnan_obs[i]
                dataset['{:}_{:}_perc_nan'.format(name, season)][0] = perc_nans[i]
                dataset['{:}_{:}_var'.format(name, season)][0] = nanvars[i]
                dataset['{:}_{:}_mean'.format(name, season)][0] = nanmeans[i]
                dataset['{:}_{:}_max'.format(name, season)][0] = nanmaxs[i]
                dataset['{:}_{:}_min'.format(name, season)][0] = nanmins[i]
                dataset['{:}_{:}_snr'.format(name, season)][0] = snr[0][i]
                dataset['{:}_{:}_err_var'.format(name, season)][0] = err[0][i]
                dataset['{:}_{:}_beta'.format(name, season)][0] = beta[0][i]
                dataset['{:}_{:}_w'.format(name, season)][0] = w[i]

            for tds_name_e in self.tds_names_extend:
                w = w_pairs_dict[tds_name_e]
                split_tds_name = tds_name_e.split('_and_')
                tds_name_key = "{:}_{:}".format(self.ds_names_lut[
                                                    split_tds_name[0]],
                                                self.ds_names_lut[
                                                    split_tds_name[1]])
                dataset['{:}_{:}_w'.format(tds_name_key, season)][0] = w

            # calculate Pearson correlation
            r"""
            The function df_metrics.pearsonr a wrapper around
            scipy.stats.pearsonr, which provides the correlation R and
            the p value of the related hypothesis test, assuming no
            correlation between the data.
            
            The correlation coefficient is calculated as follows:
                
            .. math::
                
                r = \frac{\sum (x - m_x) (y - m_y)}
                        {\sqrt{\sum (x - m_x)^2 \sum (y - m_y)^2}}
                
            where :math:`m_x` is the mean of the vector :math:`x` and :math:`m_y` is
            the mean of the vector :math:`y`.
            
            The p-value roughly indicates the probability of an uncorrelated system
            producing datasets that have a Pearson correlation at least as extreme
            as the one computed from these datasets. The p-value calculated here is
            two-tailed.
            
            """
            pearson_R, pearson_p = df_metrics.pearsonr(data)
            pearson_R = pearson_R._asdict()
            pearson_p = pearson_p._asdict()

            # calculate Spearman correlation
            r"""
            To calculate Spearman correlation, one has to use the
            ``df_metrics.spearmanr`` method, which calculates a Spearman
            correlation coefficient with associated two-tailed p-value.
            The Spearman correlation coefficient is a measure of the
            monotonicity of the relationship between two datasets. It
            does not imply that the data are normally distributed, nor
            it implies a linear correlation between them.
            
            """
            spea_rho, spea_p = df_metrics.spearmanr(data)
            spea_rho = spea_rho._asdict()
            spea_p = spea_p._asdict()

            # calculate bias
            r"""
            Bias of two datasets is just the difference between
            the mean values of the two.
            """
            bias_nT = df_metrics.bias(data)
            bias_dict = bias_nT._asdict()

            # calculate ubRMSD
            r"""
            Unbiased root-mean-square deviation (uRMSD).
            This corresponds to RMSD with mean biases removed beforehand, that is

            ..math::

                ubRMSD = \sqrt{\frac{1}{n}\sum\limits_{i=1}^n
                                   \left((x - \bar{x}) - (y - \bar{y}))^2}

            NOTE: If you are scaling the data beforehand to have zero mean bias, this
            is exactly the same as RMSD.
            
            """
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

            # info end season
            dataset["status"][0] = eh.OK
            logging.info(' ------> Season "' + season + '" ... DONE')

        # info end metrics
        logging.info(' ------> Compute metrics ... DONE')

        return dataset

    @staticmethod
    def calc_weights(snr):
        r"""Calculate the weights of the datasets for merging them.

        Given a tuple of SNR values, calculated following covariance
        notation [Gruber2015]_, returns the respective weights as
        in [Gruber2017]_ (Eq. 4),

        .. math::

            w_i = \frac{SNR_i}{\sum{SNR_j}}

        where i, j = x,y,z.

        Parameters
        ----------
        snr: numpy.ndarray
            Array of ordered SNR (Signal to Noise Ratio) of the datasets.

        Returns
        -------
        w: numpy.ndarray
            Array of ordered weights of the datasets for merging.

        References
        ----------
        .. [Gruber2015] Gruber, A., C. -H. Su, S. Zwieback, W. Crow, W. Dorigo,
            and W. Wagner. “Recent Advances in (Soil Moisture) Triple Collocation
            Analysis.” International Journal of Applied Earth Observation and
            Geoinformation, Advances in the Validation and Application of Remotely
            Sensed Soil Moisture - Part 1, 45 (March 1, 2016): 200–211.
            https://doi.org/10.1016/j.jag.2015.09.002.
        .. [Gruber2017] A. Gruber, W. A. Dorigo, W. Crow and W. Wagner, "Triple
            Collocation-Based Merging of Satellite Soil Moisture Retrievals," in
            IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 12,
            pp. 6780-6792, Dec. 2017, doi: 10.1109/TGRS.2017.2734070

        """
        # snr values are set to linear scale
        db2lin = lambda x : 10**(x/10)
        snr_lin = [db2lin(s) for s in snr]

        # calculate w[i]
        w_dict = {name: snr_lin[i]/np.sum(snr_lin) for i, name in enumerate(snr._fields)}
        w = collections.namedtuple('w', w_dict.keys())(**w_dict)

        # check if each w[i] is in 0-1 range
        for we, name in zip(w, w._fields):
            if we < 0. or we > 1.1:
                logging.warning('===> Weight of {:} is negative or exceeding 1. Check your stats.'.format(name))

        # check if sum of w[i] is 1
        if np.sum(w) > 1.1:
            logging.warning('===> Sum of weights exceeds 1. Check your stats.')

        return w

    @staticmethod
    def calc_weights_pairs(snr: collections.namedtuple):

        w2 = lambda x, y: x / (x + y)
        # snr values are set to linear scale
        db2lin = lambda x : 10**(x/10)
        snr_dict = {key:db2lin(s) for key, s in zip(snr._fields, snr)}

        w_pairs_dict = {}
        for combi in itertools.combinations(snr_dict.keys(), 2):
            e0, e1 = combi
            key1 = "{:}_and_{:}".format(e0, e1)
            w_pairs_dict[key1] = w2(snr_dict[e0], snr_dict[e1])
            key2 = "{:}_and_{:}".format(e1, e0)
            w_pairs_dict[key2] = w2(snr_dict[e1], snr_dict[e0])

        return collections.namedtuple('w_pairs', w_pairs_dict.keys())(**w_pairs_dict)

    @staticmethod
    def check_data(data_vectors, names):

        # info start check metrics
        logging.info(' -------> Datasets check ... ')
        # Initialize variables to track checks
        all_numeric = True
        consistent_shape = True
        minimum_lenght = True

        # Check data types, shapes, and for non-numeric values
        first_shape = None
        for e, name in zip(data_vectors, names):

            # check if data are not empty
            if e.shape[0] > 0:
                if not isinstance(e, (np.ndarray, list)):
                    logging.warning(f' --------> An element is not a NumPy array or list: "{type(e)}". '
                                    f'INCONSISTENT DATA TYPE')
                    all_numeric = False
                elif isinstance(e, np.ndarray):
                    if not np.issubdtype(e.dtype, np.number):
                        logging.warning(f' --------> An element contains a non-numeric data type: "{e.dtype}". '
                                        f'NON-NUMERIC VALUES')
                        all_numeric = False

                if first_shape is None:
                    first_shape = e.shape
                elif e.shape != first_shape:
                    logging.warning(f' --------> Elements have inconsistent shapes: "{first_shape}" and "{e.shape}". '
                                    f'INCONSISTENT SHAPES')
                    consistent_shape = False

                if np.sum(~np.isnan(e)) <= 10:
                    logging.warning(f' --------> No-nan values in "{name}" slice are <= 10. '
                                    f'NOT ENOUGH VALID DATA')
                    minimum_lenght = False
            else:
                logging.warning(f' --------> No data in "{name}" slice. NO DATA')
                all_numeric, consistent_shape, minimum_lenght = False, False, False

        # Check if all checks passed
        if all_numeric and consistent_shape and minimum_lenght:
            # info end check datasets
            logging.info(' --------> All datasets check passed')
            logging.info(' -------> Datasets check ... DONE')
            return 1
        else:
            # info end check datasets
            logging.warning(' --------> Datasets for this point are not eligible to compute metrics. Skip this point')
            logging.info(' -------> Datasets check ... FAILED')
            return 0

    @staticmethod
    def calc_swi(sm_ts, t):
        """wrapper around pyswi.calc_swi_ts"""
        jd = np.squeeze(sm_ts['sm_jd'])
        dtype = np.dtype([('sm_jd', np.float64), ('sm', np.float32)])
        stack = np.hstack((jd[:, np.newaxis], sm_ts['sm'][:, np.newaxis]))
        data = unstructured_to_structured(stack, dtype=dtype)
        swi_ts, gain_out = calc_swi_ts(data, jd, t_value=t)
        return swi_ts

    @staticmethod
    def errfun(t, *args):
        sm_ts, ref_ts = args
        swi_ts = ExtendedMetrics.calc_swi(sm_ts, t)
        swi_ts = swi_ts[f'swi_{np.squeeze(t)}'].astype(np.float32)

        # maximize correlation coefficient
        # errsqr = -np.corrcoef(x=swi_ts, y=ref_ts)[0,1]

        # maximize Nash-Sutcliff efficiency
        errsqr = -he.evaluator(he.nse, swi_ts, ref_ts)
        return errsqr


# -------------------------------------------------------------------------------------
