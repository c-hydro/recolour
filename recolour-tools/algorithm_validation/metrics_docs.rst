==============================================
Metrics Calculation in `calc_metrics` Function
==============================================

The `calc_metrics` function calculates various statistics and metrics for a
given dataset and grid point.

Number of Observations (n_obs)
-------------------------------
The function calculates the number of observations (n_obs) for each dataset
based on the provided data. If the number of observations is less than 10, the
calculation for the current season is skipped.

Signal-to-Noise Ratio (SNR)
---------------------------
The signal-to-noise ratio (SNR) is computed using the `tcol_snr` function from
the `metrics` module. It takes three datasets (x, y, and z) and calculates the
SNR, error variance (err), and beta values. These values are computed for each
dataset and each season.

This function estimates the triple collocation errors, the scaling
parameter :math:`\\beta` and the signal to noise ratio directly from the
covariances of the dataset. For a general overview and how this function
and :py:func:`pytesmo.metrics.tcol_error` are related please see
[Gruber2015]_.

Estimation of the error variances from the covariances of the datasets
(e.g. :math:`\\sigma_{XY}` for the covariance between :math:`x` and
:math:`y`) is done using the following formula:

.. math::

   \\sigma_{\\varepsilon_x}^2 =
       \\sigma_{X}^2 - \\frac{\\sigma_{XY}\\sigma_{XZ}}{\\sigma_{YZ}}
.. math::

   \\sigma_{\\varepsilon_y}^2 =
       \\sigma_{Y}^2 - \\frac{\\sigma_{YX}\\sigma_{YZ}}{\\sigma_{XZ}}
.. math::

   \\sigma_{\\varepsilon_z}^2 =
       \\sigma_{Z}^2 - \\frac{\\sigma_{ZY}\\sigma_{ZX}}{\\sigma_{YX}}

:math:`\\beta` can also be estimated from the covariances:

.. math:: \\beta_x = 1
.. math:: \\beta_y = \\frac{\\sigma_{XZ}}{\\sigma_{YZ}}
.. math:: \\beta_z=\\frac{\\sigma_{XY}}{\\sigma_{ZY}}

The signal to noise ratio (SNR) is also calculated from the variances
and covariances:

.. math::

   \\text{SNR}_X[dB] = -10\\log\\left(\\frac{\\sigma_{X}^2\\sigma_{YZ}}
                                     {\\sigma_{XY}\\sigma_{XZ}}-1\\right)
.. math::

   \\text{SNR}_Y[dB] = -10\\log\\left(\\frac{\\sigma_{Y}^2\\sigma_{XZ}}
                                     {\\sigma_{YX}\\sigma_{YZ}}-1\\right)
.. math::

   \\text{SNR}_Z[dB] = -10\\log\\left(\\frac{\\sigma_{Z}^2\\sigma_{XY}}
                                     {\\sigma_{ZX}\\sigma_{ZY}}-1\\right)

It is given in dB to make it symmetric around zero. If the value is zero
it means that the signal variance and the noise variance are equal. +3dB
means that the signal variance is twice as high as the noise variance.


Pearson Correlation (pearson_R, pearson_p)
-------------------------------------------
Pearson correlation coefficients (pearson_R) and their corresponding p-values
(pearson_p) are calculated using the `df_metrics.pearsonr` method. The Pearson
correlation measures the linear relationship between pairs of datasets. The
correlation coefficient and p-value are calculated for each dataset and season.

.. math::

    r = \frac{\sum (x - m_x) (y - m_y)}
            {\sqrt{\sum (x - m_x)^2 \sum (y - m_y)^2}}

The p-value roughly indicates the probability of an uncorrelated system
producing datasets that have a Pearson correlation at least as extreme as the
one computed from these datasets. The p-value calculated here is two-tailed.

Spearman Correlation (spea_rho, spea_p)
---------------------------------------
Spearman correlation coefficients (spea_rho) and their corresponding p-values
(spea_p) are calculated using the `df_metrics.spearmanr` method. Spearman
correlation assesses the monotonic relationship between pairs of datasets.
The coefficients and p-values are computed for each dataset and season.

Bias
----
The bias between pairs of datasets is calculated as the difference between their
mean values.

Unbiased Root-Mean-Square Deviation (ubRMSD)
--------------------------------------------
The unbiased root-mean-square deviation (ubRMSD) is calculated using the
`df_metrics.ubrmsd` method. This is the root-mean-square deviation with mean
biases removed beforehand. It's a measure of the spread between two datasets
after removing the mean bias.

.. math::

    ubRMSD = \sqrt{\frac{1}{n}\sum_{i=1}^n
                   \left((x - \bar{x}) - (y - \bar{y})\right)^2}

Result Storage
--------------
All calculated metrics are stored in the `dataset` dictionary with appropriate
keys that include the dataset names, season, and metric type.

Function Return
---------------
The function returns the `dataset` dictionary containing the computed metrics
for the provided data and grid point.

Notes
----
Some of the metric calculation methods, such as `tcol_snr`,
`df_metrics.pearsonr`, `df_metrics.spearmanr`, `df_metrics.bias`, and
`df_metrics.ubrmsd`, are assumed to be part of a separate module or library
(`metrics`) that contains implementations for these specific calculations.


Metrics' calculation
--------------------

At the moment, the code computes the following metrics, which are implemented
in the `lib_utils_metrics` package from the `ExtendedMetrics` class:

In `__init__`:
    - `metrics_sds`: dict, contains snr, err_var, beta
    - `metrics_tds`: dict, contains R (Pearson correlation coefficient) and
      related two-tailed p-value; rho (Spearman correlation coefficient) and
      related two-tailed p-value; bias; ubRMSD (unbiased root mean standard
      deviation), i.e. RMSD as if bias between the data was removed.

In `calc_metrics`:
    - `snr`, `err`, `beta`: TCA outputs, calculated via
      `pytesmo.metrics.tcol_snr`
    - WARNING: this function is deprecated, and
      `pytesmo.df_metrics.tcol_metrics` is suggested to be used in its place.

NOTE: the function `pytesmo.df_metrics.tcol_metrics` is a wrapper for
triplet-wise calculation around the basic function `pytesmo.tcol.tcol_metrics`.
The latter provides "estimation of signal-to-noise ratio, absolute errors, and
rescaling coefficients" by applying the covariance notation of triple
collocation [Gruber2015]_.

There's another function, which is `pytesmo.tcol.ecol`, which provides
"extended collocation analysis to obtain estimates of
- signal variances
- error variances
- signal-to-noise ratios [dB] (of each single product, _ndr_)
- error cross-covariances (and -correlations)" (pairwise, _ndr_)

Notice that the covariance notation for triple collocation does not require any
matching (e.g. CDF-matching) or bias removal for the datasets involved, unlike
the difference notation.


Proposed changes
----------------
There's no real need to substitute the actual function `metrics.tcol_snr` with
`metrics.tcol_metrics` which belongs to the `tcol` package, since the
"deprecated" statement returns the correct function and is a decorator which is
supposedly used for maintenance purposes only.
Anyways, using `df_metrics.tcol_metrics` is more correct.


Weight estimation
-----------------
Suppose that the assumptions of triple collocation analysis (TCA) are satisfied
[Gruber2015]_. Then TCA can be applied to estimate the optimal weights to use for
merging a triplet of temporally collocated soil moisture products. The weights
are optimal in the least square sense, i.e. that the variance of residual random
errors shall be minimized [Gruber2017]_.
To correctly compute



References
----------
.. [Gruber2015] Gruber, A., Su, C., Zwieback, S., Crow, W., Dorigo, W., Wagner,
  W. (2015). Recent advances in (soil moisture) triple collocation analysis.
  International Journal of Applied Earth Observation and Geoinformation
.. [Gruber2017] A. Gruber, W. A. Dorigo, W. Crow and W. Wagner, "Triple
  Collocation-Based Merging of Satellite Soil Moisture Retrievals," in IEEE Transactions on Geoscience and Remote
Sensing, vol. 55, no. 12, pp. 6780-6792, Dec. 2017, doi: 10.1109/TGRS.2017.2734070.
