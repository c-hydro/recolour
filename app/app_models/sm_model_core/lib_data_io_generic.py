"""
Library Features:

Name:          lib_data_io_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260615'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import pandas as pd

from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to combine data over the expected time range
def combine_data_point_by_time(
        dframe_k1,
        dframe_k2,
        dframe_k3,
        time_tag='time',
        time_frequency='h',
        time_reverse=True,
        time_ref=None,
        fill_value_missing=-9997,
        no_data_k1=-9999,
        no_data_k2=-9999,
        no_data_k3=-9999,
        scale_factor_k1=1,
        scale_factor_k2=1,
        scale_factor_k3=0.01,
        status_tag='status'):

    """
    Combine three time-series dataframes over a common expected time range.

    The final data values remain NaN when unavailable. A bit-mask column
    describes which datasets are unavailable and whether a timestamp was
    introduced while building the expected time range.

    Status flags:
        0  : all datasets are available
        1  : dataset k1 is unavailable
        2  : dataset k2 is unavailable
        4  : dataset k3 is unavailable
        8  : timestamp is absent from all original datasets

    Flags can be combined:
        3  : k1 and k2 unavailable
        5  : k1 and k3 unavailable
        6  : k2 and k3 unavailable
        7  : all datasets unavailable, but timestamp existed in at least
             one original dataframe
        15 : timestamp absent from every original dataframe and all values
             are unavailable

    Parameters
    ----------
    fill_value_missing : int or float
        Retained for backward compatibility. Final missing values are stored
        as NaN; the original configured value is saved in dataframe attrs.
    """

    # -------------------------------------------------------------------------
    # Status bit-mask definitions
    status_ok = 0
    status_k1_missing = 1
    status_k2_missing = 2
    status_k3_missing = 4
    status_time_missing = 8

    value_tag_k1 = 'values_k1'
    value_tag_k2 = 'values_k2'
    value_tag_k3 = 'values_k3'

    value_tags = [
        value_tag_k1,
        value_tag_k2,
        value_tag_k3
    ]

    time_frequency = str(time_frequency).lower()

    # -------------------------------------------------------------------------
    # Check input dataframes
    if (
            dframe_k1 is None or
            dframe_k2 is None or
            dframe_k3 is None):

        log_stream.warning(
            ' ===> Dataframes are not defined; '
            'one or more dataframes are None'
        )
        return None

    if (
            dframe_k1.empty or
            dframe_k2.empty or
            dframe_k3.empty):

        log_stream.warning(
            ' ===> One or more input dataframes are empty'
        )
        return None

    # Work on copies to avoid changing the original dataframes
    dframe_k1 = dframe_k1.copy()
    dframe_k2 = dframe_k2.copy()
    dframe_k3 = dframe_k3.copy()

    # Copy attributes from the first dataframe
    attrs_common = dict(
        getattr(dframe_k1, 'attrs', {}) or {}
    )

    # -------------------------------------------------------------------------
    # Ensure datetime indexes
    try:
        dframe_k1.index = pd.DatetimeIndex(dframe_k1.index)
        dframe_k2.index = pd.DatetimeIndex(dframe_k2.index)
        dframe_k3.index = pd.DatetimeIndex(dframe_k3.index)
    except Exception as exc:
        log_stream.error(
            ' ===> Unable to convert dataframe indexes to datetime'
        )
        raise RuntimeError(
            'All dataframe indexes must contain valid datetime values'
        ) from exc

    # Remove invalid timestamps
    dframe_k1 = dframe_k1.loc[~dframe_k1.index.isna()]
    dframe_k2 = dframe_k2.loc[~dframe_k2.index.isna()]
    dframe_k3 = dframe_k3.loc[~dframe_k3.index.isna()]

    if (
            dframe_k1.empty or
            dframe_k2.empty or
            dframe_k3.empty):

        log_stream.warning(
            ' ===> One or more dataframes have no valid timestamps'
        )
        return None

    # Sort indexes
    dframe_k1 = dframe_k1.sort_index()
    dframe_k2 = dframe_k2.sort_index()
    dframe_k3 = dframe_k3.sort_index()

    # -------------------------------------------------------------------------
    # Check duplicate timestamps
    if dframe_k1.index.has_duplicates:
        log_stream.warning(
            ' ===> Duplicate timestamps found in dataframe k1; '
            'keeping the last occurrence'
        )
        dframe_k1 = dframe_k1[
            ~dframe_k1.index.duplicated(keep='last')
        ]

    if dframe_k2.index.has_duplicates:
        log_stream.warning(
            ' ===> Duplicate timestamps found in dataframe k2; '
            'keeping the last occurrence'
        )
        dframe_k2 = dframe_k2[
            ~dframe_k2.index.duplicated(keep='last')
        ]

    if dframe_k3.index.has_duplicates:
        log_stream.warning(
            ' ===> Duplicate timestamps found in dataframe k3; '
            'keeping the last occurrence'
        )
        dframe_k3 = dframe_k3[
            ~dframe_k3.index.duplicated(keep='last')
        ]

    # Save the original timestamp indexes. These are used to distinguish
    # unavailable values from timestamps absent from all source datasets.
    time_index_k1 = dframe_k1.index.copy()
    time_index_k2 = dframe_k2.index.copy()
    time_index_k3 = dframe_k3.index.copy()

    # -------------------------------------------------------------------------
    # Define common time range
    time_start_k1 = dframe_k1.index.min()
    time_end_k1 = dframe_k1.index.max()

    time_start_k2 = dframe_k2.index.min()
    time_end_k2 = dframe_k2.index.max()

    time_start_k3 = dframe_k3.index.min()
    time_end_k3 = dframe_k3.index.max()

    time_start_common = min(
        time_start_k1,
        time_start_k2,
        time_start_k3
    )

    time_end_common = max(
        time_end_k1,
        time_end_k2,
        time_end_k3
    )

    if time_ref is not None:
        time_ref = pd.Timestamp(time_ref)
        time_end_common = max(time_end_common, time_ref)

    time_range_common = pd.date_range(
        start=time_start_common,
        end=time_end_common,
        freq=time_frequency
    )

    # -------------------------------------------------------------------------
    # Organize dataset k1
    dframe_k1 = dframe_k1.drop(
        columns=[time_tag],
        errors='ignore'
    )

    if value_tag_k1 not in dframe_k1.columns:
        log_stream.error(
            ' ===> Dataframe 1 does not have the column "values_k1"'
        )
        raise RuntimeError(
            'Column "values_k1" must be included in dataframe 1'
        )

    dframe_k1 = dframe_k1[[value_tag_k1]]

    dframe_k1[value_tag_k1] = pd.to_numeric(
        dframe_k1[value_tag_k1],
        errors='coerce'
    )

    dframe_k1.loc[
        dframe_k1[value_tag_k1] == no_data_k1,
        value_tag_k1
    ] = np.nan

    dframe_k1[value_tag_k1] = (
        dframe_k1[value_tag_k1] * scale_factor_k1
    )

    dframe_k1.loc[
        ~np.isfinite(dframe_k1[value_tag_k1]),
        value_tag_k1
    ] = np.nan

    # -------------------------------------------------------------------------
    # Organize dataset k2
    dframe_k2 = dframe_k2.drop(
        columns=[time_tag],
        errors='ignore'
    )

    if value_tag_k2 not in dframe_k2.columns:
        log_stream.error(
            ' ===> Dataframe 2 does not have the column "values_k2"'
        )
        raise RuntimeError(
            'Column "values_k2" must be included in dataframe 2'
        )

    dframe_k2 = dframe_k2[[value_tag_k2]]

    dframe_k2[value_tag_k2] = pd.to_numeric(
        dframe_k2[value_tag_k2],
        errors='coerce'
    )

    dframe_k2.loc[
        dframe_k2[value_tag_k2] == no_data_k2,
        value_tag_k2
    ] = np.nan

    dframe_k2[value_tag_k2] = (
        dframe_k2[value_tag_k2] * scale_factor_k2
    )

    dframe_k2.loc[
        ~np.isfinite(dframe_k2[value_tag_k2]),
        value_tag_k2
    ] = np.nan

    # -------------------------------------------------------------------------
    # Organize dataset k3
    dframe_k3 = dframe_k3.drop(
        columns=[time_tag],
        errors='ignore'
    )

    if value_tag_k3 not in dframe_k3.columns:
        log_stream.error(
            ' ===> Dataframe 3 does not have the column "values_k3"'
        )
        raise RuntimeError(
            'Column "values_k3" must be included in dataframe 3'
        )

    dframe_k3 = dframe_k3[[value_tag_k3]]

    dframe_k3[value_tag_k3] = pd.to_numeric(
        dframe_k3[value_tag_k3],
        errors='coerce'
    )

    dframe_k3.loc[
        dframe_k3[value_tag_k3] == no_data_k3,
        value_tag_k3
    ] = np.nan

    dframe_k3[value_tag_k3] = (
        dframe_k3[value_tag_k3] * scale_factor_k3
    )

    dframe_k3.loc[
        ~np.isfinite(dframe_k3[value_tag_k3]),
        value_tag_k3
    ] = np.nan

    # -------------------------------------------------------------------------
    # Join datasets over the expected time range
    dframe_common = pd.DataFrame(
        index=time_range_common
    )
    dframe_common.index.name = time_tag

    dframe_common = dframe_common.join(
        dframe_k1,
        how='left'
    )
    dframe_common = dframe_common.join(
        dframe_k2,
        how='left'
    )
    dframe_common = dframe_common.join(
        dframe_k3,
        how='left'
    )

    # -------------------------------------------------------------------------
    # Limit output to the reference time
    if time_ref is not None:
        dframe_common = dframe_common.loc[
            dframe_common.index <= time_ref
        ]

    # -------------------------------------------------------------------------
    # Create status bit mask
    status_values = np.zeros(
        dframe_common.shape[0],
        dtype=np.uint8
    )

    mask_k1_missing = (
        dframe_common[value_tag_k1]
        .isna()
        .to_numpy()
    )
    mask_k2_missing = (
        dframe_common[value_tag_k2]
        .isna()
        .to_numpy()
    )
    mask_k3_missing = (
        dframe_common[value_tag_k3]
        .isna()
        .to_numpy()
    )

    status_values[mask_k1_missing] |= status_k1_missing
    status_values[mask_k2_missing] |= status_k2_missing
    status_values[mask_k3_missing] |= status_k3_missing

    # A time is marked as missing only when it did not exist in any of
    # the three original dataframe indexes.
    mask_time_available = (
        dframe_common.index.isin(time_index_k1) |
        dframe_common.index.isin(time_index_k2) |
        dframe_common.index.isin(time_index_k3)
    )

    mask_time_missing = ~mask_time_available

    status_values[mask_time_missing] |= status_time_missing

    dframe_common[status_tag] = status_values

    # -------------------------------------------------------------------------
    # Add time as a regular column
    if time_tag not in dframe_common.columns:
        dframe_common[time_tag] = dframe_common.index

    # Define final column order
    dframe_common = dframe_common[
        [
            time_tag,
            value_tag_k1,
            value_tag_k2,
            value_tag_k3,
            status_tag
        ]
    ]

    # -------------------------------------------------------------------------
    # Sort output
    if time_reverse:
        dframe_common = dframe_common.sort_index(
            ascending=False
        )
    else:
        dframe_common = dframe_common.sort_index(
            ascending=True
        )

    # -------------------------------------------------------------------------
    # Add descriptive attributes
    attrs_common['time_column'] = time_tag
    attrs_common['time_frequency'] = time_frequency

    attrs_common['status_column'] = status_tag
    attrs_common['status_dtype'] = 'uint8'
    attrs_common['status_encoding'] = 'bit_mask'
    attrs_common['status_description'] = (
        'Bit mask describing dataset availability. '
        'Multiple flags are combined using bitwise OR.'
    )

    attrs_common['status_flags'] = {
        'ok': status_ok,
        'k1_missing': status_k1_missing,
        'k2_missing': status_k2_missing,
        'k3_missing': status_k3_missing,
        'time_missing': status_time_missing
    }

    attrs_common['status_flag_meanings'] = {
        status_ok: 'All datasets are available',
        status_k1_missing: 'Dataset k1 is unavailable',
        status_k2_missing: 'Dataset k2 is unavailable',
        status_k3_missing: 'Dataset k3 is unavailable',
        status_time_missing: (
            'Timestamp is absent from all original datasets'
        )
    }

    attrs_common['status_codes'] = {
        0: 'All datasets available',
        1: 'Dataset k1 unavailable',
        2: 'Dataset k2 unavailable',
        3: 'Datasets k1 and k2 unavailable',
        4: 'Dataset k3 unavailable',
        5: 'Datasets k1 and k3 unavailable',
        6: 'Datasets k2 and k3 unavailable',
        7: (
            'All datasets unavailable, but timestamp exists '
            'in at least one original dataset'
        ),
        8: 'Timestamp absent from all original datasets',
        9: 'Timestamp absent and dataset k1 unavailable',
        10: 'Timestamp absent and dataset k2 unavailable',
        11: 'Timestamp absent and datasets k1 and k2 unavailable',
        12: 'Timestamp absent and dataset k3 unavailable',
        13: 'Timestamp absent and datasets k1 and k3 unavailable',
        14: 'Timestamp absent and datasets k2 and k3 unavailable',
        15: (
            'Timestamp absent from all original datasets and '
            'all dataset values unavailable'
        )
    }

    attrs_common['no_data_values_original'] = {
        value_tag_k1: no_data_k1,
        value_tag_k2: no_data_k2,
        value_tag_k3: no_data_k3
    }

    attrs_common['scale_factors'] = {
        value_tag_k1: scale_factor_k1,
        value_tag_k2: scale_factor_k2,
        value_tag_k3: scale_factor_k3
    }

    attrs_common['missing_value_output'] = 'NaN'
    attrs_common['legacy_fill_value_missing'] = fill_value_missing

    dframe_common.attrs = attrs_common

    return dframe_common, attrs_common
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to combine data over the expected time range
def combine_data_point_by_time_OLD(dframe_k1, dframe_k2, dframe_k3,
                               time_tag='time', time_frequency='h', time_reverse=True,
                               time_ref=None, fill_value_missing=-9997,
                               no_data_k1=-9999, no_data_k2=-9999, no_data_k3=-9999,
                               scale_factor_k1=1, scale_factor_k2=1, scale_factor_k3=0.01):

    time_frequency = time_frequency.lower()

    if (dframe_k1 is None) or (dframe_k2 is None) or (dframe_k3 is None):
        log_stream.warning(' ===> Dataframes are not defined; One or more dataframes are defined by None')
        return None

    time_start_k1, time_end_k1 = dframe_k1.index.min(), dframe_k1.index.max()
    time_start_k2, time_end_k2 = dframe_k2.index.min(), dframe_k2.index.max()
    time_start_k3, time_end_k3 = dframe_k3.index.min(), dframe_k3.index.max()

    time_start_common = pd.DatetimeIndex([time_start_k1, time_start_k2, time_start_k3]).min()
    time_end_common = pd.DatetimeIndex([time_end_k1, time_end_k2, time_end_k3]).max()

    if time_ref is not None:
        time_ref = pd.Timestamp(time_ref)
        time_end_common = max(time_end_common, time_ref)

    time_range_common = pd.date_range(time_start_common, time_end_common, freq=time_frequency)

    attrs_common = dframe_k1.attrs

    dframe_k1 = dframe_k1.drop(columns=[time_tag], errors='ignore')
    if 'values_k1' not in dframe_k1.columns:
        log_stream.error(' ===> Dataframe 1 does not have the column "values_k1"')
        raise RuntimeError('Column "values_k1" must be included in the dataframe.')

    dframe_k1.loc[dframe_k1['values_k1'] == no_data_k1, 'values_k1'] = np.nan
    dframe_k1['values_k1'] = dframe_k1['values_k1'].values * scale_factor_k1
    dframe_k1.loc[np.isnan(dframe_k1['values_k1']), 'values_k1'] = no_data_k1

    dframe_k2 = dframe_k2.drop(columns=[time_tag], errors='ignore')
    if 'values_k2' not in dframe_k2.columns:
        log_stream.error(' ===> Dataframe 2 does not have the column "values_k2"')
        raise RuntimeError('Column "values_k2" must be included in the dataframe.')

    dframe_k2.loc[dframe_k2['values_k2'] == no_data_k2, 'values_k2'] = np.nan
    dframe_k2['values_k2'] = dframe_k2['values_k2'].values * scale_factor_k2
    dframe_k2.loc[np.isnan(dframe_k2['values_k2']), 'values_k2'] = no_data_k2

    dframe_k3 = dframe_k3.drop(columns=[time_tag], errors='ignore')
    if 'values_k3' not in dframe_k3.columns:
        log_stream.error(' ===> Dataframe 3 does not have the column "values_k3"')
        raise RuntimeError('Column "values_k3" must be included in the dataframe.')

    dframe_k3.loc[dframe_k3['values_k3'] == no_data_k3, 'values_k3'] = np.nan
    dframe_k3['values_k3'] = dframe_k3['values_k3'].values * scale_factor_k3
    dframe_k3.loc[np.isnan(dframe_k3['values_k3']), 'values_k3'] = no_data_k3

    dframe_common = pd.DataFrame(index=time_range_common)
    dframe_common.index.name = time_tag

    dframe_common = dframe_common.join(dframe_k1)
    dframe_common = dframe_common.join(dframe_k2)
    dframe_common = dframe_common.join(dframe_k3)

    if time_ref is not None:
        dframe_common = dframe_common.loc[dframe_common.index <= time_ref]

        time_range_expected = pd.date_range(
            start=dframe_common.index.min(),
            end=time_ref,
            freq=time_frequency
        )

        dframe_common = dframe_common.reindex(time_range_expected)
        dframe_common.index.name = time_tag

        value_cols = ['values_k1', 'values_k2', 'values_k3']

        missing_rows = dframe_common[value_cols].isna().all(axis=1)
        dframe_common.loc[missing_rows, value_cols] = fill_value_missing

    if time_tag not in dframe_common.columns:
        dframe_common[time_tag] = dframe_common.index

    if time_reverse:
        dframe_common = dframe_common.sort_index(ascending=False)

    dframe_common.attrs = attrs_common

    return dframe_common
# ----------------------------------------------------------------------------------------------------------------------
