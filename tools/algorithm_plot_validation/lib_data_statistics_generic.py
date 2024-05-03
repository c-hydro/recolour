"""
Library Features:

Name:          lib_data_statistics_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240411
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import warnings

import numpy as np
import pandas as pd

from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

logging.getLogger('pandas').setLevel(logging.WARNING)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter dataframe nan(s) for stats
def filter_dataframe_nan(
        df_obj: pd.DataFrame,
        var_carea_name: str = 'committed_area', var_carea_value: int = 1,
        var_data_name: str = 'xy_pr'):

    sub_obj = ((df_obj[var_carea_name] == var_carea_value) & (np.isnan(df_obj[var_data_name])))
    value_filtered_nan = sub_obj.sum()

    return value_filtered_nan
# ----------------------------------------------------------------------------------------------------------------------
