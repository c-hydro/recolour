"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file obj
def filter_dframe_by_vars(dframe_obj_in: pd.DataFrame,
                          dframe_parameter: str = None,
                          dframe_variables_data: list = None, dframe_variables_grid: list = None):
    if dframe_variables_data is None:
        dframe_variables_data = ['gpi', 'lon', 'lat', 'obs']
    if dframe_variables_grid is None:
        dframe_variables_grid = ['committed_area', 'land_flag']

    dframe_variables_in = dframe_variables_data + dframe_variables_grid
    dframe_variables_out = []
    for variable_name in dframe_variables_in:
        if variable_name in list(dframe_obj_in.columns):
            dframe_variables_out.append(variable_name)
        else:
            logging.warning(' ===> Variable "' + variable_name +
                            ' is not available in the DataFrame. Errors could be occurred')

    dframe_obj_out = dframe_obj_in[dframe_variables_out]

    if dframe_parameter not in dframe_obj_out.columns:
        logging.error(' ===> Variable "' + dframe_parameter + '" must be in the DataFrame')
        raise RuntimeError('Check your DataFrame to correctly set the columns names or data')
    return dframe_obj_out
# ----------------------------------------------------------------------------------------------------------------------


