"""
Library Features:

Name:          lib_utils_plot
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231212'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

from lib_info_args import logger_name

# set logger obj
alg_logger = logging.getLogger(logger_name)

# debugging
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot data 2d
def plot_data_2d(var_data, var_geo_x=None, var_geo_y=None, var_min=0, var_max=1):

    plt.figure(figsize=(10, 10))
    plt.imshow(var_data)
    plt.colorbar()
    plt.clim(var_min, var_max)

    plt.show()
# ----------------------------------------------------------------------------------------------------------------------
