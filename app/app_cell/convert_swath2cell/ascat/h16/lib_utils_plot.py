"""
Library Features:

Name:          lib_utils_plot
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240423'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import numpy as np

from matplotlib import cm
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to plot points, lons and lats
def plot_points(values, lons, lats, cmap='coolwarm', value_min=None, value_max=None, value_no_data=np.nan):

    if isinstance(values, np.ma.MaskedArray):
        values = values.data

    if value_min is not None:
        values[values < value_min] = value_no_data
    if value_max is not None:
        values[values > value_max] = value_no_data

    fig = plt.figure()
    colors = cm.get_cmap(cmap)(values)[..., :3]
    plt.scatter(lons, lats, c=colors)
    plt.show()

# ----------------------------------------------------------------------------------------------------------------------
