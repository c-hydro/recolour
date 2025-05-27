"""
Library Features:

Name:          lib_utils_method_fx
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210408'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# Libraries
import logging
import pandas as pd
import xarray as xr
import numpy as np

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to convert volume to soil moisture
def convert_volume2sm(values_vtot: np.array, values_vmax: np.array = None, values_geo: np.array = None,
                      values_min: float = 0.0, values_max: float = 100.0) -> np.array:

    values_sm = None
    if values_vmax is not None:
        assert values_vtot.shape == values_vmax.shape, 'volume tot values and volume max values have different dims'

        values_sm = values_vtot / values_vmax
        values_sm = values_sm * 100

        values_sm[np.isnan(values_vmax)] = np.nan
        values_sm[values_sm < values_min] = values_min
        values_sm[values_sm > values_max] = values_max

        if values_geo is not None:
            assert values_vtot.shape == values_geo.shape, 'volume tot values and terrain values have different dims'
            values_sm[values_geo < 0] = np.nan

        """ plot debug
        plt.figure()
        plt.imshow(values_vtot)
        plt.colorbar()

        plt.figure()
        plt.imshow(values_vmax)
        plt.colorbar()

        plt.figure()
        plt.imshow(values_geo)
        plt.colorbar()

        plt.figure()
        plt.imshow(values_sm)
        plt.colorbar()

        plt.show()
        """

    return values_sm
 # -------------------------------------------------------------------------------------
