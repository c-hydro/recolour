"""
Library Features:

Name:          lib_fx_method
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230727'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import numpy as np

import pytesmo.scaling as default_scaling_methods
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to scale data using a mean-std scaling method (between three datasets)
def mean_std(ref_nrt, ref_dr, other_dr):

    return ((ref_nrt - np.mean(ref_dr)) /
            np.std(ref_dr)) * np.std(other_dr) + np.mean(other_dr)
# ----------------------------------------------------------------------------------------------------------------------
