
"""
Class Features

Name:          driver_data_io_dynamic_datasets_obs
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210225'
Version:       '1.0.0'

Deps: conda install -y -c conda-forge gdal
"""

######################################################################################
# Library
import logging

import os
import sys
import numpy as np
import scipy.ndimage as ndimage
import subprocess
######################################################################################

#sys.path.append(r'C:\Users\*****\AppData\Local\conda\conda\envs\****\Scripts')

def filter_isolated_cells(array, struct):
    """ Return array with completely isolated single cells removed
    :param array: Array with completely isolated single cells
    :param struct: Structure array for generating unique regions
    :return: Array with minimum region size > 1
    """

    filtered_array = np.copy(array)
    id_regions, num_ids = ndimage.label(filtered_array, structure=struct)
    id_sizes = np.array(ndimage.sum(array, id_regions, range(num_ids + 1)))
    area_mask = (id_sizes == 1)
    filtered_array[area_mask[id_regions]] = 0
    return filtered_array




# -------------------------------------------------------------------------------------
# Method to define sieve command
def compose_sieve_cmd(file_name_in, file_name_out='image_sieved.tiff',
                      sieve_pixel_thr=1, format_out='GTiff',
                      sieve_command='/usr/bin/gdal_sieve.py', sieve_interpreter='python3'):

    sieve_command_list = \
        [sieve_interpreter, sieve_command,
         '-st', str(sieve_pixel_thr), '-8', '-nomask',
         '-of', format_out,
         file_name_in, file_name_out]

    sieve_command_string = ' '.join(sieve_command_list)

    return sieve_command_string
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to execute cmd
def execute_cmd(cmd_string):
    subprocess.call(cmd_string, shell=True)
# -------------------------------------------------------------------------------------ù
