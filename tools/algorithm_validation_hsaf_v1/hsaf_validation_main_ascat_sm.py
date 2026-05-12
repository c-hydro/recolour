# ----------------------------------------------------------------------------------------------------------------------
# HSAF validation procedure
# https://pytesmo.readthedocs.io/en/latest/examples.html#the-pytesmo-validation-framework
# TEST over Italy domain --> [1394, 1395]
#
# H16 (time series in indexed format)
# -configfile algorithm_configuration_h16.json
# H113 (data record in contiguous format)
# -configfile algorithm_configuration_h113.json

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import multiprocessing as mp
import argparse

from itertools import repeat

from hsaf_validation_utils import parser_config_file, get_grid_cells, get_group_cells, write_group_cells
from hsaf_validation_process import wrap_process
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def GetArgs():

    # Parser algorithm arg(s)
    parser_obj = argparse.ArgumentParser()
    parser_obj.add_argument('-configfile', action="store", dest="config_file")
    parser_value = parser_obj.parse_args()

    if parser_value.config_file:
        config_file = parser_value.config_file
    else:
        config_file = 'algorithm_configuration.json'

    return config_file

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script Main
if __name__ == "__main__":

    # -------------------------------------------------------------------------------------
    # Get script argument(s)
    config_file = GetArgs()
    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Parser configuration file to get validation parameter(s)
    val_param = parser_config_file(config_file)
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Generate cells list using land flag if grid file exists
    cells, gpis = get_grid_cells(cell_start=int(val_param['cell_start']),
                                 cell_end=int(val_param['cell_end']),
                                 path_grid=val_param['ascat_path_grid'], file_grid='TUW_WARP5_grid_info_2_3.nc')
    # ----------------------------------------------------------------------------------------------------------------------

    # Debug lines
    #print('DEBUG ATTIVO!!!')
    #cells = [1394] # [2128135]
    #val_param['mp_alg'] = False

    # ----------------------------------------------------------------------------------------------------------------------
    # Condition for mp process or n
    if val_param['mp_alg']:

        # ----------------------------------------------------------------------------------------------------------------------
        # Multiprocess option(s)
        if 'mp_cpu' in val_param:
            try:
                mp_cpu = int(val_param['mp_cpu'])
            except:
                mp_cpu = mp.cpu_count()
        else:
            mp_cpu = mp.cpu_count()

        if mp_cpu < 1:
            mp_cpu = 1
        # ----------------------------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------------------------
        # Define group of cells
        groups = get_group_cells(cells, mp_cpu)
        # Write group cells to file
        write_group_cells(groups, val_param['logs_path'])

        # Define parallel processes using pool function
        mp_pool = mp.Pool(processes=mp_cpu)
        # Iterate over cell groups
        for group, cells_array in enumerate(groups):

            # Define cell list
            cell_list = cells_array.tolist()
            # Define args
            args = zip(repeat(group), cell_list, repeat(val_param))
            # Run process(es) in multiprocessing mode
            results = mp_pool.map_async(wrap_process, args)

        mp_pool.close()
        mp_pool.join()
        #print(mp_results)
        # ----------------------------------------------------------------------------------------------------------------------

    else:

        # ----------------------------------------------------------------------------------------------------------------------
        # Iterate over cell(s)
        for group, cell in enumerate(cells):

            # ----------------------------------------------------------------------------------------------------------------------
            # Define args
            args = [group, cell, val_param]
            # Run process(es) in single mode
            results = wrap_process(args)
            # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
