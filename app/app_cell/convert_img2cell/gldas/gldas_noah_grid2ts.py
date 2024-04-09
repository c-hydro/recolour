# ----------------------------------------------------------------------------------------------------------------------
# Libraries
from os.path import join
from gldas.reshuffle import main
# ----------------------------------------------------------------------------------------------------------------------

"""
Ridefizione dei parametri con la nuova versione di GLDAS
086_L1 	SoilMoi0_10cm_inst
085_L1 	SoilTMP0_10cm_inst
065 	SWE_inst
"""

# ----------------------------------------------------------------------------------------------------------------------
# Define date(s) and variable(s)
date_from = '2020-06-01'
date_to = '2021-05-31'
var = ["SoilMoi0_10cm_inst", "SoilTMP0_10cm_inst", "SWE_inst"]

# Define base path
path_base = '/share/HSAF_Data_Validation/datasets/'
# Define product folder(s)
folder_grid = 'gldas_noah025_3h_v2.1_grid/'
folder_ts = 'gldas_noah025_3h_v2.1_ts_2021_20200601_20210531/'
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# Create path of grid and ts data
path_grid = join(path_base, folder_grid)
path_ts = join(path_base, folder_ts)

# Define args
args = [path_grid, path_ts, date_from, date_to] + var

# Call and run procedure to convert GLDAS Noah data from grid to ts
main(args)
# ----------------------------------------------------------------------------------------------------------------------
