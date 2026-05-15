
from resampler import main

# Settings (H16, H101, H102, H103)
product = 'H16'

path_swath = '/share/HSAF_SM/ascat/nrt/h16/'
path_ts = '/share/VALIDATION_SM_HSAF/cell/tmp_idx/h16/'

path_grid = '/share/VALIDATION_SM_HSAF/auxiliary/reference/ascat/grid/'
filename_grid = 'TUW_WARP5_grid_info_2_3.nc'

date_from = '2023-06-01T00:00'
date_to = '2025-06-01T00:00'

# HSAF H16 method to convert data from image to ts
main([product, path_swath, path_ts, path_grid, date_from, date_to])






