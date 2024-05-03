"""
Library Features:

Name:          lib_interface_ascat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import glob
from ascat.read_native.cdr import AscatGriddedNcTs
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to interface ascat data record
class AscatCdr(AscatGriddedNcTs):

    """
    Class reading Metop ASCAT soil moisture Climate Data Record (CDR).

    Parameters
    ----------
    cdr_path : str
        Path to Climate Data Record (CDR) data set.
    grid_path : str
        Path to grid file.
    grid_filename : str
        Name of grid file.
    static_layer_path : str
        Path to static layer files.

    Attributes
    ----------
    grid : pygeogrids.CellGrid
        Cell grid.
    """

    def __init__(self, cdr_path, grid_path,
                 grid_filename='TUW_WARP5_grid_info_2_3.nc',
                 static_layer_path=None, **kwargs):

        list_file = glob.glob(os.path.join(cdr_path, '*.nc'))
        if list_file.__len__() == 0:
            logging.error(' ===> ASCAT files not found in the folder "' + cdr_path + '"')
            raise FileNotFoundError('File must be available in the selected folder')

        first_file = list_file[0]
        version = os.path.basename(first_file).rsplit('_', 1)[0]
        fn_format = '{:}_{{:04d}}'.format(version)
        grid_filename = os.path.join(grid_path, grid_filename)

        super(AscatCdr, self).__init__(cdr_path, fn_format, grid_filename, static_layer_path, **kwargs)

# -------------------------------------------------------------------------------------
