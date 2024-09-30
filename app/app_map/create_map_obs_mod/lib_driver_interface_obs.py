"""
Library Features:

Name:          lib_driver_interface_obs
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240209'
Version:       '2.0.8'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import glob
import re

from ascat.read_native.cdr import AscatGriddedNcTs as AscatNc
#from ascat.h_saf import AscatNrtBufrFileList as AscatL2SsmBufrChunked

from pygeogrids.grids import CellGrid
from pynetcf.time_series import IndexedRaggedTs, GriddedNcTs

from lib_info_args import logger_name
from lib_driver_interface_hsaf import AscatNrtBufrFileList, AscatEpsBufrFileList

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to interface ascat h16 swath files
class AscatDriverH16_OLD(AscatNrtBufrFileList):
    """
    Parameters
    ----------
    path: string
        path where the data is stored
    month_path_str: string, optional
        if the files are stored in folders by month as is the standard on
        the H SAF FTP Server then please specify the string that should be
        used in datetime.datetime.strftime
        Default: 'h16_%Y%m_buf'
    """

    def __init__(self, path,
                 file_search_str='h16_{datetime}*.buf',
                 month_path_str='h16_%Y%m_buf', chunk_minutes=50):

        day_search_str = 'h16_%Y%m%d_*.buf'
        datetime_format = '%Y%m%d_%H%M%S'
        filename_datetime_format = (4, 19, '%Y%m%d_%H%M%S')

        super(AscatDriverH16, self).__init__(
            root_path=path,
            product_id='h16',
            filename_template='{product_id}_{date}*.buf',
            subfolder_template='%Y/%m/%d/%HH')


        '''
        super(AscatDriverH16, self).__init__(path,
                                             month_path_str=month_path_str,
                                             day_search_str=day_search_str,
                                             file_search_str=file_search_str,
                                             datetime_format=datetime_format,
                                             filename_datetime_format=filename_datetime_format,
                                             chunk_minutes=chunk_minutes)
        '''

    # ------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to interface ascat h101 swath files
class AscatDriverH101_OLD(AscatNrtBufrFileList):
    """
    Parameters
    ----------
    path: string
        path where the data is stored
    month_path_str: string, optional
        if the files are stored in folders by month as is the standard on
        the H SAF FTP Server then please specify the string that should be
        used in datetime.datetime.strftime
        Default: 'h101_%Y%m_buf'
    """

    def __init__(self, path,
                 file_search_str='h101_{datetime}*.buf',
                 month_path_str='h101_%Y%m_buf',
                 chunk_minutes=50):

        day_search_str = 'h101_%Y%m%d_*.buf'
        datetime_format = '%Y%m%d_%H%M%S'
        filename_datetime_format = (5, 20, '%Y%m%d_%H%M%S')
        '''
        super(AscatDriverH101, self).__init__(path,
                                              month_path_str=month_path_str,
                                              day_search_str=day_search_str,
                                              file_search_str=file_search_str,
                                              datetime_format=datetime_format,
                                              filename_datetime_format=filename_datetime_format,
                                              chunk_minutes=chunk_minutes)
        '''
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to interface ascat data record files
class AscatDriverSsmCdr(AscatNc):

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

    def __init__(self, cdr_path, grid_path, cdr_tmpl_filename=None,
                 grid_filename='TUW_WARP5_grid_info_2_2.nc',
                 static_layer_path=None, **kwargs):

        first_file = glob.glob(os.path.join(cdr_path, '*.nc'))[0]

        root_path, root_file = os.path.split(first_file)
        root_name, root_ext  = os.path.splitext(root_file)
        root_cell = re.findall(r'\d\d\d\d', root_file)[0]
        fn_format = root_name.replace(root_cell, '{:04d}')

        # version = os.path.basename(first_file).rsplit('_', 1)[0]
        # fn_format = '{:}_{{:04d}}'.format(version)

        grid_filename = os.path.join(grid_path, grid_filename)

        super(AscatDriverSsmCdr, self).__init__(cdr_path, fn_format,
                                                grid_filename,
                                                static_layer_path, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Class to write data in a cell indexed file
class AscatDriverIndexed(GriddedNcTs):

    def __init__(self, *args, **kwargs):
        kwargs['ioclass'] = IndexedRaggedTs
        super(AscatDriverIndexed, self).__init__(*args, **kwargs)

    def write_cell(self, cell, gpi, data, datefield='jd', filename_tmpl='%04d.nc', filename_upd=True):
        """
        Write complete data set into cell file.

        Parameters
        ----------
        cell : int
            Cell number.
        gpi : numpy.ndarray
            Location ids.
        data : dict or numpy record array
            dictionary with variable names as keys and numpy.arrays as values
        datefield: string
            field in the data dict that contains dates in correct format
        """
        if isinstance(self.grid, CellGrid) is False:
            raise TypeError("Associated grid is not of type "
                            "pygeogrids.CellGrid.")

        if datefield is None:
            log_stream.warning(' ===> Datefield not set, using default "jd"')
            datefield = 'jd'

        if self.mode != 'w':
            log_stream.error(' ===> File not opened in write mode.')
            raise ValueError("File not opened in write mode.")

        cell_all = self.grid.gpi2cell(gpi)
        cell_unique = np.unique(cell_all)

        select_idx = np.argwhere(cell_all == cell)[:, 0]
        select_gpis = gpi[select_idx]
        select_data = data[select_idx]

        idx_all = np.argwhere(self.grid.arrcell == cell)[:, 0]
        tmp_cell = np.unique(self.grid.arrcell[idx_all])

        select_lons, select_lats = self.grid.gpi2lonlat(select_gpis)

        if type(select_data) != dict:
            select_data = {key: select_data[key] for key in select_data.dtype.names}

        select_times = select_data[datefield]
        del select_data[datefield]

        if tmp_cell.size > 1 or tmp_cell != cell:
            log_stream.error(' ===> Cell is not unique or not equal to the requested cell')
            raise ValueError('Cell is not unique or not equal to the requested cell')

        filename_def = filename_tmpl.format(cell=cell)
        filename = os.path.join(self.path, filename_def)

        if filename_upd:
            if os.path.exists(filename):
                os.remove(filename)

        if os.path.isfile(filename):
            mode = 'a'
        else:
            mode = 'w'

        if self.previous_cell != cell:
            self.flush()
            self.close()
            self.previous_cell = cell
            if self.mode == 'w':
                if 'n_loc' not in self.ioclass_kws:
                    n_loc = self.grid.grid_points_for_cell(cell)[0].size
                    self.ioclass_kws['n_loc'] = n_loc
            self.fid = self.ioclass(filename, mode=mode,
                                    **self.ioclass_kws)
            self.ioclass_kws.pop('n_loc', None)

        self.fid.write_ts(select_gpis, select_data, select_times, lon=select_lons, lat=select_lats, dates_direct=True)

        self.flush()
        self.close()
# ----------------------------------------------------------------------------------------------------------------------

    def write_cell_OLD(self, cell, gpi, data, datefield, filename_tmpl='%04d.nc', filename_upd=True):
        """
        Write complete data set into cell file.

        Parameters
        ----------
        cell : int
            Cell number.
        gpi : numpy.ndarray
            Location ids.
        data : dict or numpy record array
            dictionary with variable names as keys and numpy.arrays as values
        datefield: string
            field in the data dict that contains dates in correct format
        """
        if isinstance(self.grid, CellGrid) is False:
            raise TypeError("Associated grid is not of type "
                            "pygeogrids.CellGrid.")

        if self.mode != 'w':
            raise ValueError("File not opened in write mode.")

        cell_all = self.grid.gpi2cell(gpi)
        cell_unique = np.unique(cell_all)

        idx_all = np.argwhere(cell_all == cell)[:, 0]

        tmp_cell = np.unique(self.grid.arrcell[gpi])

        lons_all = self.grid.arrlon[gpi]
        lats_all = self.grid.arrlat[gpi]
        gpis_all = self.grid.gpis[gpi].data

        cell_all = self.grid.gpi2cell(gpis_all)

        if tmp_cell.size > 1 or tmp_cell != cell:

            select_idx = np.argwhere(cell_all == cell)[:, 0]

            gpis_select = gpi[select_idx]
            data_select = data[select_idx]
            lons_select = lons_all[select_idx]
            lats_select = lats_all[select_idx]

        else:
            gpis_select = gpi
            data_select = data
            lons_select = lons_all
            lats_select = lats_all

        if type(data) != dict:
            data_select = {key: data_select[key] for key in data_select.dtype.names}

        dates_select = data_select[datefield]
        del data_select[datefield]

        # lons = self.grid.arrlon[gpi]
        # lats = self.grid.arrlat[gpi]
        # gpis = self.grid.gpis[gpi]
        # gpis = gpis.data

        filename_def = filename_tmpl.format(cell=cell)
        filename = os.path.join(self.path, filename_def)

        if filename_upd:
            if os.path.exists(filename):
                os.remove(filename)

        if os.path.isfile(filename):
            mode = 'a'
        else:
            mode = 'w'

        if self.previous_cell != cell:
            self.flush()
            self.close()
            self.previous_cell = cell
            if self.mode == 'w':
                if 'n_loc' not in self.ioclass_kws:
                    n_loc = self.grid.grid_points_for_cell(cell)[0].size
                    self.ioclass_kws['n_loc'] = n_loc
            self.fid = self.ioclass(filename, mode=mode,
                                    **self.ioclass_kws)
            self.ioclass_kws.pop('n_loc', None)

        self.fid.write_ts(gpis_select, data_select, dates_select, lon=lons_select, lat=lats_select,
                          dates_direct=True)

        self.flush()
        self.close()
# ----------------------------------------------------------------------------------------------------------------------