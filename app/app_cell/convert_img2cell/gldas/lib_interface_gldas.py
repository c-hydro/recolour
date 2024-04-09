# -*- coding: utf-8 -*-

"""
Library Features:

Name:           lib_interface_gldas
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20230711'
Version:        '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np
import os

from gldas.interface import GLDAS_Noah_v21_025Img

from pygeobase.io_base import MultiTemporalImageBase
from pynetcf.time_series import GriddedNcOrthoMultiTs

from datetime import timedelta

from lib_grid_gldas import cell_grid
from pygeogrids.netcdf import load_grid

# debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class for reading a dataset of gldas images
class gldas_ds(MultiTemporalImageBase):
    """
    Class for reading cci images in nc format.

    Parameters
    ----------
    data_path : string
        path to the nc files
    parameter : string or list, optional
        one or list of parameters to read
        Default : 'SoilMoi0_10cm_inst'
    array_1D: boolean, optional
        if set then the data is read into 1D arrays. Needed for some legacy code.
    """

    def __init__(self,
                 data_path, data_sub_path=None,
                 product='gldas_v2.1',
                 file_name_tmpl=None,
                 datetime_format=None,
                 parameter='SoilMoi0_10cm_inst',
                 grid_path=None, subgrid=None, array_1D=False,
                 start_step=0, end_step=0,
                 ):

        ioclass_kws = {
            'parameter': parameter,
            'array_1D': array_1D,
            'subgrid': subgrid,
            #"grid_path": grid_path
        }

        if product == 'gldas_v2.1':
            if data_sub_path is None:
                data_sub_path = ["%Y", "%j"]
            if file_name_tmpl is None:
                file_name_tmpl = "GLDAS_NOAH025_3H*.A{datetime}.*.nc4"
            if datetime_format is None:
                datetime_format = "%Y%m%d.%H%M"
        else:
            logging.error(' ===> Product name is not supported')
            raise NotImplemented('Case not implemented yet')

        if '{datetime}' in file_name_tmpl:
            dtime_placeholder = 'datetime'
        elif '{datetime_source}' in file_name_tmpl:
            dtime_placeholder = 'datetime_source'
        else:
            logging.error(' ===> Datetime placeholder format is not supported')
            raise NotImplemented('Case not implemented yet')

        self.start_step = start_step
        self.end_step = end_step
        self.product = product

        super(gldas_ds, self).__init__(
            data_path,
            GLDAS_Noah_v21_025Img,
            fname_templ=file_name_tmpl,
            datetime_format=datetime_format,
            subpath_templ=data_sub_path,
            exact_templ=False,
            dtime_placeholder=dtime_placeholder,
            ioclass_kws=ioclass_kws,
        )

    def tstamps_for_daterange(self, start_date, end_date):
        """
        return timestamps for daterange,

        Parameters
        ----------
        start_date: datetime
            start of date range
        end_date: datetime
            end of date range

        Returns
        -------
        timestamps : list
            list of datetime objects of each available image between
            start_date and end_date
        """
        img_offsets = np.array(
            [
                timedelta(hours=0),
                timedelta(hours=3),
                timedelta(hours=6),
                timedelta(hours=9),
                timedelta(hours=12),
                timedelta(hours=15),
                timedelta(hours=18),
                timedelta(hours=21),
            ]
        )

        timestamps = []
        diff = end_date - start_date
        for i in range(diff.days + 1):
            daily_dates = start_date + timedelta(days=i) + img_offsets
            timestamps.extend(daily_dates.tolist())

        return timestamps
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to execute cci ts reader
class gldas_ts(GriddedNcOrthoMultiTs):

    def __init__(self, ts_path, grid_path=None, **kwargs):

        if grid_path is None:
            grid_path = os.path.join(ts_path, "grid.nc")

        grid = load_grid(grid_path)
        super(gldas_ts, self).__init__(ts_path, grid, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------
