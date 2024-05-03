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
import warnings
import xarray as xr
import numpy as np
import os

from gldas.interface import GLDAS_Noah_v21_025Img

from pygeobase.io_base import ImageBase, MultiTemporalImageBase
from pygeobase.object_base import Image
from pynetcf.time_series import GriddedNcOrthoMultiTs

from datetime import timedelta

# from lib_grid_gldas import cell_grid
from pygeogrids.netcdf import load_grid

from gldas.grid import GLDAS025Cellgrid
# from netCDF4 import Dataset

# debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
class GLDAS_Noah_v2_025Img(ImageBase):
    """
    Class for reading one GLDAS Noah v2.1 nc file in 0.25 deg grid.
    """

    def __init__( self,
                  filename, mode="r", parameter="SoilMoi0_10cm_inst",
                  grid_path=None, subgrid=None, array_1D=False):
        """
        Parameters
        ----------
        filename: str
            filename of the GLDAS nc file
        mode: string, optional
            mode of opening the file, only 'r' is implemented at the moment
        parameter : string or list, optional
            one or list of parameters to read, see GLDAS v2.1 documentation
            for more information (default: 'SoilMoi0_10cm_inst').
        subgrid : Cell Grid
            Subgrid of the global GLDAS Grid to use for reading image data (e.g only land points)
        array_1D: boolean, optional
            if set then the data is read into 1D arrays.
            Needed for some legacy code.
        """

        super(GLDAS_Noah_v2_025Img, self).__init__(filename, mode=mode)

        if type(parameter) != list:
            parameter = [parameter]

        self.grid_path = grid_path
        self.parameters = parameter
        self.fill_values = np.repeat(-9999.0, 1440 * 120)
        self.grid = GLDAS025Cellgrid() if not subgrid else subgrid
        self.array_1D = array_1D

    def read(self, timestamp=None):

        # print 'read file: %s' %self.filename
        # Returns the selected parameters for a gldas image and
        # according metadata

        return_img, return_metadata = {}, {}
        try:
            # dataset_test = Dataset(self.filename)
            dataset = xr.open_dataset(self.filename)  # add to use xarray library
        except IOError:
            raise IOError(f"Error opening file {self.filename}")

        param_names = []
        for parameter in self.parameters:
            param_names.append(parameter)

        # check the lats (not needed)

        # lat_1d, lon_1d = dataset.coords['lat'], dataset.coords['lon']
        # lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)
        # lat_max, lat_min = lat_2d[0, 0], lat_2d[-1, 0]

        # if lat_max < lat_min:
        #    dataset.coords['lat'] = dataset.coords['lat'][::-1]
        #    dataset = dataset.sortby(dataset['lat'])

        for parameter, variable in dataset.variables.items():
            if parameter in param_names:
                param_metadata, param_data = {}, {}
                # for attrname in variable.ncattrs():
                for attrname in variable.attrs:  # add to use xarray library
                    if attrname in ["long_name", "units"]:
                        #param_metadata.update({str(attrname): getattr(variable, attrname)})
                        param_metadata.update({str(attrname): variable.attrs[attrname]})

                param_data = dataset.variables[parameter][:]
                param_data = param_data.values  # add to use xarray library
                param_data[np.isnan(param_data)] = -9999.0

                if param_data.shape.__len__() == 3:
                    param_data = param_data[0, :, :]
                elif param_data.shape.__len__() == 4:
                    param_data = param_data[0, 0, :, :]
                else:
                    raise IOError('bad dimensions definition')

                ''' debug
                lon_1d, lat_1d = dataset.variables['lon'], dataset.variables['lat']
                lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d)

                import matplotlib.pylab as plt
                plt.figure()
                plt.imshow(param_data)
                plt.colorbar()
                plt.figure()
                plt.imshow(lat_2d)
                plt.colorbar()
                plt.show()
                '''

                np.ma.set_fill_value(param_data, -9999)

                param_data = np.concatenate((param_data.flatten(), self.fill_values))

                return_img.update({str(parameter): param_data[self.grid.activegpis]})
                return_metadata.update({str(parameter): param_metadata})

                # Check for corrupt files
                try:
                    return_img[parameter]
                except KeyError:
                    path, thefile = os.path.split(self.filename)
                    print(
                        "%s in %s is corrupt - filling"
                        "image with NaN values" % (parameter, thefile)
                    )
                    return_img[parameter] = np.empty(self.grid.n_gpi).fill(
                        np.nan
                    )

                    return_metadata["corrupt_parameters"].append()

        dataset.close()

        if self.array_1D:
            return Image(
                self.grid.activearrlon,
                self.grid.activearrlat,
                return_img,
                return_metadata,
                timestamp,
            )
        else:
            raise NotImplementedError('array 2d condition not implemented yet')

    def write(self, data):
        raise NotImplementedError()

    def flush(self):
        pass

    def close(self):
        pass
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to wrap image reader
class WrapImage(GLDAS_Noah_v2_025Img):

    def __init__(
        self,
        filename, mode="r",
        grid_path=None, parameter="SoilMoi0_10cm_inst",
        subgrid=None, array_1D=False,
    ):

        super(WrapImage, self).__init__(
            filename=filename, mode=mode,
            grid_path=grid_path, parameter=parameter,
            subgrid=subgrid, array_1D=array_1D,
        )
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
            "grid_path": grid_path
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
            WrapImage,
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
# class to execute gldas ts reader
class gldas_ts(GriddedNcOrthoMultiTs):

    def __init__(self, ts_path, grid_path=None, **kwargs):

        if grid_path is None:
            pass
            # grid_path = os.path.join(ts_path, "grid.nc")

        # grid = load_grid(grid_path)
        grid = GLDAS025Cellgrid()
        super(gldas_ts, self).__init__(ts_path, grid, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------
