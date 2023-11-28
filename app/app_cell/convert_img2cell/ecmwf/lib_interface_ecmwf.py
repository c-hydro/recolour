# -*- coding: utf-8 -*-

"""
Library Features:

Name:           lib_interface_ecmwf
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230711'
Version:        '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import warnings
import numpy as np
import xarray as xr
import os

try:
    import pygrib
except ImportError:
    warnings.warn("pygrib has not been imported")

from pygeobase.io_base import ImageBase, MultiTemporalImageBase
from pygeobase.object_base import Image
from pynetcf.time_series import GriddedNcOrthoMultiTs

from datetime import timedelta

from lib_grid_ecmwf import cell_grid
from netCDF4 import Dataset
from pygeogrids.netcdf import load_grid

# debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to handle image
class BaseImage(ImageBase):
    """
    Class for reading one ECMWF file.

    Parameters
    ----------
    filename: string
        filename of the ECMWF nc file
    mode: string, optional
        mode of opening the file, only 'r' is implemented at the moment
    parameter : string or list, optional
        one or list of parameters to read, see ASCAT documentation for more information
        Default : 'var40'
    array_1D: boolean, optional
        if set then the data is read into 1D arrays. Needed for some legacy code.
    """
    
    def __init__(self,
                 filename, mode='r', parameter='var40',
                 grid_path=None, subgrid=None, array_1D=False):

        super(BaseImage, self).__init__(filename, mode=mode)
        
        if type(parameter) != list:
            parameter = [parameter]

        self.grid_path = grid_path
        self.parameters = parameter
        self.fill_values = np.repeat(9999., 0)
        self.grid = cell_grid(self.grid_path) if not subgrid else subgrid
        self.array_1D = array_1D

    def read(self, timestamp=None):
        
        # returns the selected parameters for a gldas image and according metadata
        return_img, return_metadata = {}, {}
        
        try:
            # dataset = Dataset(self.filename)
            dataset = xr.open_dataset(self.filename) # add to use xarray library
        except IOError as e:
            print(e)
            print(" ".join([self.filename, "can not be opened"]))
            raise e

        param_names = []
        for parameter in self.parameters:
            param_names.append(parameter)

        # check the lons because are from 0 to 360; according to the rest of
        # algorithm at this point the dataset will be order from 0-180 and -180-0
        # the grid loaded in the init must be from -180-0 and 0-180
        lon_max = np.nanmax(dataset.coords['lon'])
        if lon_max > 180:
            dataset.coords['lon'] = (dataset.coords['lon'] + 180) % 360 - 180 # correct!
            dataset = dataset.sortby(dataset['lon'])

        for parameter, variable in dataset.variables.items():
            if parameter in param_names:
                param_metadata, param_data = {}, {}
                # for attrname in variable.ncattrs():
                for attrname in variable.attrs: # add to use xarray library
                    if attrname in ['long_name', 'units']:
                        param_metadata.update({str(attrname): getattr(variable, attrname)})

                param_data = dataset.variables[parameter][:]
                param_data = param_data.values  # add to use xarray library
                param_data[np.isnan(param_data)] = -9999.0 # add to use xarray library

                if param_data.shape.__len__() == 3:
                    param_data = param_data[0, :, :]
                elif param_data.shape.__len__() == 4:
                    param_data = param_data[0, 0, :, :]
                else:
                    raise IOError('bad dimensions definition')
                np.ma.set_fill_value(param_data, 9999)

                ''' debug
                data_tmp = param_data.copy()
                data_tmp[data_tmp < 0] = np.nan
                plt.figure()
                plt.imshow(data_tmp)
                plt.colorbar()
                plt.show()
                '''

                param_data = np.concatenate((self.fill_values, param_data.flatten()))

                return_img.update({str(parameter): param_data[self.grid.activegpis]})
                return_metadata.update({str(parameter): param_metadata})
                        
            # check for corrupt files
                try:
                    return_img[parameter]
                except KeyError:

                    folder, file = os.path.split(self.filename)
                    logging.warning(" ===> %s in %s is corrupt - filling image with NaN values" % (parameter, file))

                    return_img[parameter] = np.empty(self.grid.n_gpi).fill(np.nan)
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
            for key in return_img:
                return_img[key] = np.flipud(
                    return_img[key].reshape((720, 1440))
                )

            return Image(
                np.flipud(self.grid.activearrlon.reshape((720, 1440))),
                np.flipud(self.grid.activearrlat.reshape((720, 1440))),
                return_img,
                return_metadata,
                timestamp,
            )

    def write(self):
        raise NotImplementedError()
        
    def flush(self):
        pass
    
    def close(self):
        pass
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to wrap image reader
class WrapImage(BaseImage):

    def __init__(
        self,
        filename, mode="r",
        grid_path=None, parameter="var40",
        subgrid=None, array_1D=False,
    ):

        super(WrapImage, self).__init__(
            filename=filename, mode=mode,
            grid_path=grid_path,parameter=parameter,
            subgrid=subgrid, array_1D=array_1D,
        )
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# class for reading a dataset of ECMWF images
class ecmwf_ds(MultiTemporalImageBase):
    """
    Class for reading ECMWF images in nc format.

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
                 product='h14',
                 file_name_tmpl=None,
                 datetime_format=None,
                 parameter='var40',
                 grid_path=None, subgrid=None, array_1D=False,
                 start_step=0, end_step=0,
                 ):

        ioclass_kws = {
            'parameter': parameter,
            'array_1D': array_1D,
            'subgrid': subgrid,
            "grid_path": grid_path
        }

        if product == 'h14':
            if data_sub_path is None:
                data_sub_path = ['%Y', '%m', '%d']
            if file_name_tmpl is None:
                file_name_tmpl = "rzsm_{datetime}_h14.nc"
            if datetime_format is None:
                datetime_format = "%Y%m%d00"
        elif product == 'h26':
            if data_sub_path is None:
                data_sub_path = []
            if file_name_tmpl is None:
                file_name_tmpl = "h26_{datetime}_R01.nc"
            if datetime_format is None:
                datetime_format = "%Y%m%d00"
        elif product == 'h141':
            if data_sub_path is None:
                data_sub_path = ['%Y', '%m', '%d']
            if file_name_tmpl is None:
                file_name_tmpl = "h141_{datetime}_R01.nc"
            if datetime_format is None:
                datetime_format = "%Y%m%d00"
        elif product == 'h142':
            if data_sub_path is None:
                data_sub_path = ['%Y', '%m', '%d']
            if file_name_tmpl is None:
                file_name_tmpl = "h142_{datetime}_R01.nc"
            if datetime_format is None:
                datetime_format = "%Y%m%d00"
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

        super(ecmwf_ds, self).__init__(
            data_path, WrapImage,
            fname_templ=file_name_tmpl,
            datetime_format=datetime_format,
            subpath_templ=data_sub_path,
            exact_templ=True,
            dtime_placeholder=dtime_placeholder,
            ioclass_kws=ioclass_kws
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
        img_offsets = np.array([timedelta(hours=0)])

        timestamps = []
        diff = end_date - start_date
        for i in range(diff.days + 1):
            daily_dates = start_date + timedelta(days=i) + img_offsets
            timestamps.extend(daily_dates.tolist())

        return timestamps
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to execute ecmwf ts reader
class ecmwf_ts(GriddedNcOrthoMultiTs):

    def __init__(self, ts_path, grid_path=None, **kwargs):

        if grid_path is None:
            grid_path = os.path.join(ts_path, "grid.nc")

        grid = load_grid(grid_path)
        super(ecmwf_ts, self).__init__(ts_path, grid, **kwargs)
# ----------------------------------------------------------------------------------------------------------------------
