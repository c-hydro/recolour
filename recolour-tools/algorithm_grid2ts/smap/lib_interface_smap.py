"""
Library Features:

Name:           lib_interface_smap
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org),
                Martina Natali (martina01.natali@edu.unife.it)
Date:           '20230711'
Version:        '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import logging
import rasterio
import numpy as np
import os

from pygeobase.io_base import ImageBase, MultiTemporalImageBase
from pygeobase.object_base import Image
from pynetcf.time_series import GriddedNcOrthoMultiTs

from datetime import timedelta

from lib_grid_smap import cell_grid
from pygeogrids.netcdf import load_grid
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to execute smap image reader
class BaseImage(ImageBase):
    """
    Class for reading one smap file in 0.005 deg grid.
    """

    def __init__(self,
                 filename, mode="r", parameter="soil_moisture",
                 grid_path=None, subgrid=None, array_1D=False,):
        """
        Parameters
        ----------
        filename: str
            filename of the SMAP tiff file
        mode: string, optional
            mode of opening the file, only 'r' is implemented at the moment
        parameter : string or list, optional
            one or list of parameters to read, see SMAP documentation
            for more information (default: 'soil_moisture').
        subgrid : Cell Grid
            Subgrid of the global SMAP Grid to use for reading image data (e.g only land points)
        array_1D: boolean, optional
            if set then the data is read into 1D arrays.
            Needed for some legacy code.
        """

        super(BaseImage, self).__init__(filename, mode=mode)

        if type(parameter) != list:
            parameter = [parameter]

        self.parameters = parameter
        # self.fill_values = np.repeat(9999.0, 1440 * 120)
        self.fill_values = np.repeat(9999.0, 0)
        self.grid = cell_grid(grid_path=grid_path) if not subgrid else subgrid
        self.array_1D = array_1D

    def read(self, timestamp=None):

        return_img, return_metadata = {}, {}
        try:
            dataset = rasterio.open(self.filename)
            data_raw = dataset.read()
            data_values = np.float64(data_raw[0, :, :])

        except IOError:
            logging.error(' ===> Open "' + self.filename + '" failed for unknown reason.')
            raise IOError('Error in opening file')

        param_metadata, param_data = {}, {}
        for parameter in self.parameters:
            if parameter == 'soil_moisture':

                param_metadata['long_name'] = parameter
                param_metadata['units'] = '-'

                param_data = np.ma.masked_array(data_values)

                np.ma.set_fill_value(param_data, 9999)
                param_data = np.concatenate((self.fill_values, np.ma.getdata(param_data.filled()).flatten(),))

                return_img.update({str(parameter): param_data[self.grid.activegpis]})
                return_metadata.update({str(parameter): param_metadata})

                # check the corrupted file
                try:
                    return_img[parameter]
                except KeyError:

                    folder, file = os.path.split(self.filename)
                    logging.warning(" ===> %s in %s is corrupt - filling image with NaN values" % (parameter, file))

                    return_img[parameter] = np.empty(self.grid.n_gpi).fill(np.nan)
                    return_metadata["corrupt_parameters"].append()

            else:
                logging.error(' ===> Parameter "' + parameter + '" is not supported')
                raise NotImplemented('Parameter not implemented yet. Only "soil_moisture" parameter is supported')

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
# class to wrap smap image reader
class WrapImage(BaseImage):

    def __init__(
        self,
        filename, mode="r",
        grid_path=None, parameter="soil_moisture",
        subgrid=None, array_1D=False,
    ):

        super(WrapImage, self).__init__(
            filename=filename, mode=mode,
            grid_path=grid_path, parameter=parameter,
            subgrid=subgrid, array_1D=array_1D,
        )
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to set smap image reader
class smap_ds(MultiTemporalImageBase):
    """
    Class for reading SMAP images in tiff format.

    Parameters
    ----------
    data_path : string
        Path to the nc files
    parameter : string or list, optional
        one or list of parameters to read, see GLDAS v2.1 documentation
        for more information (default: 'SoilMoi0_10cm_inst').
    subgrid : Cell Grid
        Subgrid of the global GLDAS Grid to use for reading image data (e.g only land points)
    array_1D: boolean, optional
        If set then the data is read into 1D arrays.
        Needed for some legacy code.
    """

    def __init__(self,
                 data_path, data_sub_path=None,
                 product='smap',
                 file_name_tmpl=None, datetime_tmpl=None,
                 parameter=None,
                 grid_path=None, subgrid=None, array_1D=False,
                 start_step=6, end_step=18,
                 ):

        ioclass_kws = {
            "parameter": parameter,
            "subgrid": subgrid,
            "array_1D": array_1D,
            "grid_path": grid_path
        }

        if data_sub_path is None:
            data_sub_path = ['%Y', '%m', '%d']
        if file_name_tmpl is None:
            file_name_tmpl = 'soil_moisture_{datetime}.tiff'
        if datetime_tmpl is None:
            datetime_tmpl = "%Y%m%d_%H%M"

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

        super(smap_ds, self).__init__(
            data_path,
            WrapImage,
            fname_templ=file_name_tmpl,
            datetime_format=datetime_tmpl,
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

        start_step = self.start_step
        end_step = self.end_step
        delta = end_step - start_step

        img_offsets = np.array([timedelta(hours=step) for step in range(start_step, end_step +1, delta)])

        timestamps = []
        diff = end_date - start_date
        for i in range(diff.days + 1):
            daily_dates = start_date + timedelta(days=i) + img_offsets
            timestamps.extend(daily_dates.tolist())

        return timestamps
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to execute smap ts reader
class smap_ts(GriddedNcOrthoMultiTs):

    def __init__(self, ts_path, grid_path=None, **kwargs):
        """
        Class for reading SMAP time series after reshuffling.

        Parameters
        ----------
        ts_path : str
            Directory where the netcdf time series files are stored
        grid_path : str, optional (default: None)
            Path to grid file, that is used to organize the location of time
            series to read. If None is passed, grid.nc is searched for in the
            ts_path.

        Optional keyword arguments that are passed to the Gridded Base:
        ------------------------------------------------------------------------
            parameters : list, optional (default: None)
                Specific variable names to read, if None are selected, all are read.
            offsets : dict, optional (default:None)
                Offsets (values) that are added to the parameters (keys)
            scale_factors : dict, optional (default:None)
                Offset (value) that the parameters (key) is multiplied with
            ioclass_kws: dict
                Optional keyword arguments to pass to OrthoMultiTs class:
                ----------------------------------------------------------------
                    read_bulk : boolean, optional (default:False)
                        if set to True the data of all locations is read into memory,
                        and subsequent calls to read_ts read from the cache and not from disk
                        this makes reading complete files faster#
                    read_dates : boolean, optional (default:False)
                        if false dates will not be read automatically but only on specific
                        request useable for bulk reading because currently the netCDF
                        num2date routine is very slow for big datasets
        """

        if grid_path is None:
            grid_path = os.path.join(ts_path, "grid.nc")

        grid = load_grid(grid_path)
        super(smap_ts, self).__init__(ts_path, grid, **kwargs)
# -------------------------------------------------------------------------------------
