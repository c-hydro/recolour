"""
Library Features:

Name:          lib_interface_soilgrids
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260521'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# libraries
import os
import logging
import rasterio
import numpy as np

from pygeobase.io_base import ImageBase
from pygeobase.object_base import Image

from lib_grid_hmc import cell_grid
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to handle image
class BaseImage(ImageBase):
    """
    Class for reading one static SoilGrids file.
    """

    def __init__(
            self,
            filename,
            mode="r",
            parameter="soil_variable",
            grid_path=None,
            subgrid=None,
            array_1D=True):

        super(BaseImage, self).__init__(filename, mode=mode)

        if type(parameter) != list:
            parameter = [parameter]

        self.parameters = parameter
        self.grid = cell_grid(grid_path=grid_path) if not subgrid else subgrid
        self.array_1D = array_1D

    def read(self, timestamp=None):

        return_img = {}
        return_metadata = {}

        try:
            dataset = rasterio.open(self.filename)
            data_raw = dataset.read(1).astype(np.float32)

            if dataset.nodata is not None:
                data_raw[data_raw == dataset.nodata] = np.nan

        except IOError:
            logging.error(' ===> Open "' + self.filename + '" failed for unknown reason.')
            raise IOError('Error in opening file')

        for parameter in self.parameters:

            param_metadata = {
                "long_name": parameter,
                "units": "-"
            }

            param_data = np.ma.masked_invalid(data_raw)
            np.ma.set_fill_value(param_data, 9999.0)

            param_data = np.ma.getdata(
                param_data.filled()
            ).flatten()

            return_img.update({
                str(parameter): param_data[self.grid.activegpis]
            })

            return_metadata.update({
                str(parameter): param_metadata
            })

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
                return_img[key] = return_img[key].reshape(
                    self.grid.activearrlat.shape
                )

            return Image(
                self.grid.activearrlon,
                self.grid.activearrlat,
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
            filename,
            mode="r",
            grid_path=None,
            parameter="soil_variable",
            subgrid=None,
            array_1D=True):

        super(WrapImage, self).__init__(
            filename=filename,
            mode=mode,
            grid_path=grid_path,
            parameter=parameter,
            subgrid=subgrid,
            array_1D=array_1D,
        )
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class to set soilgrids static image reader
class soilgrids_ds(object):
    """
    Class for reading static SoilGrids maps in GeoTIFF format.
    """

    def __init__(
            self,
            data_path,
            file_name_tmpl="{variable}.tif",
            parameter=None,
            grid_path=None,
            subgrid=None,
            array_1D=True):

        if parameter is None:
            parameter = []

        if type(parameter) != list:
            parameter = [parameter]

        self.data_path = data_path
        self.file_name_tmpl = file_name_tmpl
        self.parameters = parameter
        self.grid_path = grid_path
        self.subgrid = subgrid
        self.array_1D = array_1D

    def read(self, variable_name):

        file_name = self.file_name_tmpl.format(
            variable=variable_name
        )

        file_path = os.path.join(
            self.data_path,
            file_name
        )

        image_obj = WrapImage(
            filename=file_path,
            grid_path=self.grid_path,
            parameter=variable_name,
            subgrid=self.subgrid,
            array_1D=self.array_1D
        )

        return image_obj.read(timestamp=None)
# -------------------------------------------------------------------------------------