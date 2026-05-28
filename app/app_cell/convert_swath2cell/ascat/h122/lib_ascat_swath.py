# Copyright (c) 2025, TU Wien
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of TU Wien, Department of Geodesy and Geoinformation
#      nor the names of its contributors may be used to endorse or promote
#      products derived from this software without specific prior written
#      permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL TU WIEN DEPARTMENT OF GEODESY AND
# GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import logging
from datetime import timedelta
from functools import partial

import dask
import numpy as np
import xarray as xr

from pathlib import Path

from pyresample import kd_tree
from pyresample.geometry import AreaDefinition
from pyresample.geometry import SwathDefinition

from ascat.grids import GridRegistry

from ascat.utils import get_grid_gpis
from ascat.file_handling import Filenames
from ascat.file_handling import ChronFiles
from ascat.file_handling import FilenameTemplate
import warnings

from config_info import LOGGER_NAME

registry = GridRegistry()
logger = logging.getLogger(LOGGER_NAME)


class Swath(Filenames):
    """
    Class to read and merge swath files given one or more file paths.
    """

    def _read(self,
              filename,
              generic=True,
              preprocessor=None,
              **xarray_kwargs):
        """
        Open one swath file as an xarray.Dataset and preprocess it if necessary.

        Parameters
        ----------
        filename : str
            File to read.
        generic : bool, optional
            Not yet implemented, kept to match the signature.
        preprocessor : callable, optional
            Function to preprocess the dataset after opening.
        xarray_kwargs : dict
            Additional keyword arguments passed to xarray.open_dataset.

        Returns
        -------
        ds : xarray.Dataset
            Dataset.
        """
        ds = xr.open_dataset(
            filename,
            engine="h5netcdf",
            **xarray_kwargs,
        )
        if ds["location_id"].dtype != np.int32:
            ds["location_id"] = ds["location_id"].astype(np.int32)
        if preprocessor is not None:
            ds = preprocessor(ds)

        return ds

    def read(self, parallel=False, mask_and_scale=True, **kwargs):
        """
        Read the file or a subset of it.

        Parameters
        ----------
        parallel : bool, optional
            If True, read files in parallel.
        mask_and_scale : bool, optional
            If True, mask and scale the data.
        kwargs : dict
            Additional keyword arguments passed to `Filenames.read`.

        Returns
        -------
        ds : xarray.Dataset
            Dataset.
        """

        ds, closers = super().read(
            closer_attr="_close",
            parallel=parallel,
            mask_and_scale=mask_and_scale,
            **kwargs)
        if ds is not None:
            ds.set_close(partial(super()._multi_file_closer, closers))
            return ds

    @staticmethod
    def _nbytes(ds):
        return ds.nbytes

    def _merge(self, data):
        """
        Merge datasets.

        Parameters
        ----------
        data : list of xarray.Dataset
            Datasets to merge.

        Returns
        -------
        xarray.Dataset
            Merged dataset.
        """
        if data == []:
            return None

        merged_ds = xr.concat(
            [ds for ds in data if ds is not None and len(ds["obs"]) > 0],
            dim="obs",
            combine_attrs=self.combine_attributes,
            data_vars="minimal",
            coords="minimal",
        )

        return merged_ds

    @staticmethod
    def _ensure_obs(ds):
        """
        Makes sure that the sample dimension is named `obs`.
        """
        ds = ds.cf_geom.set_sample_dimension("obs")
        return ds

    @staticmethod
    def combine_attributes(attrs_list, context):
        """
        Decides which attributes to keep when merging swath files.

        Parameters
        ----------
        attrs_list : list of dict
            List of attributes dictionaries.
        context : None
            This currently is None, but will eventually be passed information about
            the context in which this was called.
            (see https://github.com/pydata/xarray/issues/6679#issuecomment-1150946521)

        Returns
        -------
        """
        # we don't need to pass on anything from global attributes, except for these
        global_attributes_to_pass_on_merge = [
            "grid_mapping_name", "featureType"
        ]
        if "global_attributes_flag" in attrs_list[0].keys():
            attrs_list[0].pop("global_attributes_flag")
            result = {}
            for attr in global_attributes_to_pass_on_merge:

                val = attrs_list[0].get(attr, False)
                if val:
                    result[attr] = val
            return result

        variable_attrs = attrs_list

        # this code taken straight from xarray/core/merge.py
        # Replicates the functionality of "drop_conflicts"
        # but just for variable attributes.
        result = {}
        dropped_keys = set()
        for attrs in variable_attrs:
            result.update({
                key: value
                for key, value in attrs.items()
                if key not in result and key not in dropped_keys
            })
            result = {
                key: value
                for key, value in result.items()
                if key not in attrs or
                xr.core.utils.equivalent(attrs[key], value)
            }
            dropped_keys |= {key for key in attrs if key not in result}
        return result


class SwathGridFiles(ChronFiles):
    """
    Class to manage chronological swath files with a date field in the filename.
    """

    def __init__(
        self,
        root_path_in,
        root_path_out,
        fn_templ_in,
        sf_templ_in,
        grid_name,
        date_field_fmt_in="%Y%m%d*",
        date_field_fmt_out='%Y%m%d_%H00',
        cell_fn_format=None,
        cls_kwargs=None,
        err=True,
        fn_templ_out=None,
        sf_templ_out=None,
        fn_read_fmt=None,
        sf_read_fmt=None,
        fn_write_fmt=None,
        sf_write_fmt=None,
        preprocessor=None,
        postprocessor=None,
        cache_size=0,
    ):
        """
        Initialize SwathFiles class.

        Parameters
        ----------
        root_path : str
            Root path.
        fn_templ : str
            Filename template (e.g. "{date}_ascat.nc").
        sf_templ : dict, optional
            Subfolder template defined as dictionary (default: None).
        grid_name : str
            Name of the grid - must be registered in the grid registry.
        date_field_fmt : str
            Date field format (e.g. "%Y%m%d").
        cell_fn_format : str, optional
            String to use to format cell file names (e.g. "{:04d}.nc").
        cls_kwargs : dict, optional
            Class keyword arguments (default: None).
        err : bool, optional
            Set true if a file error should be re-raised instead of
            reporting a warning.
            Default: True
        fn_read_fmt : str or function, optional
            Filename format for read operation.
        sf_read_fmt : str or function, optional
            Subfolder format for read operation.
        fn_write_fmt : str or function, optional
            Filename format for write operation.
        sf_write_fmt : str or function, optional
            Subfolder format for write operation.
        preprocessor : callable, optional
            Function to preprocess datasets after opening.
        postprocessor : callable, optional
            Function to pass to the `postprocessor` argument of `ascat.cell.RaggedArrayTs.write`
            when stacking to cell files.
        cache_size : int, optional
            Number of files to keep in memory (default=0).
        """
        # first check if any files directly under root_path contain the ending (make
        # sure not to iterate through every file - just stop after the first one).
        # This allows the user to set the root path either at the place necessitated by
        # the sf_templ or directly at the level of the files. However, the user still
        # cannot set the root path anywhere else in the directory structure (e.g. within
        # a satellite but above a year). In order to choose a specific satellite, must
        # pass that as a fmt_kwarg
        ending = fn_templ_in.split(".")[-1]
        for f in Path(root_path_in).glob(f"*.{ending}"):
            if f.is_file():
                sf_templ = None
                sf_read_fmt = None
                break

        super().__init__(root_path_in, Swath, fn_templ_in, sf_templ_in, cls_kwargs, err,
                         fn_read_fmt, sf_read_fmt, fn_write_fmt, sf_write_fmt,
                         cache_size)

        self.date_field_fmt_in = date_field_fmt_in
        self.date_field_fmt_out = date_field_fmt_out
        self.grid_name = grid_name
        self.grid = registry.get(grid_name)

        self.cell_fn_format = cell_fn_format
        self.preprocessor = preprocessor
        self.postprocessor = postprocessor

        self.root_path_out = root_path_out
        self.fn_templ_out = fn_templ_out
        self.sf_templ_out = sf_templ_out

    @classmethod
    def from_product_id(
            cls,
            path_in, filename_in,
            path_out, filename_out,
            product_id,
            date_fmt_field_in='%Y%m%d%H%M%S', date_field_fmt_out='%Y%m%d',
            cell_fn_format='{:4d}', grid_name='fibgrid_6.25',
            **kwargs
    ):

        from lib_ascat_product_info import swath_io_catalog

        logger.info(f" ----> Select product {product_id} ... ")

        product_id = product_id.upper()
        if product_id in swath_io_catalog:

            default_class = swath_io_catalog[product_id]

            product_class = default_class(
                path_in, filename_in, path_out, filename_out,
                date_fmt_field_in=date_fmt_field_in, date_field_fmt_out=date_field_fmt_out,
                cell_fn_format=cell_fn_format, grid_name=grid_name)

        else:
            error_str = f"Product {product_id} not recognized. Valid products are"
            error_str += f" {', '.join(swath_io_catalog.keys())}."
            raise ValueError(error_str)

        # get root path in
        path_in = product_class.root_path_in
        # get root path in
        path_out = product_class.root_path_out

        logger.info(f" ----> Select product {product_id} ... DONE")

        return cls.from_product_class(path_in, path_out, product_class)

    @classmethod
    def from_product_class(
        cls,
        path_in, path_out,
        product_class,
    ):

        return cls(
            root_path_in=path_in, root_path_out=path_out,
            fn_templ_in=product_class.fn_pattern_in, sf_templ_in=product_class.sf_pattern_in,
            fn_templ_out=product_class.fn_pattern_out, sf_templ_out=product_class.sf_pattern_out,
            grid_name=product_class.grid_name, cell_fn_format=product_class.cell_fn_format,
            date_field_fmt_in=product_class.date_field_fmt_in,
            date_field_fmt_out=product_class.date_field_fmt_out,
            fn_read_fmt=product_class.fn_read_fmt, sf_read_fmt=product_class.sf_read_fmt,
            fn_write_fmt=product_class.fn_write_fmt, sf_write_fmt=product_class.sf_write_fmt,
            preprocessor=product_class.preprocess_,
        )

    def _spatial_filter(
        self,
        filenames,
        cell=None,
        location_id=None,
        coords=None,
        bbox=None,
        geom=None,
    ):
        """
        Filter a search result for cells matching a spatial criterion.

        Parameters
        ----------
        cell : int or list of int
            Grid cell number to read.
        location_id : int or list of int
            Location id.
        coords : tuple of numeric or tuple of iterable of numeric
            Tuple of (lon, lat) coordinates.
        bbox : tuple
            Tuple of (latmin, latmax, lonmin, lonmax) coordinates.

        Returns
        -------
        filenames : list of str
            Filenames.
        """

        if cell is not None:
            gpis = get_grid_gpis(self.grid, cell=cell)
            spatial = SwathDefinition(
                lats=self.grid.arrlat[gpis],
                lons=self.grid.arrlon[gpis],
            )
        elif location_id is not None:
            gpis = get_grid_gpis(self.grid, location_id=location_id)
            spatial = SwathDefinition(
                lats=self.grid.arrlat[gpis],
                lons=self.grid.arrlon[gpis],
            )
        elif coords is not None:
            spatial = SwathDefinition(
                lats=[coords[1]],
                lons=[coords[0]],
            )
        elif (bbox or geom) is not None:
            if bbox is not None:
                # AreaDefinition expects (lonmin, latmin, lonmax, latmax)
                # but bbox is (latmin, latmax, lonmin, lonmax)
                bbox = (bbox[2], bbox[0], bbox[3], bbox[1])
            else:
                # If we get a geometry just take its bounding box and check
                # that intersection.
                #
                # shapely.geometry.bounds is already in the correct order
                bbox = geom.bounds
            spatial = AreaDefinition(
                "bbox",
                "",
                "EPSG:4326",
                {
                    "proj": "latlong",
                    "datum": "WGS84"
                },
                1000,
                1000,
                bbox,
            )
        else:
            spatial = None

        if spatial is None:
            return filenames

        filtered_filenames = []
        for filename in filenames:
            lazy_result = dask.delayed(self._check_intersection)(filename,
                                                                 spatial)
            filtered_filenames.append(lazy_result)

        def none_filter(fname_list):
            return [l for l in fname_list if l is not None]

        filtered_filenames = dask.delayed(none_filter)(
            filtered_filenames).compute()

        return filtered_filenames

    def _check_intersection(self, filename, spatial):
        """
        Check if a file intersects with a pyresample SwathDefinition or AreaDefinition.

        Parameters
        ----------
        filename : str
            Filename.
        gpis : list of int
            List of gpis.

        Returns
        -------
        bool
            True if the file intersects with the gpis.
        """
        f = self.cls(filename)
        ds = f.read()
        with f.read() as ds:
            lons, lats = ds["longitude"].values, ds["latitude"].values
            swath_def = SwathDefinition(lats=lats, lons=lons)
            n_info = kd_tree.get_neighbour_info(
                swath_def,
                spatial,
                radius_of_influence=15000,
                neighbours=1,
            )
            valid_input_index, _, _ = n_info[:3]
        if np.any(valid_input_index):
            return filename
        return None


    def cell_search(
            self,
            dt_run,
            date_field_fmt="%Y%m%d_%H%M",
            date_field="date", fmt_kwargs={},):

        date_field_fmt = self.date_field_fmt_out

        file_tmpl = FilenameTemplate(
            root_path=self.root_path_out, fn_templ=self.fn_templ_out, sf_templ=self.sf_templ_out)

        _, _, fn_write_fmt, sf_write_fmt = self._fmt(dt_run, **fmt_kwargs)
        fn_write_fmt[date_field] = dt_run.strftime(date_field_fmt)

        filenames = file_tmpl.build_filename(fn_fmt=fn_write_fmt, sf_fmt=sf_write_fmt)

        return filenames

    def swath_search(
        self,
        dt_start,
        dt_end,
        dt_delta=None,
        search_date_fmt="%Y%m%d*",
        date_field="date",
        end_inclusive=True,
        cell=None,
        location_id=None,
        coords=None,
        bbox=None,
        geom=None,
        **fmt_kwargs,
    ):
        """
        Search for swath files within a time range and spatial criterion.

        Parameters
        ----------
        dt_start : datetime
            Start date.
        dt_end : datetime
            End date.
        dt_delta : timedelta
            Time delta.
        search_date_fmt : str
            Search date format.
        date_field : str
            Date field.
        end_inclusive : bool
            End date inclusive.
        cell : int or list of int
            Grid cell number to read.
        location_id : int or list of int
            Location id.
        coords : tuple of numeric or tuple of iterable of numeric
            Tuple of (lon, lat) coordinates.
        bbox : tuple
            Tuple of (latmin, latmax, lonmin, lonmax) coordinates.
        geom : shapely.geometry
            Geometry.
        fmt_kwargs : dict
            Additional keyword arguments passed to ascat.file_handling.ChronFiles.search_period.

        Returns
        -------
        list of str
            Filenames.
        """
        dt_delta = dt_delta or timedelta(days=1)

        filenames = self.search_period(
            dt_start,
            dt_end,
            dt_delta,
            search_date_fmt,
            date_field,
            date_field_fmt=self.date_field_fmt_in,
            end_inclusive=end_inclusive,
            **fmt_kwargs,
        )

        filtered_filenames = self._spatial_filter(
            filenames,
            cell=cell,
            location_id=location_id,
            coords=coords,
            bbox=bbox,
            geom=geom,
        )

        if len(filtered_filenames) == 0:
            raise FileNotFoundError(
                f"No files found for {dt_start} to {dt_end} with the given spatial criteria."
            )

        return filtered_filenames

    def read(
        self,
        date_range,
        dt_delta=None,
        search_date_fmt="%Y%m%d*",
        date_field="date",
        end_inclusive=True,
        cell=None,
        location_id=None,
        coords=None,
        max_coord_dist=None,
        bbox=None,
        geom=None,
        read_kwargs=None,
        **fmt_kwargs,
    ):
        """
        Extract data from swath files within a time range and spatial criterion.

        Parameters
        ----------
        date_range : tuple of datetime.datetime
            Start and end date.
        dt_delta : timedelta
            Time delta.
        search_date_fmt : str
            Search date format.
        date_field : str
            Date field.
        end_inclusive : bool
            If True (default), include data from the end date in the result. Otherwise,
            exclude it.
        cell : int or list of int
            Grid cell number to read.
        location_id : int or list of int
            Location id to read.
        coords : tuple of numeric or tuple of iterable of numeric
            Tuple of (lon, lat) coordinates to read.
        max_coord_dist : float
            Maximum distance in meters to search for grid points near the given
            coordinates. If None, the default is np.inf.
        bbox : tuple
            Tuple of (latmin, latmax, lonmin, lonmax) coordinates to bound the data.
        geom : shapely.geometry
            Geometry to bound the data.

        Returns
        -------
        xarray.Dataset
            Dataset.
        """
        dt_start, dt_end = date_range
        filenames = self.swath_search(
            dt_start,
            dt_end,
            dt_delta,
            search_date_fmt,
            date_field,
            end_inclusive,
            **fmt_kwargs,
        )

        date_range = (np.datetime64(dt_start), np.datetime64(dt_end))

        read_kwargs = read_kwargs or {}

        def filter_ds_spatial(ds):
            if cell is not None:
                ds = ds.pgg.sel_cells(cell)
            elif location_id is not None:
                ds = ds.pgg.sel_gpis(location_id, gpi_var="location_id")
            elif coords is not None:
                ds = ds.pgg.sel_coords(coords, max_coord_dist=max_coord_dist)
            elif bbox is not None:
                ds = ds.pgg.sel_bbox(bbox)
            elif geom is not None:
                ds = ds.pgg.sel_geom(geom)
            return ds

        def preprocessor(ds):
            if self.preprocessor is not None:
                ds = self.preprocessor(ds)
            ds = filter_ds_spatial(ds)
            return ds

        read_kwargs["preprocessor"] = preprocessor

        data = self.cls(filenames).read(**read_kwargs)

        if data:
            if date_range is not None:
                mask = (data["time"] >= date_range[0]) & (
                    data["time"] <= date_range[1])
                data = data.sel(obs=mask.compute())

            return data
        warning_str = (
            "No data found for specified criteria, returning None:\n"
            f"date_range={date_range}\n"
            f"cell={cell}, location_id={location_id}, coords={coords}, bbox={bbox},\n"
            f"geom={geom}, max_coord_dist={max_coord_dist}, \n")
        warnings.warn(warning_str, UserWarning, 2)

    def stack_to_cell_files(
            self,
            dt_run,
            max_nbytes_mb,
            date_range=None,
            fmt_kwargs=None,
            cells=None,
            print_progress=True,
            parallel=True,
            **kwargs
    ):
        """
        Stack all swath files to cell files, writing them in parallel.

        Parameters
        ----------
        out_dir : str
            Output directory.
        max_nbytes : int
            Maximum number of bytes to open as xarray datasets before dumping to disk.
        date_range : tuple of datetime.datetime, optional
            Start and end date for the search.
        fmt_kwargs : dict, optional
            Additional keyword arguments passed to ascat.file_handling.ChronFiles.search_period.
        cells : list of int, optional
            List of grid cell numbers to read. If None (default), all cells are read.
        print_progress : bool, optional
            If True (default), print progress bars.
        parallel: bool, optional
            If True, write data to files in parallel (use all available resources).
        """
        from ascat.cell import RaggedArrayTs

        logger.info(f" ----> Stack datasets to cell files ... ")

        max_nbytes = max_nbytes_mb * 1024**2

        # manage filenames in
        fmt_kwargs = fmt_kwargs or {}
        if date_range is not None:
            dt_start, dt_end = date_range
            filenames_in = self.swath_search(
                dt_start, dt_end, cell=cells, **fmt_kwargs)
        else:
            filenames_in = list(Path(self.root_path).glob("**/*.nc"))

        swath = self.cls(filenames_in)

        # manage filenames out
        filenames_out = self.cell_search(dt_run)

        # iterate over periods
        for ds in swath.iter_read_nbytes(
                max_nbytes,
                preprocessor=self.preprocessor,
                print_progress=print_progress,
                chunks=-1):
            ds_cells = self.grid.gpi2cell(ds["location_id"])
            if isinstance(ds_cells, np.ma.MaskedArray):
                ds_cells = ds_cells.compressed()
            ds_cells = xr.DataArray(ds_cells, dims="obs", name="cell")

            # sorting here enables us to manually select each cell's data much faster
            # than using a .groupby
            ds = ds.sortby(ds_cells)

            unique_cells, cell_counts = np.unique(ds_cells, return_counts=True)
            cell_counts = np.hstack([0, np.cumsum(cell_counts)])

            # for each cell in unique cells, isel the slice from the dataarray corresponding to it
            ds_list = []
            cell_fnames = []
            for i, c in enumerate(unique_cells):
                if (cells is None) or (c in cells):
                    cell_ds = ds.isel(
                        obs=slice(cell_counts[i], cell_counts[i + 1]))
                    if len(cell_ds) == 0:
                        continue
                    ds_list.append(cell_ds)

                    # format only if placeholder exists
                    if '{' in filenames_out and '}' in filenames_out:
                        cell_fname = filenames_out.format(c)
                    else:
                        cell_fname = filenames_out
                    # ensure folder creation
                    Path(cell_fname).parent.mkdir(parents=True, exist_ok=True)

                    # save cell filenames (related to the list of celle selected or over the domain)
                    cell_fnames.append(cell_fname)

            writer_class = RaggedArrayTs(cell_fnames)
            writer_class.write(
                ds_list,
                parallel=parallel,
                postprocessor=self.postprocessor,
                ra_type="point",
                mode="a",
                print_progress=print_progress)



        if print_progress:
            print("\n")

        logger.info(f" ----> Stack datasets to cell files ... DONE")