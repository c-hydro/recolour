"""
Class Features

Name:          cpl_data_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230824'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

import numpy as np
import xarray as xr
import pandas as pd

from copy import deepcopy

from lib_info_args import logger_name, time_format_algorithm
from lib_data_io_cell import create_grid_cell
from lib_data_io_nc import read_file_nc, organize_file_nc

from lib_utils_time import set_time_file
from lib_utils_io import fill_path_with_tags

# set option(s)
np.set_printoptions(suppress=True)
# set logger obj
alg_logger = logging.getLogger(logger_name)

# debugging
import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class coupler dataset
class CouplerDataset:

    # method to initialize class
    def __init__(self, time_data, cell_data,
                 folder_name, file_name,
                 var_name_data=None, var_name_geo_x='longitude', var_name_geo_y='latitude',
                 time_period=0, time_frequency='D', time_rounding='D',
                 dataset_name=None, grid_ref=None,
                 dataset_template=None, time_template=None, **kwargs):

        # set class variable(s)
        self.dataset_name = dataset_name
        self.time_data = time_data
        self.time_period = time_period
        self.time_frequency = time_frequency
        self.time_rounding = time_rounding

        self.cell_data = cell_data
        self.grid_ref = grid_ref

        self.dataset_template, self.time_template = dataset_template, time_template

        self.file_name = file_name
        self.folder_name = folder_name
        self.file_path = os.path.join(self.folder_name, self.file_name)

        # variable name(s)
        self.var_name_data = var_name_data
        self.var_name_geo_x, self.var_name_geo_y = var_name_geo_x, var_name_geo_y

    # method to read data (internal to the class)
    def read(self, time_step=None):

        # method to fill the filename(s)
        file_path = fill_path_with_tags(
            self.file_path, time_step, tmpl_tags_time=self.time_template,
            tmpl_dset_obj={'cell_n': str(self.cell_data)}, tmpl_tags_dset=self.dataset_template
        )
        # read data
        dset_src = read_file_nc(file_path)

        return dset_src

    # method to organize source
    def get_data(self):

        # info start method
        alg_logger.info(' --------> Get data ... ')

        # compute time set
        time_reference = set_time_file(
            time_start=self.time_data, time_end=None,
            time_frequency=self.time_frequency, time_rounding=self.time_rounding,
            time_period=self.time_period, time_reverse=True)

        # iterate over time(s)
        dset_src_collections, dset_src_time = None, None
        for time_step in time_reference:

            # call read method
            dset_src_step = self.read(time_step)

            if dset_src_collections is None:
                if dset_src_step is not None:
                    # alg_logger.error(f' ====> No data for day {time_step.strftime("%Y/%m/%d")}. SKIPPING ... ')
                    # raise SystemExit
                    dset_src_collections = dset_src_step.copy()
                    dset_src_time = deepcopy(time_step)
            else:
                alg_logger.error(' ===> Time steps must be equal to 1')
                raise NotImplemented('Case not implemented yet')

        # info end method
        alg_logger.info(' --------> Get data ... DONE')

        return dset_src_collections, dset_src_time

    # method to organize data
    def organize_data(self, dset_src_in, dset_src_time=None):

        # info start method
        alg_logger.info(' --------> Organize data ... ')

        # organize data
        dset_src_out = organize_file_nc(
            dset_src_in, obj_time=dset_src_time, obj_variable=self.var_name_data, obj_cell=self.cell_data)

        # info end method
        alg_logger.info(' --------> Organize data ... DONE')

        return dset_src_out

    # method to organize time
    def organize_time(self, dset_src_in,
                      time_reference_dset, time_reference_group,
                      time_tolerance_period=12, time_tolerance_frequency='H'):

        # info start method
        alg_logger.info(' --------> Organize time ... ')

        # check datasets obj
        if dset_src_in is not None:

            # check time dimension availability
            if 'time' not in list(dset_src_in.dims):
                alg_logger.error(' ===> Time dimension is not available in the datasets')
                raise RuntimeError('Time dimension is needed by the algorithm to properly run')

            # order dataset (to avoid misleading in time indexes)
            dset_src_sort = dset_src_in.sortby("time", ascending=True)

            # check time sorting
            time_sort_start = pd.to_datetime(dset_src_sort['time'].values)[0]
            time_sort_end = pd.to_datetime(dset_src_sort['time'].values)[-1]

            # info datasets time period
            alg_logger.info('  (1) Datasets time period :: time_start "' +
                            time_sort_start.strftime(time_format_algorithm) + '" - time_end "' +
                            time_sort_end.strftime(time_format_algorithm) + '"')

            # check time reference
            time_reference_start = pd.date_range(
                end=time_reference_group, periods=time_tolerance_period + 1, freq=time_tolerance_frequency)[0]
            time_reference_end = pd.date_range(
                start=time_reference_group, periods=time_tolerance_period + 1, freq=time_tolerance_frequency)[-1]
            time_reference_range = pd.date_range(start=time_reference_start, end=time_reference_end,
                                                 freq=time_tolerance_frequency)

            # info search time period
            alg_logger.info('  (2) Search time period :: time_start "' +
                            time_reference_start.strftime(time_format_algorithm) + '" - time_end "' +
                            time_reference_end.strftime(time_format_algorithm) + '"')

            # define time difference between time reference and time start/end
            time_difference_start = time_reference_group - time_reference_start
            time_difference_end = time_reference_group - time_reference_end

            hours_difference_start = int(time_difference_start.days * 24) + int(time_difference_start.seconds / 3600)
            hours_difference_end = int(time_difference_end.days * 24) + int(time_difference_end.seconds / 3600)

            # select sliced times
            if time_sort_start in time_reference_range:
                time_start = time_sort_start
            else:
                time_start = time_reference_range[0]
            if time_sort_end in time_reference_range:
                time_end = time_sort_end
            else:
                time_end = time_reference_range[-1]

            # check times expected with times available
            if (time_start <= time_sort_start) or (time_sort_end <= time_end):

                # dataset selection
                dset_src_select = dset_src_sort.sel(time=slice(time_start, time_end))

                # check time selection
                if dset_src_select['time'].shape[0] >= 1:

                    # dataset last
                    dset_src_last = dset_src_select.isel(time=[-1])
                    # get last available time step
                    time_last = pd.to_datetime(dset_src_last['time'].values)[0]

                    # compute time tolerance
                    time_range_tolerance_sx = pd.date_range(
                        end=time_reference_dset, periods=time_tolerance_period + 1, freq=time_tolerance_frequency)
                    time_range_tolerance_dx = pd.date_range(
                        start=time_reference_dset, periods=time_tolerance_period + 1, freq=time_tolerance_frequency)
                    time_range_tolerance = time_range_tolerance_sx.append(time_range_tolerance_dx)

                    # check time last in time tolerance
                    if time_last in time_range_tolerance:

                        # select data and time according with time tolerance
                        dset_src_out = deepcopy(dset_src_last)
                        dset_time = deepcopy(time_last)

                        # info select time reference and time selection
                        alg_logger.info('  (3) Select time data :: time_reference "' +
                                        time_reference_group.strftime(time_format_algorithm) + '" - time_select "' +
                                        dset_time.strftime(time_format_algorithm) + '"')

                    else:
                        alg_logger.warning(' ===> Time datasets steps are not included in the time period tolerance')
                        dset_src_out, dset_time = None, None
                else:
                    alg_logger.warning(' ===> Time datasets steps are not included in the sliced times')
                    dset_src_out, dset_time = None, None
            else:
                alg_logger.warning(' ===> Time group steps are not included in the time period tolerance')
                dset_src_out, dset_time = None, None
        else:
            alg_logger.warning(' ===> Dataset obj is defined by NoneType or empty. Time information are not available')
            dset_src_out, dset_time = None, None

        # info end method
        alg_logger.info(' --------> Organize time ... DONE')

        return dset_src_out, dset_time

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class coupler ancillary
class CouplerAncillary(CouplerDataset):

    # method to initialize class
    def __init__(self, cell_data,
                 folder_name, file_name,
                 var_name_data=None, var_name_geo_x='longitude', var_name_geo_y='latitude',
                 **kwargs):

        time_data = None
        super(CouplerAncillary, self).__init__(
            time_data=time_data, cell_data=cell_data, folder_name=folder_name, file_name=file_name,
            var_name_data=var_name_data, **kwargs)

    def get_ancillary(self):

        # info start method
        alg_logger.info(' -------> Get ancillary ... ')
        # call read method
        anc_src = CouplerDataset.read(self, time_step=None)
        # info end method
        alg_logger.info(' -------> Get ancillary ... DONE')

        return anc_src

    # method to organize ancillary
    def organize_ancillary(self, dset_src_in, dset_src_time=None):

        # info start method
        alg_logger.info(' -------> Organize ancillary ... ')
        # organize ancillary
        dset_src_out = organize_file_nc(dset_src_in,
                                        obj_time=dset_src_time, obj_variable=self.var_name_data, obj_cell=self.cell_data)
        # info end method
        alg_logger.info(' -------> Organize ancillary ... DONE')

        return dset_src_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to select datasets fields
def select_datasets_fields(alg_datasets_fields, tag_variable_obj=None):

    if tag_variable_obj is None:
        tag_variable_obj = {
            'folder_name': 'folder_name', 'file_name': 'file_name',
            'variable': 'var_name_data',
            'time_rounding': 'time_rounding', 'time_frequency': 'time_frequency', 'time_period': 'time_period'
        }

    alg_datasets_args = {}
    for var_name_in, var_name_out in tag_variable_obj.items():
        if var_name_in in list(alg_datasets_fields.keys()):
            alg_datasets_args[var_name_out] = alg_datasets_fields[var_name_in]
        else:
            alg_logger.warning(' ===> Datasets fields "' + var_name_in + '" does not exist in datasets obj')
            alg_datasets_args[var_name_out] = None

    return alg_datasets_args

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get datasets geo info
def get_dataset_geo(dset_obj, var_geo_x='longitude', var_geo_y='latitude'):
    if dset_obj is not None:

        if isinstance(dset_obj, xr.Dataset):
            dset_geo_x = dset_obj[var_geo_x].values
            dset_geo_y = dset_obj[var_geo_y].values
        elif isinstance(dset_obj, dict):
            dset_geo_x = dset_obj[var_geo_x]
            dset_geo_y = dset_obj[var_geo_y]
        else:
            alg_logger.error(' ===> Geographical object format is not supported')
            raise NotImplemented('Case not implemented yet')
        dset_grid = create_grid_cell(dset_geo_x, dset_geo_y)
    else:
        alg_logger.warning(' ===> Geographical information are not available. Datasets object is defined by NoneType')
        dset_geo_x, dset_geo_y, dset_grid = None, None, None

    return dset_geo_x, dset_geo_y, dset_grid
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to find datasets lut
def get_datasets_lut(grid_ref, grid_other=None, max_dist=25000):
    if grid_other is not None:
        idxs_ref_other = grid_ref.calc_lut(other=grid_other, max_dist=max_dist)
    else:
        alg_logger.warning(' ===> LUT information are not available. Datasets object is defined by NoneType')
        idxs_ref_other = None
    return idxs_ref_other
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to filter datasets obj
def filter_datasets_obj(dset_obj_in, idxs_obj,
                        remove_no_data=False, value_no_data=-9999, fill_no_data=np.nan):

    # set variable list excluded
    var_list_excluded = ['time', 'alt', 'location_id']

    # get variable list generic
    if isinstance(dset_obj_in, xr.Dataset):
        var_list_generic = dset_obj_in.variables
    elif isinstance(dset_obj_in, dict):
        var_list_generic = list(dset_obj_in.keys())
    else:
        alg_logger.error(' ===> Object data format is not supported')
        raise NotImplemented('Case not implemented yet')

    # iterate over variable(s)
    dset_obj_out = {}
    for var_name_data in var_list_generic:
        # check if variable must be included in the output obj
        if var_name_data not in var_list_excluded:
            # get data
            if isinstance(dset_obj_in, xr.Dataset):
                var_data_tmp = dset_obj_in[var_name_data].values
                var_data_tmp = np.squeeze(var_data_tmp)
            elif isinstance(dset_obj_in, dict):
                var_data_tmp = dset_obj_in[var_name_data]
            else:
                alg_logger.error(' ===> Object data format is not supported')
                raise NotImplemented('Case not implemented yet')

            # filter data
            var_data_filter = var_data_tmp[idxs_obj]

            # check that the dimension is correct: multiple values for the same timestamp are not allowed
            if len(var_data_filter.shape) > 1:
                # more than 1 dimension: check if they are duplicates
                if np.unique(np.diff(var_data_filter)) == 0:  # this means that the array is made of duplicates
                    alg_logger.warning(
                        f' ====> More than 1 layer for single time for product {dset_obj_in.product}')
                    var_data_filter = var_data_filter[:, 0] # take only the first values
                else:
                    alg_logger.error(
                        f' ====> More than 1 layer for single time for product {dset_obj_in.product}. Cannot hash to single product.')
                    raise ValueError('More than one layer of values.')

            if remove_no_data:
                if (not isinstance(var_data_filter[0], (int, np.integer))) and (np.isnan(fill_no_data)):
                    var_data_filter[var_data_filter == value_no_data] = fill_no_data
                else:
                    pass

            # store data
            dset_obj_out[var_name_data] = var_data_filter

    return dset_obj_out
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remove datasets fields
def remove_datasets_fields(dset_obj, field_list=None):

    if field_list is None:
        field_list = ['longitude', 'latitude']

    for field_name in field_list:
        if field_name in list(dset_obj.keys()):
            dset_obj.pop(field_name)
    return dset_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to organize datasets obj
def organize_datasets_obj(dset_obj_ref, time_ref,
                          dset_obj_anc=None, dset_obj_k1=None, dset_obj_k2=None, dset_obj_composite=None,
                          var_geo_x='longitude', var_geo_y='latitude'):

    # get geographical information
    geo_x_values, geo_y_values = dset_obj_ref[var_geo_x], dset_obj_ref[var_geo_y]
    dset_geo = {var_geo_x: geo_x_values, var_geo_y: geo_y_values}

    # organize reference information
    if dset_obj_ref is not None:
        dset_obj_ref = remove_datasets_fields(dset_obj_ref, field_list=[var_geo_x, var_geo_y])
    else:
        dset_obj_ref = {}
    # organize ancillary information
    if dset_obj_anc is not None:
        dset_obj_anc = remove_datasets_fields(dset_obj_anc, field_list=[var_geo_x, var_geo_y])
    else:
        dset_obj_anc = {}
    # organize k1 information
    if dset_obj_k1 is not None:
        dset_obj_k1 = remove_datasets_fields(dset_obj_k1, field_list=[var_geo_x, var_geo_y])
    else:
        dset_obj_k1 = {}
    # organize k2 information
    if dset_obj_k2 is not None:
        dset_obj_k2 = remove_datasets_fields(dset_obj_k2, field_list=[var_geo_x, var_geo_y])
    else:
        dset_obj_k2 = {}
    # organize composite information
    if dset_obj_composite is not None:
        dset_obj_composite = remove_datasets_fields(dset_obj_composite, field_list=[var_geo_x, var_geo_y])
    else:
        dset_obj_composite = {}

    dict_collections = {**dset_geo, **dset_obj_ref, **dset_obj_k1, **dset_obj_k2, **dset_obj_composite, **dset_obj_anc}

    dframe_collections = pd.DataFrame(data=dict_collections)
    dframe_collections.attrs = {'time': time_ref}

    return dframe_collections
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to remap dataset to reference
def remap_dataset_2_ref(dset_ref, dset_other, var_other='soil_moisture_k1',
                        search_rad=25000, min_neigh=1, neigh=8, min_value=0.01, max_value=1,
                        filter_st_dev=5, filter_mode='center'):

    from repurpose.resample import resample_to_grid
    from astropy.convolution import convolve, Gaussian2DKernel

    ref_geo_x, ref_geo_y = dset_ref['longitude'].values, dset_ref['latitude'].values
    ref_geo_x_1d, ref_geo_y_1d = np.unique(ref_geo_x), np.unique(ref_geo_y)
    ref_geo_x_2d, ref_geo_y_2d = np.meshgrid(ref_geo_x_1d, ref_geo_y_1d)

    other_values_source = np.squeeze(dset_other[var_other].values)
    other_geo_x, other_geo_y = np.squeeze(dset_other['longitude'].values), np.squeeze(dset_other['latitude'].values)

    # check if values are not nan(s)
    if not all(np.isnan(other_values_source)):

        if other_values_source.ndim > 1:
            other_values_source = other_values_source[:,0]

        obj_resample = resample_to_grid(
            {'tmp': other_values_source},
            other_geo_x, other_geo_y, ref_geo_x_2d, ref_geo_y_2d,
            search_rad=search_rad,
            min_neighbours=min_neigh, neighbours=neigh)

        other_values_resampled = obj_resample['tmp']
        other_values_resampled[other_values_resampled < min_value] = np.nan
        other_values_resampled[other_values_resampled > max_value] = np.nan

        obj_kernel = Gaussian2DKernel(x_stddev=filter_st_dev, mode=filter_mode)
        other_values_filtered = convolve(other_values_resampled, obj_kernel)

        ''' debug
        plt.figure()
        plt.imshow(other_values_resampled)
        plt.colorbar(); plt.clim(0, 1)
        plt.figure()
        plt.imshow(other_values_filtered)
        plt.colorbar(); plt.clim(0, 1)
        plt.show()
        '''

    else:
        # maps other values to reference (nans values)
        other_values_filtered = np.full_like(ref_geo_x_2d, np.nan)

    obj_data = {var_other: other_values_filtered.ravel(),
                'longitude': ref_geo_x_2d.ravel(), 'latitude': ref_geo_y_2d.ravel()}

    return obj_data
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to join datasets
def join_dataset_obj(dset_collections_in, time_collections_in, ancillary_collections_in,
                     max_dist_k1=25000, max_dist_k2=25000, max_dist_anc=25000,
                     tag_name_ref='ref', tag_name_k1='k1', tag_name_k2='k2',
                     remap_k1_2_ref=True, remap_k2_2_ref=False,
                     bnd_value_min=0.0, bnd_value_max=1):

    # info start method
    alg_logger.info(' -------> Join data ... ')

    # info start organize data
    alg_logger.info(' --------> Organize data ... ')
    # define datasets
    dset_ref, time_ref = dset_collections_in[tag_name_ref], time_collections_in[tag_name_ref]
    dset_k1, time_k1 = dset_collections_in[tag_name_k1], time_collections_in[tag_name_k1]
    dset_k2, time_k2 = dset_collections_in[tag_name_k2], time_collections_in[tag_name_k2]

    # debug
    if dset_k1 is None:
        alg_logger.warning(' ===> Data "k1" defined by NoneType. Remap based on data "ref" initialize with nans')
        time_k1 = time_ref
        dset_k1 = create_dataset_null(dset_ref, 'soil_moisture_ref', 'soil_moisture_k1', np.nan)
    if dset_k2 is None:
        alg_logger.warning(' ===> Data "k2" defined by NoneType. Remap based on data "ref" initialize with nans')
        time_k2 = time_ref
        dset_k2 = create_dataset_null(dset_ref, 'soil_moisture_ref', 'soil_moisture_k2', np.nan)

    dset_anc = deepcopy(ancillary_collections_in)

    # info start organize data
    alg_logger.info(' --------> Organize data ... DONE')

    # remap datasets k1 to reference
    if remap_k1_2_ref:
        if (dset_ref is not None) and (dset_k1 is not None):
            dset_k1_remap = remap_dataset_2_ref(dset_ref, dset_k1,
                                                var_other='soil_moisture_k1', search_rad=max_dist_k1)
        else:
            dset_k1_remap = None
    else:
        dset_k1_remap = dset_k1

    # info start remap data
    alg_logger.info(' --------> Remap data ... ')

    # remap datasets k2 to reference
    if remap_k2_2_ref:
        if (dset_ref is not None) and (dset_k2 is not None):
            dset_k2_remap = remap_dataset_2_ref(dset_ref, dset_k2,
                                                var_other='soil_moisture_k2', search_rad=max_dist_k2)
        else:
            dset_k2_remap = None
    else:
        dset_k2_remap = dset_k2

    # info end remap data
    alg_logger.info(' --------> Remap data ... DONE')

    # info start combine data
    alg_logger.info(' --------> Combine data ... ')

    # check ancillary datasets (not defined by NoneType)
    if dset_anc is not None:

        # check dynamic datasets (not defined by NoneType)
        if (dset_k1_remap is not None) and (dset_k2_remap is not None):

            # define geo
            geo_x_ref_raw, geo_y_ref_raw, grid_ref = get_dataset_geo(dset_ref)
            geo_x_anc_raw, geo_y_anc_raw, grid_anc = get_dataset_geo(dset_anc)
            geo_x_k1_raw, geo_y_k1_raw, grid_k1 = get_dataset_geo(dset_k1_remap)
            geo_x_k2_raw, geo_y_k2_raw, grid_k2 = get_dataset_geo(dset_k2_remap)

            # define generic indexes
            idxs_ref_anc = get_datasets_lut(grid_ref, grid_anc, max_dist=max_dist_anc)
            idxs_ref_k1 = get_datasets_lut(grid_ref, grid_k1, max_dist=max_dist_k1)
            idxs_ref_k2 = get_datasets_lut(grid_ref, grid_k2, max_dist=max_dist_k2)

            # filter indexes (select only the ones with all the datasets available)
            index_id = np.arange(0, idxs_ref_anc.__len__(), 1)
            index_data = {'ref_k1': idxs_ref_k1, 'ref_k2': idxs_ref_k2, 'ref_anc': idxs_ref_anc}
            index_dframe = pd.DataFrame(data=index_data, index=index_id)
            index_dframe[index_dframe <= 0] = np.nan
            index_dframe.dropna(axis=0, how='any', inplace=True)
            # define filtered indexes
            idxs_ref_filter = np.array(index_dframe.index, dtype=int).tolist()
            idxs_k1_filter = np.array(index_dframe['ref_k1'], dtype=int).tolist()
            idxs_k2_filter = np.array(index_dframe['ref_k2'], dtype=int).tolist()
            idxs_anc_filter = np.array(index_dframe['ref_anc'], dtype=int).tolist()

            '''
            # define filtered indexes
            idxs_ref_filter, idxs_k1_filter, idxs_k2_filter, idxs_anc_filter = [], [], [], []
            for ref_idx, (k1_idx, k2_idx, anc_idx) in enumerate(zip(idxs_ref_k1, idxs_ref_k2, idxs_ref_anc)):
                if (k1_idx > 0) and (k2_idx > 0) and (anc_idx > 0):
                    idxs_ref_filter.append(ref_idx)
                    idxs_k1_filter.append(k1_idx)
                    idxs_k2_filter.append(k2_idx)
                    idxs_anc_filter.append(anc_idx)
            '''

            # define filtered datasets
            dset_ref_filter = filter_datasets_obj(dset_obj_in=dset_ref, idxs_obj=idxs_ref_filter)
            dset_anc_filter = filter_datasets_obj(dset_obj_in=dset_anc, idxs_obj=idxs_anc_filter)
            dset_k1_filter = filter_datasets_obj(dset_obj_in=dset_k1_remap, idxs_obj=idxs_k1_filter)
            dset_k2_filter = filter_datasets_obj(dset_obj_in=dset_k2_remap, idxs_obj=idxs_k2_filter)

            ''' test
            print('TEST POINT - CPL_DATA_DYNAMIC L381')
            sm_tmp = dset_ref_filter['soil_moisture_ref']
    
            sm_tmp[sm_tmp < 0] = np.nan
            series_tmp = pd.Series(data=sm_tmp)
            series_filled = series_tmp.fillna(method='ffill')
            sm_filled = series_filled.values
            dset_ref_filter['soil_moisture_ref'] = sm_filled
            '''

            dset_composite_filter = create_composite_dset(dset_ref_filter, dset_k1_filter, dset_k2_filter)
            # ISSUE: filtering using the min all the datasets will be different ... and some zero spot will be put in the final map
            idxs_composite_finite = np.squeeze(np.argwhere(dset_composite_filter['soil_moisture_composite'] >= bnd_value_min)).tolist()
            idxs_ref_finite = np.squeeze(np.argwhere(dset_ref_filter['soil_moisture_ref'] >= bnd_value_min)).tolist()
            idxs_k1_finite = np.squeeze(np.argwhere(dset_k1_filter['soil_moisture_k1'] >= bnd_value_min)).tolist()
            idxs_k2_finite = np.squeeze(np.argwhere(dset_k2_filter['soil_moisture_k2'] >= bnd_value_min)).tolist()

            dset_composite_finite = filter_datasets_obj(
                dset_obj_in=dset_composite_filter, idxs_obj=idxs_composite_finite, remove_no_data=True)
            dset_ref_finite = filter_datasets_obj(
                dset_obj_in=dset_ref_filter, idxs_obj=idxs_composite_finite, remove_no_data=True)
            dset_anc_finite = filter_datasets_obj(
                dset_obj_in=dset_anc_filter, idxs_obj=idxs_composite_finite, remove_no_data=False)
            dset_k1_finite = filter_datasets_obj(
                dset_obj_in=dset_k1_filter, idxs_obj=idxs_composite_finite, remove_no_data=True)
            dset_k2_finite = filter_datasets_obj(
                dset_obj_in=dset_k2_filter, idxs_obj=idxs_composite_finite, remove_no_data=True)

            dframe_collections = organize_datasets_obj(
                dset_ref_finite, time_ref,
                dset_anc_finite, dset_k1_finite, dset_k2_finite,
                dset_composite_finite)

            if 'dim' in dframe_collections.columns:
                dframe_collections = dframe_collections.drop(['dim'], axis=1)
            if 'locations' in dframe_collections.columns:
                dframe_collections = dframe_collections.drop(['locations'], axis=1)

            # info end combine data
            alg_logger.info(' --------> Combine data ... DONE')

        else:
            # case k1 and k2 are not available
            dframe_collections = None
            alg_logger.warning(' ===> Dynamic dataset is/are defined by NoneType. Data collections are not available')
            # info end combine data
            alg_logger.info(' --------> Combine data ... FAILED. Datasets k1 and/or k2 are/is not available')
    else:
        # info end combine data
        alg_logger.info(' --------> Combine data ... FAILED. Datasets metrics are not available')
        dframe_collections = None

        # case metrics are not available but datasets are available
        if (dset_k1_remap is not None) and (dset_k2_remap is not None):
            alg_logger.error(' ===> Ancillary dataset is defined by NoneType. Data collections are not available')
            raise RuntimeError('Ancillary dataset must be defined to properly run the algorithm')

    # info end method
    alg_logger.info(' -------> Join data ... DONE')

    return dframe_collections

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to composite
def create_composite_dset(dset_ref, dset_k1, dset_k2):

    dset_composite = deepcopy(dset_ref)
    n_composite = dset_composite['soil_moisture_ref'].__len__()
    cell_composite, loc_composite = dset_composite['cell'], dset_composite['locations']
    lons_composite, lats_composite = dset_composite['longitude'], dset_composite['latitude']

    values_dict = {
        'ref': dset_ref['soil_moisture_ref'],
        'k1': dset_k1['soil_moisture_k1'],
        'k2': dset_k2['soil_moisture_k2'],
    }

    values_ref = values_dict['ref']
    values_k1 = values_dict['k1']
    values_k2 = values_dict['k2']

    values_composite = np.zeros([n_composite])
    values_composite[:] = -9999.0
    flags_composite = np.zeros([n_composite])
    flags_composite[:] = -1

    idxs_ref_finite = np.squeeze(np.argwhere(values_ref >= 0)).tolist()
    idxs_k1_finite = np.squeeze(np.argwhere(values_k1 >= 0)).tolist()
    idxs_k2_finite = np.squeeze(np.argwhere(values_k2 >= 0)).tolist()

    if not isinstance(idxs_ref_finite, list):
        idxs_ref_finite = [idxs_ref_finite]
    if not isinstance(idxs_k1_finite, list):
        idxs_k1_finite = [idxs_k1_finite]
    if not isinstance(idxs_k2_finite, list):
        idxs_k2_finite = [idxs_k2_finite]

    idxs_ref_k1_diff = list(set(idxs_k1_finite).difference(idxs_ref_finite))
    idxs_ref_k2_diff = list(set(idxs_k2_finite).difference(idxs_ref_finite))
    idxs_ref_k2_k1_diff = list(set(idxs_ref_k2_diff).difference(idxs_ref_k1_diff))

    values_composite[idxs_ref_finite] = values_ref[idxs_ref_finite]
    values_composite[idxs_ref_k1_diff] = values_k1[idxs_ref_k1_diff]
    values_composite[idxs_ref_k2_k1_diff] = values_k2[idxs_ref_k2_k1_diff]

    flags_composite[idxs_ref_finite] = 0
    flags_composite[idxs_ref_k1_diff] = 1
    flags_composite[idxs_ref_k2_k1_diff] = 2

    dset_composite = {'soil_moisture_composite': values_composite,
                      'flags_composite': flags_composite,
                      'cell': cell_composite, 'locations': loc_composite,
                      'longitude': lons_composite, 'latitude': lats_composite,
                      }

    return dset_composite
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to create a dataset of nan (or null) values
def create_dataset_null(dset_ref: xr.Dataset, dset_ref_var: str, var: str, fill_value: float):
    """Creates a dataset of nan values which is a copy of dset_ref.

    Parameters
    ----------
    :param dset_ref: reference dataset
    dset_ref_datarray: slice of the dataset which is the datarray containing dset_ref_var as a key
    :param dset_ref_var: key to substitute
    :param var: variable name which is substituted to the reference one
    :param fill_value: fillvalue for the data, default np.nan

    Returns
    -------
    xarray dataset

    """
    if fill_value is None:
        fill_value = np.nan

    dset_xr = None
    if dset_ref is not None:
        dset_ref_dict = dset_ref.to_dict()
        data_ref_len = len(dset_ref_dict['data_vars'][dset_ref_var]['data'])
        dset_ref_dict['data_vars'][dset_ref_var]['data'] = [[fill_value]] * data_ref_len
        dset_ref_dict['data_vars'][var] = dset_ref_dict['data_vars'].pop(dset_ref_var)
        dset_xr = xr.Dataset.from_dict(dset_ref_dict)

    return dset_xr
# ----------------------------------------------------------------------------------------------------------------------
