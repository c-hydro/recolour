# -------------------------------------------------------------------------------------
# Libraries
import logging
import os
import re
import time

from pytesmo.time_series.filters import exp_filter
from copy import deepcopy
from datetime import datetime

from multiprocessing import Pool, cpu_count

import numpy as np
import xarray as xr
import pandas as pd

from lib_info_args import time_format_algorithm
from lib_data_io_generic import create_darray_3d, create_dset, write_dset, write_obj, read_obj
from lib_data_zip_gzip import unzip_filename, zip_filename
from lib_utils_io import create_filename_tmp
from lib_utils_system import make_folder

import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to write model file
def write_file_model(file_path, var_data, var_ancillary,
                     terrain_data, terrain_geo_x, terrain_geo_y, terrain_attrs,
                     file_compression_flag=True, file_compression_ext='.gz', file_compression_level=9):

    if 'time' in list(var_ancillary.keys()):
        var_time = pd.DatetimeIndex([var_ancillary['time']])
    if 'west_east' in list(var_ancillary.keys()):
        var_geo_x = var_ancillary['west_east']
    if 'south_north' in list(var_ancillary.keys()):
        var_geo_y = var_ancillary['south_north']

    file_dataset = None
    for var_key, var_values in var_data.items():

        var_dataset = create_dset(var_values, terrain_data, terrain_geo_x, terrain_geo_y, var_time,
                                  coord_name_x='Longitude', coord_name_y='Latitude',
                                  var_data_name=var_key, var_geo_active=False)

        if file_dataset is None:
            file_dataset = deepcopy(var_dataset)
        else:
            file_dataset = xr.combine_by_coords([file_dataset, var_dataset])

    if terrain_attrs is not None:
        file_dataset.attrs = terrain_attrs

    folder_name, file_name = os.path.split(file_path)
    make_folder(folder_name)

    write_dset(file_path, file_dataset,
               dset_mode='w', dset_engine='h5netcdf', dset_compression=file_compression_level, dset_format='NETCDF4',
               dim_key_time='time', fill_data=-9999.0)

    if file_compression_flag:
        file_path_unzip = file_path
        file_path_zip = file_path + file_compression_ext
        zip_filename(file_path_unzip, file_path_zip, file_compression_level=file_compression_level)

        if os.path.exists(file_path_unzip):
            os.remove(file_path_unzip)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read model file
def read_file_model(file_name_obj, folder_tmp=None,
                    file_idxs=None,
                    var_file_ancillary='data.workspace', var_upd_ancillary=True,
                    var_name_in='SM', var_name_out='soil_moisture',
                    var_index_in='Index', var_index_out='index',
                    var_geo_x_in='Longitude', var_geo_x_out='longitude',
                    var_geo_y_in='Latitude', var_geo_y_out='latitude',
                    var_time_in='time', var_time_out='time', var_no_data=-9999.0):

    logging.info(' -----> Read datasets from model file(s) ... ')

    if isinstance(var_name_in, list):
        var_name_in = var_name_in[0]
    if isinstance(var_name_out, list):
        var_name_out = var_name_out[0]

    if var_upd_ancillary:
        if os.path.exists(var_file_ancillary):
            os.remove(var_file_ancillary)

    if not os.path.exists(var_file_ancillary):

        var_geo_idxs = file_idxs.ravel()

        var_ancillary = None
        var_time = []
        var_collections = np.zeros([var_geo_idxs.shape[0], file_name_obj.__len__()])
        var_collections[:, :] = np.nan
        for file_id, file_name_step in enumerate(file_name_obj):

            time_match = re.search(r'\d{4}\d{2}\d{2}\d{2}\d{2}', file_name_step)
            time_string = datetime.strptime(time_match.group(), '%Y%m%d%H%M')
            time_analysis = pd.Timestamp(time_string)

            logging.info(' ------> Analysis Time: ' + time_analysis.strftime(format=time_format_algorithm) + ' ... ')

            var_time.append(time_analysis)
            if os.path.exists(file_name_step):

                if file_name_step.endswith('.gz'):
                    make_folder(folder_tmp)
                    file_name_def = create_filename_tmp(prefix='tmp_', suffix='.nc', folder=folder_tmp)
                    unzip_filename(file_name_step, file_name_def)
                else:
                    file_name_def = file_name_step

                file_handle = xr.open_dataset(file_name_def, decode_times=False)
                file_geo_x = file_handle[var_geo_x_in].values
                file_geo_y = file_handle[var_geo_y_in].values
                file_time = file_handle[var_time_in].values[0]

                if var_ancillary is None:
                    geo_y_upper = file_geo_y[0, 0]
                    geo_y_lower = file_geo_y[-1, 0]
                    if geo_y_lower > geo_y_upper:
                        file_geo_y = np.flipud(file_geo_y)

                    var_geo_x = file_geo_x.ravel()
                    var_geo_y = file_geo_y.ravel()

                    var_ancillary = {var_geo_x_out: var_geo_x, var_geo_y_out: var_geo_y, var_index_out: var_geo_idxs}

                # if not isinstance(file_time, str):
                #    file_time = pd.Timestamp(time_string)

                if var_name_in not in [var_geo_x_in, var_geo_y_in, var_time_in, 'crs', 'times']:

                    if var_name_in in list(file_handle.keys()):
                        if geo_y_lower > geo_y_upper:
                            var_values = np.flipud(file_handle[var_name_in].values)
                        else:
                            var_values = file_handle[var_name_in]
                    var_values[var_values == var_no_data] = np.nan
                else:
                    var_values = np.zeros([1, file_name_obj.__len__()])
                    var_values[:] = np.nan

                if file_name_step != file_name_def:
                    os.remove(file_name_def)

                logging.info(' ------> Analysis Time: ' + time_analysis.strftime(format=time_format_algorithm) + ' ... DONE')

            else:
                var_values = np.zeros([1, file_name_obj.__len__()])
                var_values[:] = np.nan

                logging.info(' ------> Analysis Time: ' + time_analysis.strftime(format=time_format_algorithm)  +
                             ' ... SKIPPED. File "' + file_name_step + '" does not exists')

            # Store values in a collections variable
            var_collections[:, file_id] = var_values.ravel()

        variable_n = var_collections.shape[0]
        range_n = list(range(0, variable_n, 1))

        variable_dset = xr.Dataset(coords={'time': (['time'], pd.DatetimeIndex(var_time))})
        variable_dset[var_name_out] = xr.DataArray(
            var_collections, name=var_name_out,  dims=['n', 'time'],
            coords={'n': (['n'], range_n), 'time': (['time'], pd.DatetimeIndex(var_time))})
        variable_dset.attrs = var_ancillary

        folder_name, file_name = os.path.split(var_file_ancillary)
        make_folder(folder_name)

        write_obj(var_file_ancillary, variable_dset)

        logging.info(' -----> Read datasets from model file(s) ... DONE')

    else:

        logging.info(' -----> Read datasets from model file(s) ... SKIPPED. Datasets already read and imported')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to parse name datasets
def parse_name_datasets(var_name, var_split='_'):

    var_parts = var_name.split(var_split)

    var_root_folder = var_split.join([var_parts[0], var_parts[1]])
    var_sub_folder = var_parts[2]
    var_name = var_split.join([var_parts[0], var_parts[1], var_parts[2], var_parts[3]])

    return var_name, var_root_folder, var_sub_folder
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to dump datasets
def save_fx_datasets(results):
    global fx_dset_common

    attrs_expected = ['id', 'tag', 'longitude', 'latitude']

    file_dset = results[0]

    if fx_dset_common is None:
        fx_dset_common = deepcopy(file_dset)
    else:
        if 'var_list' in list(file_dset.attrs):
            var_list = file_dset.attrs['var_list']
        else:
            var_list = ['soil_moisture']

        time_common_index = fx_dset_common['time'].values

        variable_dict = {}
        for var_step in var_list:
            var_common_step = fx_dset_common[var_step].values
            var_file_step = file_dset[var_step].values
            # variable_dict[var_step] = (["n","time"], np.vstack([var_common_step, var_file_step]))
            variable_dict[var_step] = np.vstack([var_common_step, var_file_step])

        dict_attrs_upd = {}
        for attr_key, attr_value_file in file_dset.attrs.items():

            if attr_key in attrs_expected:
                if attr_key in list(fx_dset_common.attrs.keys()):
                    attr_value_common = fx_dset_common.attrs[attr_key]
                    if not isinstance(attr_value_common, list):
                        attr_value_common = [attr_value_common]
                    attr_value_common.append(attr_value_file)

                    dict_attrs_upd[attr_key] = attr_value_common

        variable_n = variable_dict['soil_moisture'].shape[0]
        range_n = list(range(0, variable_n, 1))

        fx_dset_upd = xr.Dataset(coords={'time': (['time'], time_common_index)})
        for variable_key, variable_data in variable_dict.items():

            variable_da = xr.DataArray(
                variable_data, name=variable_key,  dims=['n', 'time'],
                coords={'n': (['n'], range_n), 'time': (['time'], time_common_index)})

            fx_dset_upd[variable_key] = variable_da
        fx_dset_upd.attrs = dict_attrs_upd

        fx_dset_common = deepcopy(fx_dset_upd)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to merge fx datasets
def merge_fx_datasets(var_name, var_file_ancillary, var_file_collections,
                      file_chunks=None, format_chunks_element='{0:08d}', format_chunk_group='{0:05d}',
                      tag_chunks_geo_x='longitude', tag_chunks_geo_y='latitude', tag_chunks_id='id'):

    logging.info(' -----> Merge datasets ... ')

    id_chunks = file_chunks[tag_chunks_id]
    geo_x_chunks = file_chunks[tag_chunks_geo_x]
    geo_y_chunks = file_chunks[tag_chunks_geo_y]

    n_chunks = id_chunks.__len__()
    n_chunks_str = format_chunks_element.format(n_chunks)
    for var_name_step, file_path_step in zip(var_name, var_file_collections):

        logging.info(' ------> Analyze variable "' + var_name_step + '" ... ')

        if not os.path.exists(file_path_step):

            pnt_time = None
            pnt_values_collections = None
            pnt_attrs_collections = None
            for i_chunks, (id_part, geo_x_part, geo_y_part) in enumerate(zip(id_chunks, geo_x_chunks, geo_y_chunks)):

                i_chunks_str = format_chunk_group.format(i_chunks)

                logging.info(' -------> Analyze datasets "' + var_name_step +
                             '" for chunk "' + i_chunks_str + '/' + n_chunks_str + '" ... ')

                file_pnt_ancillary = var_file_ancillary.format(
                    point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)

                if os.path.exists(file_pnt_ancillary):

                    pnt_dframe = read_obj(file_pnt_ancillary)
                    pnt_attrs = pnt_dframe.attrs

                    if pnt_time is None:
                        pnt_time = pnt_dframe['time'].values

                    pnt_attr_list = list(pnt_attrs.keys())
                    pnt_var_list = list(pnt_dframe.data_vars)

                    if var_name_step in pnt_var_list:
                        var_name_dframe = var_name_step
                    else:
                        if var_name_step == 'soil_moisture':
                            var_name_dframe = 'sm'
                        else:
                            raise IOError('Variable ' + var_name_step + ' is not defined in the datasets')

                    if pnt_attrs_collections is None:
                        pnt_attrs_collections = {}

                    for pnt_attr_step in pnt_attr_list:
                        pnt_attr_values_step = pnt_attrs[pnt_attr_step]

                        if pnt_attr_step not in list(pnt_attrs_collections.keys()):
                            pnt_attrs_collections[pnt_attr_step] = pnt_attr_values_step
                        else:
                            pnt_attr_values_tmp = pnt_attrs_collections[pnt_attr_step]
                            pnt_attr_values_hstack = pnt_attr_values_tmp + pnt_attr_values_step
                            pnt_attrs_collections[pnt_attr_step] = pnt_attr_values_hstack

                    pnt_var_values_step = pnt_dframe[var_name_dframe].values

                    if pnt_values_collections is None:
                        pnt_values_collections = {}

                    if var_name_step not in list(pnt_values_collections.keys()):
                        pnt_values_collections[var_name_step] = pnt_var_values_step
                    else:
                        pnt_var_values_tmp = pnt_values_collections[var_name_step]
                        pnt_var_values_vstack = np.vstack((pnt_var_values_tmp, pnt_var_values_step))
                        pnt_values_collections[var_name_step] = pnt_var_values_vstack

                else:
                    logging.error(' ===> Error in reading file "' + file_pnt_ancillary + '". Merge datasets failed.')
                    raise IOError('File not found')

                logging.info(' -------> Analyze datasets "' + var_name_step +
                             '" for chunk "' + i_chunks_str + '/' + n_chunks_str + '" ... DONE')

            variable_n = pnt_values_collections[var_name_step].shape[0]
            range_n = list(range(0, variable_n, 1))

            variable_dset = xr.Dataset(coords={'time': (['time'], pnt_time)})
            for variable_key, variable_data in pnt_values_collections.items():
                variable_da = xr.DataArray(
                    variable_data, name=variable_key, dims=['n', 'time'],
                    coords={'n': (['n'], range_n), 'time': (['time'], pnt_time)})
                variable_dset[variable_key] = variable_da
            variable_dset.attrs = pnt_attrs_collections

            folder_name, file_name = os.path.split(file_path_step)
            make_folder(folder_name)

            write_obj(file_path_step, variable_dset)

            logging.info(' ------> Analyze variable "' + var_name_step + '" ... DONE.')

        else:
            logging.info(' ------> Analyze variable "' + var_name_step + '" ... SKIPPED. Filename "' +
                         file_path_step + '" is already available.')

    logging.info(' -----> Compute datasets ... DONE')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to wrap fx datasets
def wrap_fx_datasets(var_data, var_time, var_file_ancillary, var_upd_ancillary=True,
                     var_name_in='soil_moisture', var_name_out='soil_moisture',
                     var_filter_list=None, format_filter_list='swi_t{0:02d}',
                     cpu_n=2, cpu_max=None,
                     file_chunks=None, format_chunks_element='{0:08d}', format_chunk_group='{0:05d}',
                     tag_chunks_geo_x='longitude', tag_chunks_geo_y='latitude', tag_chunks_id='id'):

    global var_data_common
    global var_time_common
    global fx_dset_common

    logging.info(' -----> Compute datasets ... ')

    if var_filter_list is None:
        var_filter_list = [6, 12, 32]

    if cpu_max is None:
        cpu_max = cpu_count() - 1
    if cpu_n > cpu_max:
        logging.warning(' ===> Maximum of recommended processes must be less then ' + str(cpu_max))
        logging.warning(' ===> Set number of process from ' + str(cpu_n) + ' to ' + str(cpu_max))
        #cpu_n = cpu_max

    id_chunks = file_chunks[tag_chunks_id]
    geo_x_chunks = file_chunks[tag_chunks_geo_x]
    geo_y_chunks = file_chunks[tag_chunks_geo_y]

    var_datetime = var_time.to_pydatetime()
    var_time_length = var_datetime.shape[0]

    n_chunks = id_chunks.__len__()
    n_chunks_str = format_chunks_element.format(n_chunks)

    var_data_common = var_data
    var_time_common = var_time

    id_ts = 0
    for i_chunks, (id_part, geo_x_part, geo_y_part) in enumerate(zip(id_chunks, geo_x_chunks, geo_y_chunks)):

        i_chunks_str = format_chunk_group.format(i_chunks)

        logging.info(' ------> Analyze datasets for chunk ' + i_chunks_str + '/' + n_chunks_str + ' ... ')
        time_start = time.time()

        file_pnt_ancillary = var_file_ancillary.format(
            point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)

        if var_upd_ancillary:
            if os.path.exists(file_pnt_ancillary):
                os.remove(file_pnt_ancillary)

        if not os.path.exists(file_pnt_ancillary):

            fx_dset_common = None
            exec_pool = Pool(cpu_n)
            for i_pnt, (id_pnt, geo_x_pnt, geo_y_pnt) in enumerate(zip(id_part, geo_x_part, geo_y_part)):

                # id_pnt_name, id_pnt_root_folder, id_pnt_sub_folder = parse_name_datasets(id_pnt)

                fx_args = {'id': id_ts, 'tag': id_pnt, 'geo_x': geo_x_pnt, 'geo_y': geo_y_pnt,
                           'var_name': var_name_out,
                           'file_name': file_pnt_ancillary, 'flag_updating': True,
                           'filter_list': var_filter_list, 'filter_name': format_filter_list}

                # SERIAL
                # file_dset = exec_fx_datasets(fx_args)
                # save_fx_datasets(file_dset)
                # MULTIPROCESS ASYNC
                exec_pool.apply_async(exec_fx_datasets, args=(fx_args,), callback=save_fx_datasets)

                id_ts += 1

            exec_pool.close()
            exec_pool.join()

            folder_name, file_name = os.path.split(file_pnt_ancillary)
            make_folder(folder_name)

            # DEBUG
            # fx_dset_common['sm'][0].plot()
            # fx_dset_common['swi_t06'][0].plot()
            # fx_dset_common['swi_t12'][0].plot()
            # fx_dset_common['swi_t32'][0].plot()
            # plt.show()

            write_obj(file_pnt_ancillary, fx_dset_common)

            time_end = time.time()
            time_elapsed = time_end - time_start

            logging.info(' ------> Analyze datasets for chunk ' + i_chunks_str + '/' + n_chunks_str +
                         ' ... DONE [ELAPSED: ' + str(np.floor(time_elapsed)) + ' seconds]')

        else:

            logging.info(' ------> Analyze datasets for chunk ' + i_chunks_str + '/' + n_chunks_str +
                         ' ... SKIPPED. Datasets already computed and saved')

    logging.info(' -----> Compute datasets ... DONE')
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to execute fx datasets
def exec_fx_datasets(func_obj_interface):

    if isinstance(func_obj_interface, list) and (func_obj_interface.__len__() == 1):
        func_obj_interface = func_obj_interface[0]

    id = func_obj_interface['id']
    tag = func_obj_interface['tag']
    geo_x = func_obj_interface['geo_x']
    geo_y = func_obj_interface['geo_y']
    file_data = func_obj_interface['file_name']
    flag_updating = func_obj_interface['flag_updating']
    filter_list = func_obj_interface['filter_list']
    filter_name = func_obj_interface['filter_name']

    if 'var_name' in list(func_obj_interface.keys()):
        var_name = func_obj_interface['var_name']
    else:
        var_name = 'soil_moisture'

    print(id)

    if flag_updating:
        if os.path.exists(file_data):
            os.remove(file_data)

    if not os.path.exists(file_data):

        # Get data
        # var_dframe_tmp = deepcopy(var_dframe_common)
        # point_ts_data = deepcopy(var_dframe_tmp.variable[id].values)
        # point_ts_data = deepcopy(var_dframe.variable[id].values)

        point_ts_data = var_data_common[id, :]

        # Create dataframe to perform drop and jd obj
        point_ts_dframe = pd.DataFrame(data=point_ts_data, index=var_time_common)
        point_ts_dframe = point_ts_dframe.dropna()

        # Get julian dates of time series
        time_index = point_ts_dframe.index
        jd_values = time_index.to_julian_date().values
        sm_values = point_ts_dframe.values.reshape([jd_values.shape[0]])

        # Calculate SWI T=10
        dict_info_extended = {}
        var_list = [var_name]
        for filter_step in filter_list:
            filter_var = filter_name.format(filter_step)
            swi_values = exp_filter(sm_values, jd_values, ctime=filter_step)
            dict_info_extended[filter_var] = (["time"], swi_values)

            var_list.append(filter_var)

        dict_info_attrs = {'id': id, 'longitude': geo_x, 'latitude': geo_y,
                           'tag': tag, 'var_list': var_list}
        dict_info_base = {var_name: (["time"], sm_values)}
        dict_info_point = {**dict_info_base, **dict_info_extended}

        file_dset_expected = xr.Dataset(coords={'time': (['time'], var_time_common)})
        file_dset_expected.coords['time'] = file_dset_expected.coords['time'].astype('datetime64[ns]')

        file_dset_tmp = xr.Dataset(data_vars=dict_info_point, coords={'time': (['time'], time_index)})

        file_dset_expected = file_dset_expected.combine_first(file_dset_tmp)
        file_dset_expected.attrs = dict_info_attrs

    else:

        file_dset_expected = None
        file_data = None

    return file_dset_expected, file_data

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to split file datasets
def split_file_datasets(var_obj,
                        var_chunks=25, format_chunks_element='{0:08d}',
                        file_ref='chunk.workspace', file_upd=True,
                        var_geo_x='longitude', var_geo_y='latitude', var_id='id'):

    if var_chunks is None:
        var_chunks = var_obj[var_id].shape[0]

    str_chunks_element = format_chunks_element.format(var_chunks)
    file_ref = file_ref.format(point_chunk_element=str_chunks_element)

    if file_upd:
        if os.path.exists(file_ref):
            os.remove(file_ref)

    logging.info(' -----> Split datasets in chunks ... ')

    if not os.path.exists(file_ref):

        id_array = var_obj[var_id]
        geo_x_array = var_obj[var_geo_x]
        geo_y_array = var_obj[var_geo_y]

        id_chunks = list(list_chunks(list(id_array), var_chunks))
        geo_x_chunks = list(list_chunks(list(geo_x_array), var_chunks))
        geo_y_chunks = list(list_chunks(list(geo_y_array), var_chunks))

        chunks_dict = {var_id: id_chunks, var_geo_x: geo_x_chunks, var_geo_y: geo_y_chunks}

        folder_name, file_name = os.path.split(file_ref)
        make_folder(folder_name)
        write_obj(file_ref, chunks_dict)

        logging.info(' -----> Split datasets in chunks ... DONE')

    else:

        chunks_dict = read_obj(file_ref)

        logging.info(' -----> Split datasets in chunks ... SKIPPED. Datasets already chunked')

    return chunks_dict
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read datasets file
def read_file_datasets(file_name, file_obj='data.workspace', file_upd=True,
                       var_name_in='Nmean', var_name_out=None,
                       var_geo_x_in='lon', var_geo_y_in='lat', var_id_in='ID', var_time_in='date',
                       var_geo_x_out='longitude', var_geo_y_out='latitude', var_id_out='id', var_time_out='time'):

    if file_upd:
        if os.path.exists(file_obj):
            os.remove(file_obj)

    logging.info(' -----> Read datasets from "' + os.path.split(file_name)[1] + '" file ... ')

    if not os.path.exists(file_obj):
        if os.path.exists(file_name):

            if var_name_out is None:
                var_name_out = var_name_in

            file_dataset = xr.open_dataset(file_name)

            values_variable = file_dataset[var_name_in].values
            values_geo_x = file_dataset[var_geo_x_in].values
            values_geo_y = file_dataset[var_geo_y_in].values
            values_id = file_dataset[var_id_in].values
            values_time = file_dataset[var_time_in].values

            variable_dict = {var_name_out: values_variable,
                             var_geo_x_out: values_geo_x, var_geo_y_out: values_geo_y,
                             var_id_out: values_id, var_time_out: pd.DatetimeIndex(values_time)}

            folder_name, file_name = os.path.split(file_obj)
            make_folder(folder_name)

            write_obj(file_obj, variable_dict)

            logging.info(' -----> Read datasets from "' + os.path.split(file_name)[1] + '" file ... DONE')

        else:
            logging.error(' ===> File "' + os.path.split(file_name)[1] + '" is not available')
            raise FileNotFoundError('File not found')

    else:
        variable_dict = read_obj(file_obj)

        logging.info(' -----> Read datasets from "' + os.path.split(file_name)[1] +
                     '" file ... SKIPPED. Datasets already read and imported')

    return variable_dict
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read datasets file
def read_file_geo(file_name,
                  var_geo_x_in='lon', var_geo_y_in='lat', var_id_in='ID', var_time_in='date',
                  var_geo_x_out='longitude', var_geo_y_out='latitude', var_id_out='id', var_time_out='time'):

    if os.path.exists(file_name):

        file_dataset = xr.open_dataset(file_name)

        values_geo_x = file_dataset[var_geo_x_in].values
        values_geo_y = file_dataset[var_geo_y_in].values
        values_id = file_dataset[var_id_in].values
        values_time = file_dataset[var_time_in].values

        variable_dict = {var_geo_x_out: values_geo_x, var_geo_y_out: values_geo_y,
                         var_id_out: values_id, var_time_out: pd.DatetimeIndex(values_time)}

    else:
        logging.error(' ===> File "' + os.path.split(file_name)[1] + '" is not available')
        raise FileNotFoundError('File not found')

    return variable_dict
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to chunk list in n element(s)
def list_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
# -------------------------------------------------------------------------------------


"""
for time_index_step in file_index_date:

    time_str_step_start = time_index_step.strftime('%Y-%m-%d 00:00')
    time_str_step_end = time_index_step.strftime('%Y-%m-%d 23:59')

    file_dataset_sliced = file_dataset.sel(date=slice(time_str_step_start, time_str_step_end))
"""