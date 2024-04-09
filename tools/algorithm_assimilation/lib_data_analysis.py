# -------------------------------------------------------------------------------------
# Libraries
import logging
import os
import numpy as np
import pandas as pd
import xarray as xr

from multiprocessing import Pool, cpu_count

from lib_data_io_generic import create_darray_3d, read_obj, write_obj
from lib_utils_system import make_folder

from pytesmo.scaling import get_scaling_function
import pytesmo.metrics as metrics

import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to define data list
def define_fx_datasets(model_obj, dset_obj):

    if not isinstance(model_obj, list):
        model_obj = [model_obj]
    if not isinstance(dset_obj, list):
        dset_obj = [dset_obj]

    if (model_obj.__len__() == 1) and (model_obj.__len__() < dset_obj.__len__()):
        model_obj = model_obj * dset_obj.__len__()
    elif (dset_obj.__len__() == 1) and (dset_obj.__len__() < model_obj.__len__()):
        dset_obj = dset_obj * model_obj.__len__()
    else:
        raise NotImplementedError('Case of model and dset list obj not implemented yet')

    return model_obj, dset_obj
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to compute average datasets
def compute_data_average(var_data, var_time, var_geo_x, var_geo_y,
                         dim_geo_x='longitude', dim_geo_y='latitude', dim_time='time',
                         coord_geo_x='longitude', coord_geo_y='latitude', coord_time='time'):

    var_time_length = var_time.__len__()
    var_data_scaled = var_data.values / var_time_length

    var_data_expected = np.zeros([var_data.shape[0], var_data.shape[1], var_time_length])

    var_data_expected[:, :, :] = var_data_scaled[:, :, np.newaxis]

    var_data_array = create_darray_3d(var_data_expected, var_time, var_geo_x, var_geo_y,
                                      geo_1d=True, var_name=None,
                                      coord_name_x=coord_geo_x, coord_name_y=coord_geo_y, coord_name_time=coord_time,
                                      dim_name_x=dim_geo_x, dim_name_y=dim_geo_y, dim_name_time=dim_time,
                                      dims_order=[dim_geo_y, dim_geo_x, dim_time])

    return var_data_array
# -------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Method to read model and datasets file(s)
def read_fx_file(file_model, file_dst):

    global model_obj
    global datasets_obj

    if os.path.exists(file_model):
        model_obj = read_obj(file_model)
    else:
        raise IOError('File' + file_model + ' not found')
    if os.path.exists(file_dst):
        datasets_obj = read_obj(file_dst)
    else:
        raise IOError('File' + file_dst + ' not found')
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Method to wrap analysis
def wrap_fx_analysis(file_model_obj, file_dst_obj,
                     file_analysis_data, file_analysis_stats,
                     terrain_data=None, map_data=None,
                     var_model_name='soil_moisture', var_dst_name='swi_t06',
                     var_chunks=1000, format_chunks_element='{0:08d}', format_chunk_group='{0:05d}',
                     cpu_n=50, cpu_max=None, flag_cleaning_file=True):

    global fx_common_data
    global fx_common_attrs
    global fx_common_time
    global model_obj
    global datasets_obj

    logging.info(' -----> Compute analysis ... ')

    if cpu_max is None:
        cpu_max = cpu_count() - 1
    if cpu_n > cpu_max:
        logging.warning(' ===> Maximum of recommended processes must be less then ' + str(cpu_max))
        logging.warning(' ===> Set number of process from ' + str(cpu_n) + ' to ' + str(cpu_max))
        #cpu_n = cpu_max

    file_model_list, file_dst_list = define_fx_datasets(file_model_obj, file_dst_obj)
    var_model_list, var_dst_list = define_fx_datasets(var_model_name, var_dst_name)

    map_sm_id = map_data['sm_id'].values  # [:50]
    map_terrain_id = map_data['terrain_id'].values  # [:50]
    map_terrain_geo_x = map_data['terrain_geo_x'].values  # [:50]
    map_terrain_geo_y = map_data['terrain_geo_y'].values  # [:50]

    map_sm_id_parts = list(list_chunks(list(map_sm_id), var_chunks))
    map_terrain_id_parts = list(list_chunks(list(map_terrain_id), var_chunks))
    map_terrain_geo_x_parts = list(list_chunks(list(map_terrain_geo_x), var_chunks))
    map_terrain_geo_y_parts = list(list_chunks(list(map_terrain_geo_y), var_chunks))

    for var_model_step, file_model_step, var_dst_step, file_dst_step, file_anl_data_step, file_anl_stats_step in zip(
            var_model_list, file_model_list, var_dst_list, file_dst_list,
            file_analysis_data, file_analysis_stats):

        logging.info(' ------> Analyze datasets variable "' + var_dst_step +
                     '" vs model variable "' + var_model_step + '" ... ')

        model_obj = None
        datasets_obj = None
        step_k = 0
        for chunk_i, (model_id_step, dst_id_step, geo_x_step, geo_y_step) in enumerate(
                zip(map_terrain_id_parts, map_sm_id_parts, map_terrain_geo_x_parts, map_terrain_geo_y_parts)):

            i_chunks_str = format_chunk_group.format(chunk_i)
            n_chunks = map_sm_id_parts.__len__()
            n_chunks_str = format_chunks_element.format(n_chunks)

            logging.info(' -------> Analyze datasets for chunk ' + i_chunks_str + '/' +
                         n_chunks_str + ' ... ')

            file_anl_data = file_anl_data_step.format(
                point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)
            file_anl_stats = file_anl_stats_step.format(
                point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)

            if flag_cleaning_file:
                if os.path.exists(file_anl_data):
                    os.remove(file_anl_data)
                if os.path.exists(file_anl_stats):
                    os.remove(file_anl_stats)

            logging.info(' --------> Compute point time-series ... ')

            if (not os.path.exists(file_anl_data)) or (not os.path.exists(file_anl_stats)):

                if (model_obj is None) or (datasets_obj is None):
                    logging.info(' --------> Read point datasets variable "' + var_dst_step + '" and model variable "' +
                                 var_model_step + '"... ')
                    read_fx_file(file_model_step, file_dst_step)
                    logging.info(' --------> Read point datasets variable "' + var_dst_step + '" and model variable "' +
                                 var_model_step + '"... DONE')

                fx_common_data = None
                fx_common_attrs = None
                fx_common_time = None
                exec_pool = Pool(cpu_n)
                for step_j, (model_id, dst_id, geo_x, geo_y) in enumerate(
                        zip(model_id_step, dst_id_step, geo_x_step, geo_y_step)):

                    fx_args = {'model_id': model_id, 'datasets_id': dst_id, 'step': step_k,
                               'model_var_name': var_model_step, 'datasets_var_name': var_dst_step,
                               'var_limit_min': 0.0, 'var_limit_max': 1.0, 'time_sample': '12H',
                               'longitude': geo_x, 'latitude': geo_y}

                    # SERIAL
                    #fx_result = exec_fx_analysis(fx_args)
                    #save_fx_analysis(fx_result)
                    # MULTIPROCESS ASYNC
                    exec_pool.apply_async(exec_fx_analysis, args=(fx_args,), callback=save_fx_analysis)

                    step_k += 1

                exec_pool.close()
                exec_pool.join()

                logging.info(' --------> Compute point time-series ... DONE')

                logging.info(' --------> Save point time-series ... ')

                fx_common_n = model_id_step.__len__()
                fx_common_range = list(range(0, fx_common_n, 1))

                fx_common_dset = xr.Dataset(coords={'time': (['time'], fx_common_time)})
                for fx_key, fx_data in fx_common_data.items():
                    fx_common_da = xr.DataArray(
                        fx_data, name=fx_key, dims=['n', 'time'],
                        coords={'n': (['n'], fx_common_range), 'time': (['time'], fx_common_time)})
                    fx_common_dset[fx_key] = fx_common_da
                fx_common_dset.attrs = fx_common_attrs

                folder_name, file_name = os.path.split(file_anl_data)
                make_folder(folder_name)
                write_obj(file_anl_data, fx_common_dset)

                folder_name, file_name = os.path.split(file_anl_stats)
                make_folder(folder_name)
                write_obj(file_anl_stats, fx_common_attrs)

                logging.info(' --------> Save point time-series ... DONE')

                logging.info(' -------> Analyze datasets for chunk ' + i_chunks_str + '/' +
                             n_chunks_str + ' ... DONE')

            else:
                logging.info(' -------> Analyze datasets for chunk ' + i_chunks_str + '/' +
                             n_chunks_str + ' ... SKIPPED. Analysis previously computed')

                if step_k == 0:
                    step_k = step_k + model_id_step.__len__() - 1
                else:
                    step_k = step_k + model_id_step.__len__()

        logging.info(' ------> Analyze datasets variable "' + var_dst_step +
                     '" vs model variable "' + var_model_step + '" ... DONE')

    logging.info(' -----> Compute analysis ... DONE')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to merge fx analysis
def merge_fx_analysis(file_analysis_data_grp, file_analysis_stats_grp,
                      file_analysis_data_clt, file_analysis_stats_clt,
                      terrain_data=None, map_data=None,
                      var_model_name='soil_moisture', var_dst_name='swi_t06',
                      var_datasets_name=None,
                      var_chunks=1000, format_chunks_element='{0:08d}', format_chunk_group='{0:05d}',
                      cpu_n=100, cpu_max=None, flag_cleaning_file=True):

    var_model_list, var_dst_list = define_fx_datasets(var_model_name, var_dst_name)

    map_sm_id = map_data['sm_id'].values
    map_terrain_id = map_data['terrain_id'].values
    map_terrain_geo_x = map_data['terrain_geo_x'].values
    map_terrain_geo_y = map_data['terrain_geo_y'].values

    map_sm_id_parts = list(list_chunks(list(map_sm_id), var_chunks))
    map_terrain_id_parts = list(list_chunks(list(map_terrain_id), var_chunks))
    map_terrain_geo_x_parts = list(list_chunks(list(map_terrain_geo_x), var_chunks))
    map_terrain_geo_y_parts = list(list_chunks(list(map_terrain_geo_y), var_chunks))

    for var_model_step, var_dst_step, \
        file_anl_data_grp_step, file_anl_stats_grp_step, \
        file_anl_data_clt_step, file_anl_stats_clt_step in zip(
            var_model_list, var_dst_list,
            file_analysis_data_grp, file_analysis_stats_grp, file_analysis_data_clt, file_analysis_stats_clt):

        logging.info(' ------> Merge datasets scaled variable "' + var_dst_step +
                     '" vs model variable "' + var_model_step + '" ... ')

        if (not os.path.exists(file_anl_data_clt_step)) or (not os.path.exists(file_anl_stats_clt_step)):

            dset_n = None
            dset_time = None
            var_data_collections = {}
            var_attrs_collections = {}
            for chunk_i, (model_id_step, dst_id_step, geo_x_step, geo_y_step) in enumerate(
                    zip(map_terrain_id_parts, map_sm_id_parts, map_terrain_geo_x_parts, map_terrain_geo_y_parts)):

                i_chunks_str = format_chunk_group.format(chunk_i)
                n_chunks = map_sm_id_parts.__len__()
                n_chunks_str = format_chunks_element.format(n_chunks)

                logging.info(' -------> Merge datasets for chunk ' + i_chunks_str + '/' +
                             n_chunks_str + ' ... ')

                file_anl_data_grp_def = file_anl_data_grp_step.format(
                    point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)
                file_anl_stats_grp_def = file_anl_stats_grp_step.format(
                    point_chunk_element=n_chunks_str, point_chunk_group=i_chunks_str)

                dset_data_step = read_obj(file_anl_data_grp_def)
                dset_stats_step = read_obj(file_anl_stats_grp_def)

                if var_datasets_name is None:
                    dset_vars_obj = list(dset_data_step.data_vars)
                else:
                    dset_var_tmp = list(dset_data_step.data_vars)
                    if var_datasets_name in dset_var_tmp:
                        var_idx = dset_var_tmp.index(var_datasets_name)
                        dset_vars_obj = dset_var_tmp[var_idx]
                    else:
                        logging.error(' ==> Datasets variable name ' + var_datasets_name + ' is not found')
                        raise RuntimeError('Variable name is not correct')

                if not isinstance(dset_vars_obj, list):
                    dset_vars_list = [dset_vars_obj]
                else:
                    dset_vars_list = dset_vars_obj

                for var_name in dset_vars_list:
                    var_data_step = dset_data_step[var_name].values

                    if dset_time is None:
                        dset_time = pd.DatetimeIndex(dset_data_step['time'].values)

                    if var_name not in list(var_data_collections.keys()):
                        var_data_collections[var_name] = {}
                        var_data_collections[var_name] = var_data_step
                        dset_n = var_data_step.shape[0]
                    else:
                        var_data_tmp = var_data_collections[var_name]
                        var_data_vstack = np.vstack((var_data_tmp, var_data_step))
                        var_data_collections[var_name] = var_data_vstack
                        dset_n = var_data_vstack.shape[0]

                for stats_name, stats_values in dset_stats_step.items():
                    if stats_name not in list(var_attrs_collections.keys()):
                        var_attrs_collections[stats_name] = stats_values
                    else:
                        stats_tmp = var_attrs_collections[stats_name]
                        stats_tmp.extend(stats_values)
                        var_attrs_collections[stats_name] = stats_tmp

                logging.info(' -------> Merge datasets for chunk ' + i_chunks_str + '/' +
                             n_chunks_str + ' ... DONE')

            logging.info(' -------> Save datasets data ... ')
            dset_n_range = list(range(0, dset_n, 1))

            var_dset = xr.Dataset(coords={'time': (['time'], dset_time)})
            for var_key, var_data in var_data_collections.items():
                var_da = xr.DataArray(
                    var_data, name=var_key, dims=['n', 'time'],
                    coords={'n': (['n'], dset_n_range), 'time': (['time'], dset_time)})
                var_dset[var_key] = var_da
            var_dset.attrs = var_attrs_collections

            folder_name, file_name = os.path.split(file_anl_data_clt_step)
            make_folder(folder_name)
            write_obj(file_anl_data_clt_step, var_dset)
            logging.info(' -------> Save datasets data ... DONE')

            logging.info(' -------> Save datasets statistics ... ')
            folder_name, file_name = os.path.split(file_anl_stats_clt_step)
            make_folder(folder_name)
            write_obj(file_anl_stats_clt_step, var_attrs_collections)
            logging.info(' -------> Save datasets statistics ... DONE')

            logging.info(' ------> Merge datasets scaled variable "' + var_dst_step +
                         '" vs model variable "' + var_model_step + '" ... DONE')

        else:

            logging.info(' ------> Merge datasets variable "' + var_dst_step +
                         '" vs model variable "' + var_model_step + '" ... SKIPPED. Analysis previously computed.')


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to compute ts values, score and information
def exec_fx_analysis(fx_args):

    i = fx_args['step']
    model_id = fx_args['model_id']
    datasets_id = fx_args['datasets_id']
    model_var_name = fx_args['model_var_name']
    datasets_var_name = fx_args['datasets_var_name']

    var_limit_min = fx_args['var_limit_min']
    var_limit_max = fx_args['var_limit_max']
    time_sample = fx_args['time_sample']
    geo_x = fx_args['longitude']
    geo_y = fx_args['latitude']

    no_data = -9999.0

    print(i)

    scaling_method_lr = get_scaling_function('linreg')
    scaling_method_ms = get_scaling_function('mean_std')

    datasets_tags = datasets_obj.tag
    datasets_idx = datasets_tags.index(datasets_id)

    model_idxs = model_obj.index.tolist()
    model_idx = model_idxs.index(model_id)

    datasets_ts_raw = datasets_obj[datasets_var_name][datasets_idx]
    model_ts_raw = model_obj[model_var_name][model_idx]

    time_ts_commom = pd.DatetimeIndex(model_ts_raw['time'].values)
    time_n_common = time_ts_commom.shape[0]

    datasets_values_tmp = datasets_ts_raw.values
    time_values_tmp = pd.DatetimeIndex(datasets_ts_raw['time'].values).round(time_sample)

    datasets_dframe_tmp = pd.DataFrame(data=datasets_values_tmp, index=time_values_tmp)
    datasets_dframe_tmp = datasets_dframe_tmp.loc[~datasets_dframe_tmp.index.duplicated(keep='first')]

    datasets_dframe_finite = datasets_dframe_tmp.dropna()
    time_ts_finite = datasets_dframe_finite.index
    datasets_values_finite = datasets_dframe_finite.values.flatten()

    datasets_ts_finite = xr.DataArray(
        datasets_values_finite, name=datasets_var_name, dims=['time'], coords={'time': (['time'], time_ts_finite)})
    model_ts_finite = model_ts_raw.sel(time=time_ts_finite) # .groupby('time.month').mean('time')

    datasets_values_check = datasets_ts_finite.dropna(dim='time').values
    model_values_check = model_ts_finite.dropna(dim='time').values

    if datasets_values_check.shape[0] == model_values_check.shape[0]:

        datasets_ts_scaled_ms = scaling_method_ms(datasets_ts_finite, model_ts_finite)
        datasets_ts_scaled_lr = scaling_method_lr(datasets_ts_finite, model_ts_finite)

        datasets_ts_scaled_ms[datasets_ts_scaled_ms > var_limit_max] = var_limit_max
        datasets_ts_scaled_ms[datasets_ts_scaled_ms < var_limit_min] = var_limit_min
        datasets_ts_scaled_lr[datasets_ts_scaled_lr > var_limit_max] = var_limit_max
        datasets_ts_scaled_lr[datasets_ts_scaled_lr < var_limit_min] = var_limit_min

        datasets_r_ms = metrics.pearsonr(datasets_ts_scaled_ms.values, model_ts_finite.values)[0]
        datasets_rho_ms = metrics.spearmanr(datasets_ts_scaled_ms.values, model_ts_finite.values)[0]

        datasets_r_lr = metrics.pearsonr(datasets_ts_scaled_lr.values, model_ts_finite.values)[0]
        datasets_rho_lr = metrics.spearmanr(datasets_ts_scaled_lr.values, model_ts_finite.values)[0]

        datasets_values_scaled_ms = datasets_ts_scaled_ms.values
        datasets_values_scaled_lr = datasets_ts_scaled_lr.values

    elif model_values_check.shape[0] == 0:
        logging.warning(' ===> Model point ' + str(model_id) + ' has all nan(s)')
        datasets_n = datasets_ts_finite.shape[0]
        datasets_values_scaled_ms = np.zeros([datasets_n])
        datasets_values_scaled_ms[:] = no_data
        datasets_values_scaled_lr = np.zeros([datasets_n])
        datasets_values_scaled_lr[:] = no_data
        datasets_r_ms = no_data
        datasets_rho_ms = no_data
        datasets_r_lr = no_data
        datasets_rho_lr = no_data
    elif datasets_values_check.shape[0] == 0:
        logging.warning(' ===> Datasets point ' + str(datasets_id) + ' has all nan(s)')
        model_n = model_ts_finite.shape[0]
        datasets_values_scaled_ms = np.zeros([model_n])
        datasets_values_scaled_ms[:] = no_data
        datasets_values_scaled_lr = np.zeros([model_n])
        datasets_values_scaled_lr[:] = no_data
        datasets_r_ms = no_data
        datasets_rho_ms = no_data
        datasets_r_lr = no_data
        datasets_rho_lr = no_data
    else:
        logging.warning(' ===> Datasets shape ' + str(datasets_values_check.shape[0]) +
                        ' and model shape ' + str(model_values_check.shape[0]) + ' are not equal')
        datasets_values_scaled_ms = np.zeros([time_n_common])
        datasets_values_scaled_ms[:] = no_data
        datasets_values_scaled_lr = np.zeros([time_n_common])
        datasets_values_scaled_lr[:] = no_data
        datasets_r_ms = no_data
        datasets_rho_ms = no_data
        datasets_r_lr = no_data
        datasets_rho_lr = no_data

    datasets_analysis_ts = {
        'datasets_scaled_ms': (["time"], datasets_values_scaled_ms),
        'datasets_scaled_lr': (["time"], datasets_values_scaled_lr)
    }

    datasets_analysis_attrs = {
        'datasets_pearson_ms': datasets_r_ms,
        'datasets_pearson_lr': datasets_r_lr,
        'datasets_spearman_ms': datasets_rho_ms,
        'datasets_spearman_lr': datasets_rho_lr,
        'datasets_var_name': datasets_var_name, 'model_var_name': model_var_name,
        'datasets_id': datasets_id, 'model_id': model_id,
        'longitude': geo_x, 'latitude': geo_y
    }

    # Debug plot
    '''
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    ax.plot(model_ts_finite.time, model_ts_finite.values, lw=1, color='#0000FF', label='hmc')
    ax.plot(datasets_ts_scaled_ms.time, datasets_ts_scaled_ms.values, lw=1, color='#BA3723', label='sentinel scaled ms')
    ax.plot(datasets_ts_scaled_lr.time, datasets_ts_scaled_lr.values, lw=1, color='#fcba03', label='sentinel scaled lr')

    ax.set_ylim(0, 1)
    ax.set_title('Soil Moisture Time Series [Person: ' + str(datasets_r_ms) + ']')
    ax.set_ylabel('soil moisture [-]')
    ax.grid(b=True)
    plt.legend()

    plt.show()
    '''

    dataset_analysis_expected = xr.Dataset(coords={'time': (['time'], time_ts_commom)})
    dataset_analysis_expected.coords['time'] = dataset_analysis_expected.coords['time'].astype('datetime64[ns]')

    dataset_analysis_tmp = xr.Dataset(data_vars=datasets_analysis_ts, coords={'time': (['time'], time_ts_finite)})

    dataset_analysis_expected = dataset_analysis_expected.combine_first(dataset_analysis_tmp)
    dataset_analysis_expected.attrs = datasets_analysis_attrs

    return dataset_analysis_expected
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to compute ts values, score and information
def save_fx_analysis(fx_result):
    global fx_common_data
    global fx_common_attrs
    global fx_common_time

    fx_attrs_tmp = fx_result.attrs

    fx_var_list = list(fx_result.data_vars)

    if fx_common_time is None:
        fx_common_time = pd.DatetimeIndex(fx_result['time'].values)

    if (fx_common_data is None) and (fx_common_attrs is None):

        fx_common_attrs = {}
        for key_attr, value_attr in fx_attrs_tmp.items():
            if not isinstance(value_attr, list):
                value_attr = [value_attr]
            fx_common_attrs[key_attr] = value_attr

        fx_common_data = {}
        for fx_var_step in fx_var_list:
            fx_step_data = fx_result[fx_var_step].values
            fx_common_data[fx_var_step] = fx_step_data
    else:

        for key_attr, value_attr in fx_attrs_tmp.items():
            fx_tmp_attr = fx_common_attrs[key_attr]
            fx_tmp_attr.append(value_attr)
            fx_common_attrs[key_attr] = fx_tmp_attr

        for fx_var_step in fx_var_list:
            fx_step_data = fx_result[fx_var_step].values
            fx_tmp_data = fx_common_data[fx_var_step]
            fx_vstack_data = np.vstack((fx_tmp_data, fx_step_data))
            fx_common_data[fx_var_step] = fx_vstack_data

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to chunk list in n element(s)
def list_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
# -------------------------------------------------------------------------------------
