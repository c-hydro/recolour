# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import os
import netCDF4
import pprint
import logging

from datetime import datetime

from hsaf_validation_metrics import HSAF_Metrics
from hsaf_validation_utils import get_dataset_names, open_log_file, close_log_file
from hsaf_validation_datasets import GLDAS_Dataset, ASCAT_Dataset_DR, ASCAT_Dataset_NRT, CCI_Dataset, RZSM_Dataset

from pytesmo.validation_framework.validation import Validation
from pytesmo.validation_framework.results_manager import netcdf_results_manager

from pygeogrids.grids import CellGrid
from numpy import zeros
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to setup validation process
def setup_process(val_param):

    # ----------------------------------------------------------------------------------------------------------------------
    # reference reader initialization
    ref_reader = None

    # Get flag to clean tmp datasets
    if 'flag_clean_tmp' in list(val_param.keys()):
        flag_clean_tmp = val_param['flag_clean_tmp']
    else:
        flag_clean_tmp = True
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # ASCAT data record reader obj
    if 'ascat_type_ts' in val_param:
        ref_name='ASCAT'
        ref_var=['sm']
        ref_kwargs = {'mask_frozen_prob': 10, 'mask_snow_prob': 10}
        if val_param['ascat_type_ts'] == 'data_record':    # H113 (data_record)
            ref_reader = ASCAT_Dataset_DR(dr_path=val_param['ascat_path_ts'],
                                          grid_path=val_param['ascat_path_grid'],
                                          static_layer_path=val_param['ascat_path_static'],
                                          flag_clean_tmp=flag_clean_tmp)
        elif val_param['ascat_type_ts'] == 'data_nrt':     # H16, H101, H102, H103 (swath time series)

            ref_reader = ASCAT_Dataset_NRT(dr_path=val_param['ascat_path_ts'],
                                           grid_path=val_param['ascat_path_grid'],
                                           static_layer_path=val_param['ascat_path_static'],
                                           flag_clean_tmp=flag_clean_tmp)
        else:
            raise IOError(' datasets type "' + val_param['ascat_type_ts'] + '" is not supported')

    if 'rzsm_type_ts' in val_param:
        if val_param['rzsm_type_ts'] == 'data_mod_dr':    # H140 (data_record)
            ref_name = 'RZSM'
            ref_var = ['var40']
            ref_kwargs = {}
            ref_reader = RZSM_Dataset(dr_path=val_param['rzsm_path_ts'],
                                      grid_path=val_param['rzsm_path_grid'],
                                      temp_path=val_param['rzsm_path_tmp'],
                                      flag_clean_tmp=flag_clean_tmp)
        elif val_param['rzsm_type_ts'] == 'data_mod_nrt': 
            ref_name = 'RZSM'
            ref_var = ['var40']
            ref_kwargs = {}
            ref_reader = RZSM_Dataset(dr_path=val_param['rzsm_path_ts'],
                                      grid_path=val_param['rzsm_path_grid'],
                                      temp_path=val_param['rzsm_path_tmp'],
                                      flag_clean_tmp=flag_clean_tmp)
        else:
            raise IOError(' datasets type "' + val_param['ascat_type_ts'] + '" is not supported')

    if ref_reader is None:
        raise RuntimeError('reference reader is NoneType. Check your settings file to correctly define the reader')

    # GLDAS time series reader obj
    gldas_reader = GLDAS_Dataset(val_param['gldas_path_ts'])
    # GLDAS time series reader obj
    cci_reader = CCI_Dataset(val_param['passive_path_ts'])
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Get lut max dist (if available in configuration file)
    if 'lut_max_dist' in val_param:
        lut_max_dist = int(val_param.pop('lut_max_dist'))
    else:
        lut_max_dist = 35000
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set dataset(s)
    datasets = {
        ref_name: {
            'class': ref_reader, 'columns': ref_var,
            'type': 'reference',
            'args': [],
            'kwargs': ref_kwargs
        },
        'GLDAS': {
            'class': gldas_reader,
            'columns': ['SoilMoi0_10cm_inst'],
            'type': 'other',
            'grids_compatible': False,
            'use_lut': True,
            'lut_max_dist': lut_max_dist,
        },
        'PASSIVE': {
            'class': cci_reader,
            'columns': ['sm'],
            'type': 'other',
            'kwargs': {'only_valid': True},
            'grids_compatible': False,
            'use_lut': True,
            'lut_max_dist': lut_max_dist,
        }
    }
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set anomaly reference
    anomaly = val_param.pop('anomaly', False)
    # Set seasonal metrics
    seasonal_metrics = val_param.pop('seasonal_metrics', False)

    # Set validation period
    if 'time_start' in val_param and 'time_end' in val_param:
        period = [datetime.strptime(val_param['time_start'], val_param['time_format']),
                  datetime.strptime(val_param['time_end'], val_param['time_format'])]
    else:
        period = None

    # Set validation temporal window
    if 'temporal_window' in val_param:
        window = float(val_param['temporal_window'])/24.
    else:
        window = 12. / 24.

    # Set validation results path
    if 'results_path' in val_param:
        results_path = val_param['results_path']
        if os.path.exists(results_path) is False:
            os.makedirs(results_path)
    else:
        print('No results path defined.')
        raise ValueError

    # Get datasets name
    dataset_names = get_dataset_names(ref_name, datasets)
    # Get metrics
    snr_metric = HSAF_Metrics(dataset_names=dataset_names,
                              seasonal_metrics=seasonal_metrics)

    # Define results path
    save_path = results_path

    # Define process
    process = Validation(datasets, ref_name,
                         {(3, 3): snr_metric.calc_metrics},
                         temporal_window=window, scaling=None,
                         period=period)
    # Define jobs
    jobs_temp = process.get_processing_jobs()

    return ref_name, process, save_path

    # ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to wrap process for starting process
def wrap_process(setting):

    # ----------------------------------------------------------------------------------------------------------------------
    # Parser of warp process function argument(s)
    group = setting[0]
    cell = setting[1]
    val_param = setting[2]
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Method to open file log
    log_handle = open_log_file(group, cell, val_param['logs_path'])
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Info start
    print(' => ARGS: ' + str(setting))
    print(' ==> GROUP: ' + str(group) + ' ... ')
    print(' ===> CELL: ' + str(cell) + ' ... ')
    log_handle.write(' => ARGS: ' + str(setting) + ' ... ' + '\n')
    log_handle.write(' ==> GROUP: ' + str(group) + ' ... ' + '\n')
    log_handle.write(' ===> CELL: ' + str(cell) + ' ... ' + '\n')
    log_handle.write(' ====> TIMESTART: ' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '\n')

    # Get process and save_path
    ref_name, process, save_path = setup_process(val_param)
    # Get reference reader
    ds_ref_reader = process.data_manager.datasets[ref_name]['class']
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Get gpi, lons and lats of datasets
    ref_gpi, ref_lons, ref_lats = ds_ref_reader.grid.grid_points_for_cell(cell)

    # Define ascat points dataset
    ref_points = zeros((len(ref_gpi), 3))
    ref_points[:, 0] = ref_gpi
    ref_points[:, 1] = ref_lons
    ref_points[:, 2] = ref_lats

    # Get jobs (by cell)
    jobs = []
    for ref_point in ref_points:
        jobs.append((int(ref_point[0]), ref_point[1], ref_point[2]))
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Iterate computation(s) over all jobs
    for idx, job in enumerate(jobs):

        # debug
        #print(' ::::: ACTIVE DEBUG :::::')
        #job = jobs[11]

        print(' =====> CELL: ' + str(cell) + ' JOB: ' + str(idx) + '/' + str(len(jobs)) +
              ' -- VALUE: ' + str(job) + ' ... ')
        log_handle.write(' =====> CELL: ' + str(cell) + ' JOB: ' + str(idx) + '/' + str(len(jobs)) +
                          ' -- VALUE: ' + str(job) + ' ... ' + '\n')

        results = start_process(job, process, save_path)

        print(' =====> CELL: ' + str(cell) + ' JOB: ' + str(idx) + '/' + str(len(jobs)) +
              ' -- VALUE: ' + str(job) + ' ... DONE!')
        log_handle.write(' =====> CELL: ' + str(cell) + ' JOB: ' + str(idx) + '/' + str(len(jobs)) +
                         ' -- VALUE: ' + str(job) + ' ... DONE!' + '\n')

    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Info start
    log_handle.write(' ====> TIMEEND: ' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '\n')
    log_handle.write(' ===> CELL: ' + str(cell) + ' ... DONE' + '\n')
    log_handle.write(' ==> GROUP: ' + str(group) + ' ... DONE!' + '\n' )

    print(' ===> CELL: ' + str(cell) + ' ... DONE!')
    print(' ==> GROUP: ' + str(group) + ' ... DONE!')

    return results
    # ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to run validation process
def start_process(job, process, save_path, write_to_cells=True):

    """
    Function to run validation process.

    Parameters
    ----------
    job : obj
        Job object.
    """

    if write_to_cells is False:
        return process.calc(*job)
    else:
        # write results in separate files per cell
        result = process.calc(*job)
        if save_path is None:
            return result

        if bool(result) is True:

            if process.data_manager.reference_grid is not None:
                if type(process.data_manager.reference_grid) is CellGrid:
                    grid = process.data_manager.reference_grid
                    try:
                        cell = grid.gpi2cell(job[0])
                    except:
                        cell = grid.gpi2cell(job[0][0])

                    gpis, lon, lat = grid.grid_points_for_cell(cell)
                    for key in result.keys():
                        fname = '{:04d}'.format(cell) + '.nc'
                        filename = os.path.join(save_path, fname)

                        if not os.path.exists(filename):
                            ncfile = netCDF4.Dataset(filename, 'w')

                            global_attr = {}
                            s = "%Y-%m-%d %H:%M:%S"
                            global_attr['date_created'] = datetime.now().strftime(s)
                            ncfile.setncatts(global_attr)

                            ncfile.createDimension('dim', None)
                        else:
                            ncfile = netCDF4.Dataset(filename, 'a')

                        index = len(ncfile.dimensions['dim'])
                        for field in result[key]:

                            if field in ncfile.variables.keys():
                                var = ncfile.variables[field]
                            else:
                                var_type = result[key][field].dtype
                                kwargs = {'fill_value': -99999}
                                # if dtype is a object the assumption is that the data is a
                                # string
                                if var_type == object:
                                    var_type = str
                                    kwargs = {}
                                var = ncfile.createVariable(field, var_type,
                                                            'dim', **kwargs)
                            var[index:] = result[key][field]

                        ncfile.close()

                        return result
                else:
                    return result
            else:
                return result
        else:
            return result

# ----------------------------------------------------------------------------------------------------------------------
