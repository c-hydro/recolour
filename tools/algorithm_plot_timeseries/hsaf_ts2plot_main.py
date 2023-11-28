# ----------------------------------------------------------------------------------------------------------------------
# HSAF time series and plots
# https://pytesmo.readthedocs.io/en/latest/examples.html#the-pytesmo-validation-framework
# TEST over Italy domain --> [1394, 1395]
#
# General usage
# -configfile algorithm_configuration.json
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Libraries
import argparse
import progressbar

from textwrap import wrap
from copy import deepcopy
from numpy import where
from datetime import datetime
from os.path import join, exists
from pandas import read_pickle, Series

from hsaf_ts2plot_utils import parser_config_file, create_folder, open_file, compute_times
from hsaf_ts2plot_datasets import ASCAT_Dataset_DR, RISICO_Dataset

from pytesmo.validation_framework.data_manager import DataManager
from pytesmo.validation_framework.temporal_matchers import BasicTemporalMatching

import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import plt
# ----------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def GetArgs():

    # Parser algorithm arg(s)
    parser_obj = argparse.ArgumentParser()
    parser_obj.add_argument('-configfile', action="store", dest="config_file")
    parser_value = parser_obj.parse_args()

    if parser_value.config_file:
        config_file = parser_value.config_file
    else:
        config_file = 'algorithm_configuration.json'

    return config_file

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
if __name__ == "__main__":

    # -------------------------------------------------------------------------------------
    # Get script argument(s)
    config_file = GetArgs()
    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Parser configuration file to get validation parameter(s)
    val_param = parser_config_file(config_file)
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Get points workspace
    points_set = val_param['point']
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set reader(s)
    reader_ref = ASCAT_Dataset_DR(dr_path=val_param['ref_path_ts'],
                                    grid_path=val_param['ref_path_grid'],
                                    static_layer_path=val_param['ref_path_static'])

    reader_k1 = RISICO_Dataset(risico_data_folder=val_param['k1_path_ts'])
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set dataset(s)
    datasets = {
        'ASCAT': {
            'class': reader_ref,
            'columns': ['sm'],
            'type': 'reference',
            'args': [],
            'kwargs': {'mask_frozen_prob': 10, 'mask_snow_prob': 10}
        },
        'RISICO': {
            'class': reader_k1,
            'columns': ['UMB'],
            'type': 'other',
            'grids_compatible': False,
            'use_lut': True,
            'lut_max_dist': val_param['max_dist'],
        },
    }
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set date period
    period = [datetime.strptime(val_param['time_start'], val_param['time_format']),
              datetime.strptime(val_param['time_end'], val_param['time_format'])]
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Set data manager
    manager = DataManager(datasets, 'ASCAT', period)
    luts = manager.get_luts()
    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # Iterate over point(s)
    for point_id, point_data in points_set.items():

        # ---------------------------------------------------------------------------------------------------------------------
        # Get cell and gpis reference
        cell_ref = int(point_data['cell'])
        gpis_set = point_data['gpis']
        
        # Create folder(s) to save result(s)
        create_folder(val_param['ancillary_path'].format(cell_ref))
        create_folder(val_param['plot_path'].format(cell_ref))
        # ---------------------------------------------------------------------------------------------------------------------
        
        # ---------------------------------------------------------------------------------------------------------------------
        # Get gpis. lons, lats
        if gpis_set is None:
            gpis_ref = reader_ref.grid.grid_points_for_cell(cell_ref)[0]
        else:
            gpis_ref = gpis_set
        # ---------------------------------------------------------------------------------------------------------------------

        # ---------------------------------------------------------------------------------------------------------------------
        # Iterate over gpis
        pbar = progressbar.ProgressBar()
        for gpi_ref in pbar(gpis_ref):

            # ---------------------------------------------------------------------------------------------------------------------
            # Get score data
            if val_param['score_path']:
                if exists(val_param['score_path']):
                    filename = join(val_param['score_path'], '{:}.nc'.format(cell_ref))
                    filedata = open_file(filename, ['gpi', 'lon', 'lat', 'ALL_R'])

                    try:
                        index_ref = where(filedata['gpi'] == int(gpi_ref))[0][0]
                        score_ref = filedata['ALL_R'][index_ref]
                    except BaseException:
                        index_ref = None
                        score_ref = None

            # ---------------------------------------------------------------------------------------------------------------------

            # ---------------------------------------------------------------------------------------------------------------------
            # Get gpi(s), lon(s) and lat(s)
            gpi_ref = int(gpi_ref)
            gpi_k1 = luts['RISICO'][gpi_ref]

            lon_ref, lat_ref = reader_ref.grid.gpi2lonlat(gpi_ref)
            lon_k1, lat_k1 = reader_k1.grid.gpi2lonlat(gpi_k1)
            # ---------------------------------------------------------------------------------------------------------------------

            # ---------------------------------------------------------------------------------------------------------------------
            # Info starting plot
            print(' ==> GetData -- PeriodStart: ' +
                  val_param['time_start'] + ' PeriodEnd: ' + val_param['time_end'] + ' Cell: ' +
                  str(cell_ref) + ' [GPI Ref: ' + str(gpi_ref) + ' GPI k1: ' + str(gpi_k1) + ' ... ')
            # ---------------------------------------------------------------------------------------------------------------------

            # ---------------------------------------------------------------------------------------------------------------------
            # Get reference time-series
            file_ref = join(val_param['ancillary_path'].format(cell_ref), val_param['ref_ancillary_data'].format(
                gpi_ref, val_param['time_start'], val_param['time_end']))
            if not exists(file_ref):
                ts_ref = reader_ref.read_ts(gpi_ref)
                ts_ref.to_pickle(file_ref)
            else:
                ts_ref = read_pickle(file_ref)

            # Get k1 time-series
            file_k1 = join(val_param['ancillary_path'].format(cell_ref), val_param['k1_ancillary_data'].format(
                gpi_k1, val_param['time_start'], val_param['time_end']))
            if not exists(file_k1):
                ts_k1 = reader_k1.read_ts(gpi_k1)
                ts_k1.to_pickle(file_k1)
            else:
                ts_k1 = read_pickle(file_k1)

            # ----------------------------------------------------------------------------------------------------------------------

            # ---------------------------------------------------------------------------------------------------------------------
            # Info starting plot
            print(' ==> GetData -- PeriodStart: ' +
                  val_param['time_start'] + ' PeriodEnd: ' + val_param['time_end'] + ' Cell: ' +
                  str(cell_ref) + ' [GPI Ref: ' + str(gpi_ref) + ' GPI k1: ' + str(gpi_k1) + ' ... DONE!')
            # ---------------------------------------------------------------------------------------------------------------------
            
            # ----------------------------------------------------------------------------------------------------------------------
            # Generate dates group by month or year
            dates_start, dates_end = compute_times(val_param['time_start'], val_param['time_end'],
                                                   val_param['time_group'])
            # ----------------------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------------------
            # Iterate over dates start and end
            for date_start, date_end in zip(dates_start, dates_end):
                
                # ---------------------------------------------------------------------------------------------------------------------
                # Get time 
                time_start = date_start.strftime(val_param['time_format'])
                time_end = date_end.strftime(val_param['time_format'])

                # Plot starting info
                print(' ====> PlotData -- PeriodStart: ' + time_start + ' PeriodEnd: ' + time_end + ' ... ')
                # ---------------------------------------------------------------------------------------------------------------------
                
                # ---------------------------------------------------------------------------------------------------------------------
                # Get ref timeseries
                ts_ref['sm'] = ts_ref['sm'] / 100
                ts_ref_sel = ts_ref.loc[date_start:date_end]
                if ts_ref_sel.empty:
                    print(' ====> WARNING: Timeseries "ref" is empty!')
                    ts_ref_check = False
                else:
                    ts_ref_check = True
                    
                # Get k1 timeseries
                ts_k1_sel = ts_k1.loc[date_start:date_end]
                if ts_k1_sel.empty:
                    print(' ====> WARNING: Timeseries "k1" is empty!')
                    ts_k1_check = False
                else:
                    ts_k1_check = True
                # ---------------------------------------------------------------------------------------------------------------------
                
                # ---------------------------------------------------------------------------------------------------------------------
                # Check data availability
                if (ts_ref_check is True) and (ts_k1_check is True):
                
                    # ---------------------------------------------------------------------------------------------------------------------
                    # Temporal matching
                    if val_param['temporal_matching']:

                        window_match = float(val_param['temporal_matching']) / 24

                        df_dict = {}
                        df_dict['RISICO'] = ts_k1_sel
                        df_dict.update({'ASCAT': ts_ref_sel})

                        oTempMatch = BasicTemporalMatching(window=window_match)
                        results_matched = oTempMatch.combinatory_matcher(df_dict, 'ASCAT', 2)

                        data = results_matched['ASCAT', 'RISICO']

                        ts_ref_match = data['ASCAT']
                        ts_k1_match = data['RISICO']
                    else:
                        ts_ref_match = ts_ref_sel
                        ts_k1_match = ts_k1_sel
                        val_param['temporal_matching'] = None
                    # ---------------------------------------------------------------------------------------------------------------------

                    # ---------------------------------------------------------------------------------------------------------------------
                    # Data for scatter plot
                    values_ref = ts_ref_match['sm'].values
                    values_k1 = ts_k1_match['UMB'].values
                    series_data = Series(index=values_ref, data=values_k1)

                    scatter_data = series_data.to_frame()
                    scatter_data.reset_index(inplace=True)
                    scatter_data.columns = ['sm', 'UMB']
                    # ---------------------------------------------------------------------------------------------------------------------
                    
                    # ---------------------------------------------------------------------------------------------------------------------
                    # Flag active plotting
                    plot_active = val_param['plot_active']
                    # ---------------------------------------------------------------------------------------------------------------------
                
                    # ---------------------------------------------------------------------------------------------------------------------
                    # Plot results
                    if plot_active:
                    
                        # ---------------------------------------------------------------------------------------------------------------------
                        # Set and save plotting
                        file_plot = join(val_param['plot_path'],
                                         val_param['plot_ts'].format(
                                             cell_ref,
                                             gpi_ref, gpi_k1,
                                             time_start, time_end, val_param['temporal_matching']))

                        title_plot = \
                            'ASCAT [GPI: ' + str(gpi_ref) + ', Lon: ' + '{:3.3f}'.format(lon_ref) + ', Lat: ' + \
                            '{:3.3f}'.format(lat_ref) + '] vs RISICO [GPI: ' + str(gpi_k1) + ', Lon: ' + '{:3.3f}'.format(lon_k1) + \
                            ', Lat: ' + '{:3.3f}'.format(lat_k1) + '] - TemporalMatching: ' + \
                            str(val_param['temporal_matching']) + ' hours - Score PearsonR: ' + str(score_ref) + \
                            ' time_start: ' + time_start + ' date_end: ' + time_end

                        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
                        ts_ref_match['sm'].plot(ax=ax, alpha=0.4, marker='o', color='#FF0000', label='SSM')
                        ts_k1_match['UMB'].plot(ax=ax, lw=2, label='UMB')
                        plt.ylim(0, 1)
                        title = ax.set_title("\n".join(wrap(title_plot)))
                        plt.legend()

                        title.set_y(1.05)
                        fig.subplots_adjust(top=0.8)

                        fig.savefig(file_plot)
                        plt.close()
                        # ---------------------------------------------------------------------------------------------------------------------

                        # ---------------------------------------------------------------------------------------------------------------------
                        # Scatter plot
                        file_plot = join(val_param['plot_path'],
                                         val_param['plot_scatter'].format(
                                             cell_ref,
                                             gpi_ref, gpi_k1,
                                             time_start, time_end, val_param['temporal_matching']))

                        title_plot = \
                            'ASCAT [GPI: ' + str(gpi_ref) + ', Lon: ' + '{:3.3f}'.format(lon_ref) + ', Lat: ' + \
                            '{:3.3f}'.format(lat_ref) + '] vs RISICO [GPI: ' + str(gpi_k1) + ', Lon: ' + '{:3.3f}'.format(lon_k1) + \
                            ', Lat: ' + '{:3.3f}'.format(lat_k1) + '] - TemporalMatching: ' + \
                            str(val_param['temporal_matching']) + ' hours - Score PearsonR: ' + str(score_ref)  + \
                            ' time_start: ' + time_start + ' date_end: ' + time_end

                        data_line = range(1)

                        fig, ax = plt.subplots()

                        scatter_data.plot(ax=ax, kind='scatter', x='sm', y='UMB')

                        title = ax.set_title("\n".join(wrap(title_plot)))

                        plt.ylim(0, 1)
                        plt.xlim(0, 1)

                        fig.tight_layout()
                        title.set_y(1.05)
                        fig.subplots_adjust(top=0.8)

                        fig.savefig(file_plot)
                        plt.close()
                        # ---------------------------------------------------------------------------------------------------------------------

                        # ---------------------------------------------------------------------------------------------------------------------
                        # Plot ending info
                        print(' ====> PlotData -- PeriodStart: ' +
                              time_start + ' PeriodEnd: ' + time_end + ' ... OK')
                        # ---------------------------------------------------------------------------------------------------------------------

                    else:
                        #---------------------------------------------------------------------------------------------------------------------
                        # Plot ending info
                        print(' ====> PlotData -- PeriodStart: ' +
                              time_start + ' PeriodEnd: ' + time_end + ' ... SKIPPED! Plotting is not active!')
                        # ---------------------------------------------------------------------------------------------------------------------  
                else:
                    # ---------------------------------------------------------------------------------------------------------------------
                    # Plot ending info
                    print(' ====> PlotData -- PeriodStart: ' +
                          time_start + ' PeriodEnd: ' + time_end + ' ... SKIPPED! Some data are not available!')
                    # ---------------------------------------------------------------------------------------------------------------------

                # ---------------------------------------------------------------------------------------------------------------------

            # ----------------------------------------------------------------------------------------------------------------------

        # ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
