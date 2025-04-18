"""
Library Features:

Name:          lib_data_io_ascii
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231010'
Version:       '1.0.0'
"""
import glob

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import pandas as pd

from copy import deepcopy

from lib_data_io_generic import combine_data_point_by_time

from lib_utils_system import fill_tags2string
from lib_utils_obj import map_vars_dframe, sanitize_string, fill_tags_time
from lib_utils_time import replace_time_part
from lib_info_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to wrap datasets in ascii format
def wrap_datasets_ascii(file_path_template,
                        file_fields, registry_fields,
                        time_reference, time_start, time_end, time_range=None,
                        time_rounding='H', time_frequency='Y', time_format='%Y%m%d%H%M',
                        template_time_tags=None, template_datasets_tags=None,
                        file_sep=' ', file_decimal='.',
                        ascending_index=False, sort_index=True):

    # check time frequency
    if time_frequency == 'Y':

        time_now_stamp = pd.Timestamp.now()
        time_end_stamp = pd.Timestamp(time_end)
        if time_now_stamp > time_end_stamp:
            log_stream.warning(
                ' ===> "time_now" "' + str(time_now_stamp) + '" is greater than "time_end" "' +
                str(time_end_stamp) + '". Time end is set to time now. '
                                      'Change time end in the settings file to avoid this message')
            time_end_stamp = deepcopy(time_now_stamp)
            time_end_stamp = time_end_stamp.replace(month=12, day=31, hour=23, minute=0, second=0, microsecond=0)

            time_end = time_end_stamp.strftime(time_format)

        time_file_start = pd.date_range(time_start, time_end, freq='AS')
        time_file_start = replace_time_part(time_file_start, time_rounding=time_rounding, time_value=0)
        time_file_end = pd.date_range(time_start, time_end, freq='Y')
        time_file_end = replace_time_part(time_file_end, time_rounding=time_rounding, time_value=23)

    elif time_frequency == 'h' or time_frequency == 'H':

        time_step_stamp = pd.Timestamp(time_end)
        time_step_stamp = time_step_stamp.round(time_frequency.lower())
        time_file_start, time_file_end = ['*'], pd.DatetimeIndex([time_step_stamp])

    else:
        log_stream.error(' ===> Time frequency "' + time_frequency + '" is not expected')
        raise NotImplementedError('Case not implemented yet')

    # iterate over registry fields
    section_data_collections, time_start_collections, time_end_collections = {}, None, None
    for registry_row in registry_fields.iterrows():
        # get point information
        registry_code, registry_tag = registry_row[1]['code'], registry_row[1]['tag']

        # info point start
        log_stream.info(' ------> Point (1) Code "' + registry_code + '" (2) Tag "' + registry_tag + '" ... ')

        # iterate over times
        fields_data_collections = None
        for time_step_start, time_step_end in zip(time_file_start, time_file_end):

            # define time period tag
            if time_step_start is not None and time_step_end is not None:
                time_period_tag = 'from "' +str(time_step_start) + '" to "' + str(time_step_end) + '"'
            elif time_step_start is None and time_step_end is not None:
                time_period_tag = str(time_step_end)
            else:
                log_stream.error(' ===> Time period case is not expected')
                raise NotImplementedError('Case not implemented yet')

            # info time period start
            log_stream.info(' -------> Get datasets for time reference ' + time_period_tag + ' ... ')

            # fill time tags
            template_time_values = fill_tags_time(
                template_time_tags,
                time_reference=time_reference, time_start=time_step_start, time_end=time_step_end)
            template_datasets_values = {'point_name': registry_code}

            # create template tags and values
            template_generic_tags = {**template_time_tags, **template_datasets_tags}
            template_generic_values = {**template_time_values, **template_datasets_values}

            # define file path
            file_path_defined = fill_tags2string(
                file_path_template, tags_format=template_generic_tags, tags_filling=template_generic_values)[0]

            # search for undefined time start or time end file(s)
            if '*' in file_path_defined:
                file_path_list = glob.glob(file_path_defined)
                if not file_path_list:
                    file_path_selected = None
                else:
                    file_path_selected = file_path_list[0]
            else:
                file_path_selected = file_path_defined

            # check file availability
            if (file_path_selected is not None) and (os.path.exists(file_path_selected)):

                # get file fields
                fields_data_raw = pd.read_table(
                    file_path_selected, sep=file_sep, decimal=file_decimal)
                fields_data_raw.columns = fields_data_raw.columns.str.strip()
                # map file fields
                fields_data_map = map_vars_dframe(fields_data_raw, file_fields)

                # parse time field
                if 'time' in list(fields_data_map.columns):
                    time_data_raw = fields_data_map['time'].values
                    time_data_index = pd.to_datetime(time_data_raw, format=time_format)
                    fields_data_map['time'] = time_data_index
                    fields_data_map.set_index('time', inplace=True)
                    fields_data_map.index.name = 'time'
                else:
                    log_stream.error(' ===> Field "time" is not defined in the file, but it is expected')
                    raise RuntimeError('Check your settings file and set the "time" field')

                # convert data to float (if integer are found)
                fields_data_map = fields_data_map.astype(float)

                # manage data collections
                if fields_data_collections is None:
                    fields_data_collections = deepcopy(fields_data_map)
                else:
                    fields_data_collections = pd.concat([fields_data_collections, fields_data_map])

                # info time period end (done)
                log_stream.info(' -------> Get datasets for time reference  ' + time_period_tag + ' ... DONE')

            else:
                # info time period end (skipped - file is not available)
                log_stream.info(' -------> Get datasets for time reference  ' + time_period_tag +
                                ' ... SKIPPED. File "' + file_path_defined + '" does not exists.')

        # info organize dataset start
        log_stream.info(' -------> Organize datasets ... ')

        # check if data is available
        if fields_data_collections is not None:

            # check time range
            if time_range is not None:
                # Create an empty DataFrame with the time index
                fields_data_expected = pd.DataFrame(index=time_range)
                fields_data_expected.index.name = "time"
                # join data to the expected time index
                fields_data_expected = fields_data_expected.join(fields_data_collections)
            else:
                # get datasets
                fields_data_expected = deepcopy(fields_data_collections)

            # sort index
            if sort_index:
                if ascending_index:
                    fields_data_expected = fields_data_expected.sort_index(ascending=True)
                    section_time_start, section_time_end = fields_data_expected.index[0], fields_data_expected.index[-1]
                else:
                    fields_data_expected = fields_data_expected.sort_index(ascending=False)
                    section_time_start, section_time_end = fields_data_expected.index[-1], fields_data_expected.index[0]
            else:
                section_time_start, section_time_end = fields_data_expected.index[0], fields_data_expected.index[-1]

            # store time start and end
            if time_start_collections is None: time_start_collections = []
            time_start_collections.append(section_time_start)
            if time_end_collections is None: time_end_collections = []
            time_end_collections.append(section_time_end)

            # info organize dataset end (done)
            log_stream.info(' -------> Organize datasets ... DONE')

        else:
            # info organize dataset end (skipped - dataset are not available)
            fields_data_expected = None
            log_stream.info(' -------> Organize datasets ... SKIPPED. No data available for the selected time period')

        # store section data to common workspace
        section_data_collections[registry_tag] = fields_data_expected

        # info point end
        log_stream.info(' ------> Point (1) Code "' + registry_code + '" (2) Tag "' + registry_tag + '" ... DONE')

    # check collections
    section_elements = list(section_data_collections.values())
    if all(element is None for element in section_elements):
        section_data_collections = None
        section_time_start, section_time_end, section_time_file = None, None, None
    else:
        section_time_start, section_time_end = list(set(time_start_collections)), list(set(time_end_collections))
        if len(section_time_start) > 1:
            log_stream.error(' ===> Time start is not the same for all datasets')
            raise RuntimeError('Check your settings file and set the "time" field')
        elif len(section_time_start) == 1:
            section_time_start = section_time_start[0]
        else:
            log_stream.error(' ===> Time start is not defined')
            raise RuntimeError('Check your settings file and set the "time" field')
        if len(section_time_end) > 1:
            log_stream.error(' ===> Time end is not the same for all datasets')
            raise RuntimeError('Check your settings file and set the "time" field')
        elif len(section_time_end) == 1:
            section_time_end = section_time_end[0]
        else:
            log_stream.error(' ===> Time end is not defined')
            raise RuntimeError('Check your settings file and set the "time" field')

        if time_range is not None:
            if time_reference not in time_range:
                section_time_file = section_time_end
            else:
                section_time_file = time_reference
        else:
            section_time_file = time_reference

    return section_data_collections, section_time_start, section_time_end, section_time_file
# ----------------------------------------------------------------------------------------------------------------------
