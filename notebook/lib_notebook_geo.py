"""
Library Features:

Name:          lib_notebook_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220320'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read point data
def read_point_data(file_name,
                    file_delimiter=',', file_header=0, file_columns_remap=None,
                    file_units_default='%', file_description_default='NA',
                    file_location_default='NA', file_tag_default='NA'):

    # message read time-series registry start
    print(' ---> Read time-series registry ... ')

    if file_columns_remap is None:
        file_columns_remap = {"altitude": "altitude", "amm_level_1": "amm_level_1", "amm_level_2": "amm_level_2",
                              "code": "code", "longitude": "longitude", "latitude": "latitude",
                              "description": "description", "name": "name", "valid": "valid",
                              "units": "units", "location": "locations", "tag": "tag"}
    remap_keys, remap_values = list(file_columns_remap.keys()), list(file_columns_remap.values())

    df_point = pd.read_csv(file_name, delimiter=file_delimiter, header=file_header)
    df_point.columns = df_point.columns.str.replace(' ', '')
    columns_keys = list(df_point.columns)

    df_remap = {}
    for remap_key, remap_value in zip(remap_keys, remap_values):
        if remap_value in columns_keys:
            df_remap[remap_value] = remap_key
        else:
            print(' ===> Column "' + remap_value + '" is not available in the source registry')

    df_point_remap = df_point.rename(columns=df_remap)

    select_keys, select_values = list(df_remap.keys()), list(df_remap.values())
    if "units" not in columns_keys:
        df_point_remap['units'] = file_units_default
        select_values.append('units')
    if "description" not in columns_keys:
        df_point_remap['description'] = file_description_default
        select_values.append('description')
    if "location" not in columns_keys:
        df_point_remap['location'] = file_location_default
        select_values.append('location')
    if "tag" not in columns_keys:
        df_point_remap['tag'] = file_tag_default
        select_values.append('tag')

    df_point_remap['name'] = df_point_remap['name'].str.strip()
    df_point_remap['description'] = df_point_remap['description'].str.strip()
    df_point_remap['location'] = df_point_remap['location'].str.strip()
    df_point_remap['tag'] = df_point_remap['tag'].str.strip()

    select_values = list(set(select_values))

    df_point_select = df_point_remap[select_values]

    if 'valid' in list(df_point_remap.columns):
        df_point_select = df_point_select.loc[df_point_remap['valid'] == 1]

    if 'tag' in list(df_point_remap.columns):
        if all(df_point_remap['tag'] == 'NA'):

            point_name = df_point_select['name']
            list_tag = []
            for string_name in point_name.values:
                string_tag = sanitize_string(string_name)
                list_tag.append(string_tag)

            df_point_select['tag'] = list_tag

    # message read time-series registry end
    print(' ---> Read time-series registry ... DONE')

    return df_point_select
# ----------------------------------------------------------------------------------------------------------------------
