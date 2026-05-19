"""
Library Features:

Name:          lib_utils_base
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20260518'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import re
from collections import OrderedDict
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to splits path(s)
def split_path(path):

    # check dynamic template
    match = re.search(r"(.*?)(\{.*\})", path)

    if match:
        root_part = match.group(1)
        dynamic_part = match.group(2)
    else:
        root_part = path
        dynamic_part = None

    return root_part, dynamic_part
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to create sf pattern
def create_sf_pattern(dynamic_part):
    sf_pattern = OrderedDict()

    if dynamic_part is not None:

        patterns = [
            ("%Y", ("year_folder", "{year}")),
            ("%m", ("month_folder", "{month}")),
            ("%d", ("day_folder", "{day}")),
            ("%H", ("hour_folder", "{hour}")),
            ("%M", ("minute_folder", "{minute}")),
        ]

        for token, (key, value) in patterns:
            if token in dynamic_part:
                sf_pattern[key] = value

    return sf_pattern
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to get dict values
def get_dict_values(
        dict_values,
        list_keys,
        default_value=None,
        mandatory=True):

    value = dict_values

    try:

        for key in list_keys:
            value = value[key]

    except (KeyError, TypeError):

        if mandatory:
            raise RuntimeError(
                f'Path "{list_keys}" not available'
            )

        value = default_value

    return value
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to get dict values
def map_dict_values(
        dict_values,
        dict_map,
        dict_default=None,
        dict_mandatory=None):

    if dict_default is None:
        dict_default = {}

    dict_out = {}

    for field_name, list_keys in dict_map.items():

        default_value = dict_default.get(field_name, None)

        # all mandatory if dict_mandatory is None
        if dict_mandatory is None:
            mandatory = True
        else:
            mandatory = dict_mandatory.get(field_name, False)

        field_value = get_dict_values(
            dict_values,
            list_keys,
            default_value=default_value,
            mandatory=mandatory
        )
        dict_out[field_name] = field_value

    return dict_out
# ----------------------------------------------------------------------------------------------------------------------
