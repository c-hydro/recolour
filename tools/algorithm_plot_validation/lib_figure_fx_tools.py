"""
Library Features:

Name:          lib_figure_fx_tools
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20240411'
Version:       '1.1.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to configure plot kwargs
def configure_plot_kwargs(figure_settings, plot_keys=None):

    if plot_keys is None:
        plot_keys = ['cb_label', 'show_cb', 'title', 'vmin', 'vmax', 'cmap',
                     'title_fontsize', 'title_fontweight', 'cb_fontsize', 'cb_fontweight']

    plot_settings = {}
    for figure_key, figure_value in figure_settings.items():
        if figure_key in plot_keys:
            plot_settings[figure_key] = figure_value
    return plot_settings
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get field
def get_field(field_obj, field_name, field_mandatory=True, field_default=None):
    if field_name in list(field_obj.keys()):
        return field_obj[field_name]
    else:
        if field_mandatory:
            logging.error(' ===> Field name "' + field_name + '" is not available in the field obj')
            raise KeyError('Check the figure settings and add the field needed by the algorithm')
        else:
            logging.warning(' ===> Field name "' + field_name + '" is not available in the field obj')
            return field_default
# ----------------------------------------------------------------------------------------------------------------------
