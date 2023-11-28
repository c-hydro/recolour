# NOTES on algorithm_plotting usage

The main script is app_validation_publisher.py.
The configuration file is app_validation_publisher_*.json.
The main script must be run with settings:
> python app_validation_publisher.py -settings_file app_validation_publisher_*.json

Notice that the configuration file determines exactly which maps are to be plotted and with which variables.

In general, to have all the plots which are needed, one has to plot xy_pr and xz_pr in separate runs.
Find and replace xy_pr or xz_pr occurrences in the configuration file and modify the names of the products in the name of the output plot accordingly.
Be careful not to modify the [destination][variables] field, but only the [figure][*_map][variables] fields for each plot as well as the [renderer][*_map][variable_in] fields.
Always make sure the destination folders match the naming convention in the validation folder.
Always make sure that the variables inside the validation netCDF files match what you are expecting.
Always remember that in the validation files, the naming convention for variables is in alphabetical order and correlation is order-independent.