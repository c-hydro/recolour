"""
Class Features

Name:          drv_map_dynamic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20241112'
Version:       '1.1.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os

from lib_utils_generic import make_folder, reset_folder, pad_list
from lib_utils_datasets import organize_datasets_cell, organize_datasets_grid
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_nc import read_file_collection, write_file_collection, organize_file_nc, write_file_nc
from lib_data_io_tiff import organize_file_tiff, write_file_tiff

from lib_data_analysis import get_data, filter_data, interpolate_data, resample_data_args, resample_data_fx, mask_data

# debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver map
class DrvMap:

    # method to initialize class
    def __init__(self, alg_obj_static, alg_settings,
                 tag_section_flags='flags', tag_section_info='info',
                 tag_section_grid='grid',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_time='time', tag_section_log='log'):

        self.alg_obj_static = alg_obj_static

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_info = alg_settings[tag_section_info]
        self.alg_parameters = alg_settings[tag_section_params]

        self.alg_time = None

        self.alg_datasets_src = alg_settings[tag_section_datasets]['dynamic']['source']
        self.alg_datasets_anc_points = alg_settings[tag_section_datasets]['dynamic']['ancillary']['points']
        self.alg_datasets_anc_grid = alg_settings[tag_section_datasets]['dynamic']['ancillary']['grid']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['dynamic']['destination']

        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars_in = 'variables_in'
        self.tag_file_vars_out = 'variables_out'
        self.tag_format = 'format'

        self.alg_datasets_name = self.alg_info['datasets']
        self.alg_datasets_product = self.alg_info['product']

        self.reset_datasets_anc_points = self.alg_flags['reset_ancillary_datasets_points']
        self.reset_datasets_anc_grid = self.alg_flags['reset_ancillary_datasets_grid']
        self.reset_datasets_dst = self.alg_flags['reset_destination_datasets']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.variables_in_src = self.alg_datasets_src[self.tag_file_vars_in]
        self.variables_out_src = self.alg_datasets_src[self.tag_file_vars_out]
        if self.tag_format in list(self.alg_datasets_src.keys()):
            self.format_src = self.alg_datasets_src[self.tag_format]
        else:
            self.format_src = 'netcdf'
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)

        self.folder_name_anc_points = self.alg_datasets_anc_points[self.tag_folder_name]
        self.file_name_anc_points = self.alg_datasets_anc_points[self.tag_file_name]
        self.file_path_anc_points = os.path.join(self.folder_name_anc_points, self.file_name_anc_points)

        self.folder_name_anc_grid = self.alg_datasets_anc_grid[self.tag_folder_name]
        self.file_name_anc_grid = self.alg_datasets_anc_grid[self.tag_file_name]
        self.file_path_anc_grid = os.path.join(self.folder_name_anc_grid, self.file_name_anc_grid)

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        if self.tag_file_vars_in in list(self.alg_datasets_dst.keys()) and \
                self.tag_file_vars_out in list(self.alg_datasets_dst.keys()):
            self.variables_in_dst = self.alg_datasets_dst[self.tag_file_vars_in]
            self.variables_out_dst = self.alg_datasets_dst[self.tag_file_vars_out]
        else:
            self.variables_in_dst, self.variables_out_dst = None, None
        if self.tag_format in list(self.alg_datasets_dst.keys()):
            self.format_dst = self.alg_datasets_dst[self.tag_format]
        else:
            self.format_dst = 'netcdf'
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

        self.grid_cells_product = self.alg_obj_static['grid_cells_product']
        self.grid_gpis_product = self.alg_obj_static['grid_gpis_product']
        self.grid_reference_product = self.alg_obj_static['grid_reference_product']
        self.grid_reference_domain = self.alg_obj_static['grid_reference_domain']

        self.resampling_max_distance = self.alg_parameters['resampling_max_distance']
        self.resampling_grid_resolution = self.alg_parameters['resampling_grid_resolution']
        self.resampling_min_neighbours = self.alg_parameters['resampling_min_neighbours']
        self.resampling_neighbours = self.alg_parameters['resampling_neighbours']
        self.resampling_variables = self.alg_parameters['resampling_variable']

    # method to check datasets (if clean or not)
    @staticmethod
    def _clean_ancillary_folders(dset_obj,
                                 dset_key_root=None, dset_key_sub=None, dset_clean=True):

        if (dset_key_root is not None) and (dset_key_sub is not None):
            for dset_key, dset_fields in dset_obj.items():
                dset_path = dset_fields[dset_key_root][dset_key_sub]
                reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is not None):
            dset_path = dset_obj[dset_key_sub]
            reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is not None) and (dset_key_sub is None):
            dset_path = dset_obj[dset_key_root]
            reset_folder(dset_path, folder_reset=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is None):
            pass
        else:
            pass

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize datasets ... ')

        # get grid reference and cells
        grid_cells_product = self.grid_cells_product
        grid_reference_product = self.grid_reference_product
        # get file variables
        variables_in_src, variables_out_src = self.variables_in_src, self.variables_out_src
        # get file path
        file_path_src = self.file_path_src
        file_path_anc_points, file_path_anc_grid = self.file_path_anc_points, self.file_path_anc_grid
        file_path_dst = self.file_path_dst

        # clean ancillary and destination datasets (if ancillary flag is activated
        if self.reset_datasets_anc_points or self.reset_datasets_anc_grid:
            if os.path.exists(file_path_anc_points):
                os.remove(file_path_anc_points)
            if os.path.exists(file_path_anc_grid):
                os.remove(file_path_anc_grid)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)
        # clean ancillary and destination datasets if are not available together
        if (not os.path.exists(file_path_anc_points)) or (not os.path.exists(file_path_dst)):
            if os.path.exists(file_path_anc_points):
                os.remove(file_path_anc_points)
            if os.path.exists(file_path_anc_grid):
                os.remove(file_path_anc_grid)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        # check file point ancillary availability
        if not os.path.exists(file_path_anc_points):

            # check source format
            if self.format_src == 'netcdf':

                # join cell object(s)
                cell_obj_dict = organize_datasets_cell(
                    cell_list=grid_cells_product, cell_digits=4,
                    list_variable_in=variables_in_src, list_variable_out=variables_out_src,
                    folder_name_datasets=self.folder_name_src, file_name_datasets=self.file_name_src,
                    )
                # filter cell object
                cell_obj_dict = organize_datasets_grid(
                    cell_obj_dict, data_type=self.alg_datasets_name,
                    grid_dframe=grid_reference_product, file_name_grid=None)

                # save cell obj in pickle format
                # folder_name_anc_points, file_name_anc_points = os.path.split(file_path_anc_points)
                # make_folder(folder_name_anc_points)
                # write_file_obj(file_path_anc_points, cell_obj_dict)

                # save cell obj in netcdf format
                folder_name_anc_points, file_name_anc_points = os.path.split(file_path_anc_points)
                make_folder(folder_name_anc_points)
                write_file_collection(file_path_anc_points, cell_obj_dict)

                # info end method
                logging.info(' ---> Organize datasets ... DONE')

            else:
                # info end method (failed)
                logging.info(' ---> Organize datasets ... FAILED. Source format not allowed')
                raise NotImplementedError('Source format not implemented yet')

        else:
            # info end method (skipped)
            logging.info(' ---> Organize datasets ... SKIPPED. Data points previously saved')

    # method to analyze data
    def analyze_data(self):

        # info start method
        logging.info(' ---> Analyze data ... ')

        # get file path
        file_path_anc_points, file_path_anc_grid = self.file_path_anc_points, self.file_path_anc_grid

        # get resampling variables info
        res_var_name = self.resampling_variables['name']
        res_var_min_value = self.resampling_variables['min_value']
        res_var_max_value = self.resampling_variables['max_value']
        # pad list (using names as list lenght reference
        res_var_min_value, res_var_max_value = pad_list(res_var_name, [res_var_min_value, res_var_max_value])

        # get grid reference and cells
        grid_reference_product = self.grid_reference_product
        grid_reference_domain = self.grid_reference_domain
        # get variables
        variables_out_src = self.variables_out_src

        # check file ancillary points availability
        if os.path.exists(file_path_anc_points):

            # check file ancillary grid availability
            if not os.path.exists(file_path_anc_grid):

                # method to get data cell
                dframe_product_in = read_file_collection(
                    file_path_anc_points, file_obj_grid=grid_reference_product, variable_list=variables_out_src)

                # iterate over variable(s)
                variable_collections = None
                for var_name, var_min_value, var_max_value in zip(res_var_name, res_var_min_value, res_var_max_value):

                    # debug
                    # var_name = 'xy_pr'

                    # info start variable
                    logging.info(' ----> Variable "' + var_name + '" ... ')

                    # get dframe variable(s)
                    var_dframe = list(dframe_product_in.columns)

                    # check if variable is available
                    if var_name in var_dframe:

                        # method to get data
                        dframe_variable_in = get_data(dframe_product_in, variable_name=var_name)
                        # method to filter data
                        dframe_variable_filtered = filter_data(
                            dframe_variable_in,
                            variable_name=var_name, min_value=var_min_value, max_value=var_max_value)

                        # method to create resampling argument(s)
                        resample_kwargs = resample_data_args(
                            grid_reference_domain,
                            self.resampling_max_distance, self.resampling_grid_resolution)

                        # method to apply resampling data 1D to data 2D
                        da_variable_resampled = resample_data_fx(
                            dframe_variable_filtered, **resample_kwargs,
                            var_name_data=var_name, var_name_geo_x='lon', var_name_geo_y='lat',
                            coord_name_x='longitude', coord_name_y='latitude', dim_name_x='longitude', dim_name_y='latitude'
                        )

                        # method to interpolate data
                        da_variable_interpolated = interpolate_data(da_variable_resampled, grid_reference_domain,
                                                                    method='nearest')
                        # method to mask data
                        da_variable_masked = mask_data(
                            da_variable_interpolated, grid_reference_domain,
                            var_name_data=var_name, mask_value_no_land=0)

                        # save variable in a common collection
                        if variable_collections is None:
                            variable_collections = {}
                        variable_collections[var_name] = da_variable_masked

                        ''' debug
                        import matplotlib.pylab as plt
                        plt.figure(1)
                        plt.imshow(grid_reference_domain.values)
                        plt.colorbar()
                        plt.figure(2)
                        plt.imshow(da_variable_masked.values)
                        plt.colorbar()
                        plt.clim(0, 1)
                        plt.show()
                        '''

                        # info end variable
                        logging.info(' ----> Variable "' + var_name + '" ... DONE')

                    else:
                        # info end variable
                        logging.warning(' ===> Variable "' + var_name + '" not available in the dataset(s)')
                        logging.info(' ----> Variable "' + var_name + '" ... SKIPPED.')

                # check if variable collection is available
                if variable_collections is not None:

                    # save variable collection in pickle format
                    folder_name_anc_grid, file_name_anc_grid = os.path.split(file_path_anc_grid)
                    make_folder(folder_name_anc_grid)
                    write_file_obj(file_path_anc_grid, variable_collections)

                    # info end method
                    logging.info(' ---> Analyze data ... DONE')
                else:
                    # info end method
                    logging.info(' ---> Analyze data ... FAILED')
                    logging.error(' ===> Variable collection is defined by NoneType. Check your datasets path(s)')
                    raise RuntimeError('Variable collection is empty. Check the expected and the datasets variable(s)')

            else:
                # info end method
                logging.info(' ---> Analyze data ... SKIPPED. Data grid previously saved')

        else:
            # info end method
            logging.info(' ---> Analyze data... FAILED. Data points not available')

    # method to save data
    def dump_data(self):

        # info start method
        logging.info(' ---> Dump datasets ... ')

        # get file path
        file_path_anc_grid, file_path_dst = self.file_path_anc_grid, self.file_path_dst

        # check file ancillary grid availability
        if os.path.exists(file_path_anc_grid):
            if not os.path.exists(file_path_dst):

                # method to get data obj
                variable_collection = read_file_obj(file_path_anc_grid)

                # check destination format
                if self.format_dst == 'netcdf':

                    # method to organize dataset
                    variable_dset = organize_file_nc(
                        variable_collection,
                        obj_time=self.alg_time,
                        obj_variable_in=self.variables_in_dst['datasets'],
                        obj_variable_out=self.variables_out_dst['datasets'])

                    # method to write dataset
                    folder_name_dst, file_name_dst = os.path.split(file_path_dst)
                    make_folder(folder_name_dst)
                    write_file_nc(file_path_dst, variable_dset)

                elif self.format_dst == 'tiff' or self.format_dst == 'tif':

                    # method to organize dataset
                    variable_dset, attrs_dset = organize_file_tiff(
                        variable_collection,
                        obj_time=self.alg_time,
                        obj_variable_in=self.variables_in_dst['datasets'],
                        obj_variable_out=self.variables_out_dst['datasets'])

                    # method to write dataset
                    folder_name_dst, file_name_dst = os.path.split(file_path_dst)
                    make_folder(folder_name_dst)
                    write_file_tiff(file_path_dst, variable_dset, **attrs_dset)

                else:
                    # info end method
                    logging.info(' ---> Dump datasets ... FAILED. Destination format not allowed')
                    raise NotImplementedError('Destination format not implemented yet')

                # info end method
                logging.info(' ---> Dump datasets ... DONE')

            else:
                # info end method
                logging.info(' ---> Dump datasets ... SKIPPED. Data previously saved')
        else:
            # info end method
            logging.info(' ---> Dump datasets ... FAILED. Data grid not available')


# -------------------------------------------------------------------------------------




