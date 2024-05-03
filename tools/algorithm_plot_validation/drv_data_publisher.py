"""
Class Features

Name:          drv_data_viewer
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import os
import shutil

from lib_utils_generic import make_folder
from lib_utils_datasets import organize_datasets_cell, organize_datasets_grid, convert_datasets_obj
from lib_data_io_pickle import read_file_obj, write_file_obj
from lib_data_io_nc import read_file_collection, write_file_collection

from lib_utils_io import filter_dframe_by_vars

from lib_data_statistics_pie import compute_stats_pearson as compute_stats_pearson_pie
from lib_data_statistics_pie import compute_stats_snr as compute_stats_snr_pie
from lib_data_statistics_box import compute_stats_pearson as compute_stats_pearson_box
from lib_data_statistics_box import compute_stats_snr as compute_stats_snr_box

from lib_figure_settings import organize_figure_settings, organize_figure_extent
from lib_figure_fx_results_generic import plot_committed_area
from lib_figure_fx_results_data import plot_results_data
from lib_figure_fx_results_pie import plot_results_pie_pearson, plot_results_pie_snr
from lib_figure_fx_results_box import plot_results_box_pearson, plot_results_box_snr

# figure expected parameters
default_params_figure = {
    "type": None, "committed_area": None,
    "variable_in_data": None, "variable_in_p_r": None, "variable_in_r": None,
    "variable_in": None, "variable_out": None,
    "y_label": None, "x_label": None, "title_label": None,
    "lim_threshold": None, "lim_target": None, "lim_optimal": None, "lim_min": None, "lim_max": None,
    "palette_type": None,
    "cbar_label": None, "cbar_extent": None, "cbar_ticks": None, "cbar_show": None,
    "vmin": None, "vmax": None, "clip": None,
    "cmap_type": None, "cmap_n": None
}
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# class driver data
class DrvData:

    # method to initialize class
    def __init__(self, alg_cell_list, alg_cell_grid, alg_settings,
                 tag_section_flags='flags', tag_section_info='info',
                 tag_section_domain='domain', tag_section_grid='grid',
                 tag_section_params='parameters',
                 tag_section_datasets='datasets',
                 tag_section_figure='figure', tag_section_renderer='renderer',
                 tag_section_time='time', tag_section_log='log'):

        self.alg_cells = alg_cell_list
        self.alg_cell_grid = alg_cell_grid

        self.alg_flags = alg_settings[tag_section_flags]
        self.alg_info = alg_settings[tag_section_info]
        self.alg_params_datasets = alg_settings[tag_section_params]['datasets']
        self.alg_params_figure = alg_settings[tag_section_params]['figure']
        self.alg_domain = alg_settings[tag_section_domain]
        # self.alg_time = alg_settings[tag_section_time]
        self.alg_grid = alg_settings[tag_section_grid]
        self.alg_datasets_src = alg_settings[tag_section_datasets]['source']
        self.alg_datasets_anc = alg_settings[tag_section_datasets]['ancillary']
        self.alg_datasets_dst = alg_settings[tag_section_datasets]['destination']
        self.alg_figure = alg_settings[tag_section_figure]
        self.alg_renderer = alg_settings[tag_section_renderer]
        self.alg_log = alg_settings[tag_section_log]

        self.tag_folder_name = 'folder_name'
        self.tag_file_name = 'file_name'
        self.tag_file_vars = 'variables'
        self.tag_file_mode = 'active'

        self.alg_datasets_name = self.alg_info['datasets']
        self.alg_datasets_product = self.alg_info['product']

        self.reset_datasets_anc = self.alg_flags['reset_ancillary_datasets']
        self.reset_datasets_dst = self.alg_flags['reset_destination_datasets']
        self.reset_logs = self.alg_flags['reset_logs']

        self.folder_name_grid = self.alg_grid[self.tag_folder_name]
        self.file_name_grid = self.alg_grid[self.tag_file_name]
        self.file_path_grid = os.path.join(self.folder_name_grid, self.file_name_grid)

        self.folder_name_src = self.alg_datasets_src[self.tag_folder_name]
        self.file_name_src = self.alg_datasets_src[self.tag_file_name]
        self.variables_src = self.alg_datasets_src[self.tag_file_vars]
        self.file_path_src = os.path.join(self.folder_name_src, self.file_name_src)

        self.folder_name_anc = self.alg_datasets_anc[self.tag_folder_name]
        self.file_name_anc = self.alg_datasets_anc[self.tag_file_name]
        self.file_path_anc = os.path.join(self.folder_name_anc, self.file_name_anc)

        self.folder_name_dst = self.alg_datasets_dst[self.tag_folder_name]
        self.file_name_dst = self.alg_datasets_dst[self.tag_file_name]
        self.variables_dst = self.alg_datasets_dst[self.tag_file_vars]
        self.file_path_dst = os.path.join(self.folder_name_dst, self.file_name_dst)

        self.figure_collection = self._join_figure_and_renderer()

        self._clean_ancillary_folders(self.alg_log,
                                      dset_key_root='path_log', dset_key_sub=None,
                                      dset_clean=self.reset_logs)

    # method to join figure(s) and renderer(s)
    def _join_figure_and_renderer(self):

        figure_collections = {}
        for figure_key, figure_fields in self.alg_figure.items():
            if not figure_key.startswith('_'):
                folder_name = figure_fields[self.tag_folder_name]
                file_name = figure_fields[self.tag_file_name]
                figure_variables = figure_fields[self.tag_file_vars]
                figure_mode = figure_fields[self.tag_file_mode]
                figure_path = os.path.join(folder_name, file_name)

                if figure_key in list(self.alg_renderer.keys()):
                    figure_renderer = self.alg_renderer[figure_key]
                    figure_settings = {**figure_renderer, **self.alg_params_figure}
                else:
                    logging.error(' ===> Figure renderer for "' + figure_key + '" does not exist')
                    raise RuntimeError('Publish data needs the renderer obj to plot the datasets')

                # update figure settings with default parameters (if not available)
                for default_key, default_value in default_params_figure.items():
                    if default_key not in list(figure_settings.keys()):
                        figure_settings[default_key] = default_value

                figure_collections[figure_key] = {}
                figure_collections[figure_key]['figure_mode'] = figure_mode
                figure_collections[figure_key]['figure_path_file'] = figure_path
                figure_collections[figure_key]['figure_variables'] = figure_variables
                figure_collections[figure_key]['figure_settings'] = figure_settings

        return figure_collections

    # method to check datasets (if clean or not)
    def _clean_ancillary_folders(self, dset_obj,
                                 dset_key_root=None, dset_key_sub=None, dset_clean=True):

        if (dset_key_root is not None) and (dset_key_sub is not None):
            for dset_key, dset_fields in dset_obj.items():
                dset_path = dset_fields[dset_key_root][dset_key_sub]
                self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is not None):
            dset_path = dset_obj[dset_key_sub]
            self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is not None) and (dset_key_sub is None):
            dset_path = dset_obj[dset_key_root]
            self._remove_datasets(dset_path, dset_clean=dset_clean)
        elif (dset_key_root is None) and (dset_key_sub is None):
            pass
        else:
            pass

    # method to remove datasets
    @staticmethod
    def _remove_datasets(dset_path, dset_clean=True):
        if os.path.exists(dset_path):
            if dset_clean:
                shutil.rmtree(dset_path)
        os.makedirs(dset_path, exist_ok=True)

    # method to organize data
    def organize_data(self):

        # info start method
        logging.info(' ---> Organize datasets ... ')

        # get cells
        alg_cells = self.alg_cells
        # get file variables
        variables_src, variables_dst = self.variables_src, self.variables_dst
        # get file path
        file_path_src, file_path_anc, file_path_dst = self.file_path_src, self.file_path_anc, self.file_path_dst

        # clean ancillary and destination datasets (if ancillary flag is activated
        if self.reset_datasets_anc:
            if os.path.exists(file_path_anc):
                os.remove(file_path_anc)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)
        # clean ancillary and destination datasets if are not available together
        if (not os.path.exists(file_path_anc)) or (not os.path.exists(file_path_dst)):
            if os.path.exists(file_path_anc):
                os.remove(file_path_anc)
            if os.path.exists(file_path_dst):
                os.remove(file_path_dst)

        # check cell obj availability
        if not os.path.exists(file_path_anc):

            # join cell object(s)
            cell_obj_dict = organize_datasets_cell(
                cell_list=alg_cells, cell_digits=4, cell_grid=self.alg_cell_grid,
                list_variable_in=variables_src, list_variable_out=variables_dst,
                folder_name_datasets=self.folder_name_src, file_name_datasets=self.file_name_src,
                )
            # filter cell object
            cell_obj_dict = organize_datasets_grid(
                cell_obj_dict, data_type=self.alg_datasets_name,
                file_name_grid=self.file_path_grid)

            # save cell obj in pickle format
            folder_name_anc, file_name_anc = os.path.split(file_path_anc)
            make_folder(folder_name_anc)
            write_file_obj(file_path_anc, cell_obj_dict)
            # save cell obj in netcdf format
            folder_name_dst, file_name_dst = os.path.split(file_path_dst)
            make_folder(folder_name_dst)
            write_file_collection(file_path_dst, cell_obj_dict)

        else:
            # read cell obj from pickle format
            cell_obj_dict = read_file_obj(file_path_anc)

        # convert from dict to dframe
        cell_obj_dframe = convert_datasets_obj(cell_obj_dict)
        # info end method
        logging.info(' ---> Organize datasets ... DONE')

        return cell_obj_dframe

    # method to plot data
    def plot_data(self, figure_dframe_tmp):

        # info start method
        logging.info(' ---> Plot datasets ... ')

        # get cells
        alg_cells = self.alg_cells
        # get file variables
        file_path_grid = self.file_path_grid
        file_path_dst = self.file_path_dst
        datasets_name = self.alg_datasets_name
        variables_dst = self.variables_dst
        # get figure collection
        figure_collection = self.figure_collection

        # method to get figure dataframe
        variables_dframe = read_file_collection(
            file_path_dst, file_path_grid,
            variable_type=datasets_name,
            variable_list=variables_dst)

        # iterate over figure(s)
        for figure_key, figure_fields in figure_collection.items():

            # info start plot
            logging.info(' ----> Image "' + figure_key + '" ... ')

            # get figure args
            figure_path_file = figure_fields['figure_path_file']
            figure_variables = figure_fields['figure_variables']
            figure_settings = figure_fields['figure_settings']
            figure_mode = figure_fields['figure_mode']

            # check figure mode active/not active
            if figure_mode:

                # get figure settings
                figure_type = figure_settings['type']
                figure_var_in, figure_var_out = figure_settings['variable_in'], figure_settings['variable_out']
                figure_var_in_data = figure_settings['variable_in_data']
                figure_var_in_p_r, figure_var_in_r = figure_settings['variable_in_p_r'], figure_settings['variable_in_r']

                figure_title_label, figure_cbar_label = figure_settings['title_label'], figure_settings['cbar_label']
                figure_cbar_extent, figure_cbar_ticks = figure_settings['cbar_extent'], figure_settings['cbar_ticks']
                figure_vmin, figure_vmax = figure_settings['vmin'], figure_settings['vmax']
                figure_cmap_type, figure_cmap_n = figure_settings['cmap_type'], figure_settings['cmap_n']
                figure_data_extent_default = figure_settings['data_extent_default']
                figure_data_extent_over_domain = figure_settings['data_extent_over_domain']
                figure_title_fontsize = figure_settings['title_fontsize']
                figure_title_fontweight = figure_settings['title_fontweight']
                figure_cbar_fontsize = figure_settings['cb_fontsize']
                figure_cbar_fontweight = figure_settings['cb_fontweight']
                figure_cbar_show = figure_settings['cbar_show']
                figure_committed_area = figure_settings['committed_area']
                figure_lim_thr, figure_lim_target = figure_settings['lim_threshold'], figure_settings['lim_target']
                figure_lim_optimal = figure_settings['lim_optimal']
                figure_lim_min, figure_lim_max = figure_settings['lim_min'], figure_settings['lim_max']
                figure_x_label, figure_y_label = figure_settings['x_label'], figure_settings['y_label']
                figure_palette_type = figure_settings['palette_type']

                # filter dataframe using expected variable(s)
                if figure_type == 'stats_pearson_box':
                    figure_dframe = filter_dframe_by_vars(
                        variables_dframe,
                        dframe_parameter=[figure_var_in_p_r, figure_var_in_r],
                        dframe_variables_data=figure_variables, dframe_variables_grid=['land_flag', 'committed_area'])
                elif figure_type == 'stats_snr_box':
                    figure_dframe = filter_dframe_by_vars(
                        variables_dframe, dframe_parameter=[figure_var_in_data, figure_var_in_p_r, figure_var_in_r],
                        dframe_variables_data=figure_variables, dframe_variables_grid=['land_flag', 'committed_area'])
                else:
                    figure_dframe = filter_dframe_by_vars(
                        variables_dframe, dframe_parameter=figure_var_in,
                        dframe_variables_data=figure_variables, dframe_variables_grid=['land_flag', 'committed_area'])

                # organize figure extent
                figure_data_extent_def = organize_figure_extent(
                    figure_dframe, figure_data_extent_default, fig_ext_domain=figure_data_extent_over_domain)

                # select figure parameter
                if figure_type == 'data':
                    figure_parameter = figure_var_in
                elif figure_type == 'stats_pearson_pie' or figure_type == 'stats_snr_pie':
                    figure_parameter = figure_var_out
                elif figure_type == 'stats_pearson_box' or figure_type == 'stats_snr_box':
                    figure_parameter = figure_var_out
                elif figure_type == 'committed_area':
                    figure_parameter = figure_var_in
                else:
                    logging.error(' ===> Figure type "' + figure_type + '" is not expected')
                    raise NotImplemented('Case not implemented yet')

                # set figure settings as expected by the plotting fx
                figure_file_obj, figure_flag_carea_obj, figure_kwargs = organize_figure_settings(
                    figure_settings, figure_type=figure_type, figure_committed_area=figure_committed_area,
                    figure_filename_tmpl=figure_path_file,
                    figure_parameter=figure_parameter,
                    figure_parameter_data=figure_var_in_data, figure_parameter_p_r=figure_var_in_p_r,
                    figure_parameter_r=figure_var_in_r,
                    figure_title_label=figure_title_label,
                    figure_x_label=figure_x_label, figure_y_label=figure_y_label,
                    figure_data_extent=figure_data_extent_def,
                    figure_colorbar_label=figure_cbar_label, figure_colorbar_show=figure_cbar_show,
                    figure_colorbar_extent=figure_cbar_extent, figure_colorbar_ticks=figure_cbar_ticks,
                    figure_vmin=figure_vmin, figure_vmax=figure_vmax,
                    figure_cmap_type=figure_cmap_type, figure_cmap_n=figure_cmap_n,
                    figure_title_fontsize=figure_title_fontsize, figure_title_fontweight=figure_title_fontweight,
                    figure_cbar_fontsize=figure_cbar_fontsize, figure_cbar_fontweight=figure_cbar_fontweight,
                    figure_lim_min=figure_lim_min, figure_lim_max=figure_lim_max,
                    figure_lim_thr=figure_lim_thr,
                    figure_lim_target=figure_lim_target, figure_lim_optimal=figure_lim_optimal,
                    figure_palette_type=figure_palette_type
                    )

                # select plotting fx according to expected data
                if figure_type == 'data':

                    # plot figure data
                    for figure_file_name, figure_carea_flag in zip(figure_file_obj, figure_flag_carea_obj):
                        plot_results_data(figure_file_name, figure_dframe, figure_kwargs, figure_carea_flag)

                elif figure_type == 'stats_snr_pie':

                    # compute statistics snr
                    figure_dframe, figure_perc_comm, figure_perc_global = compute_stats_snr_pie(
                        figure_dframe, variable_data=figure_var_in, variable_stats=figure_var_out)
                    # plot figure stats/pie
                    plot_results_pie_snr(figure_file_obj, figure_perc_comm, figure_perc_global, figure_kwargs)

                elif figure_type == 'stats_pearson_pie':

                    # compute statistics pearson
                    figure_dframe, figure_perc_comm, figure_perc_global = compute_stats_pearson_pie(
                        figure_dframe, variable_data=figure_var_in, variable_stats=figure_var_out)
                    # plot figure stats/pie
                    plot_results_pie_pearson(figure_file_obj, figure_perc_comm, figure_perc_global, figure_kwargs)

                elif figure_type == 'stats_snr_box':

                    # compute statistics snr
                    figure_df_obj, figure_perc_obj = compute_stats_snr_box(
                        figure_dframe,
                        variable_data=figure_var_in_data, variable_r=figure_var_in_r, variable_p_r=figure_var_in_p_r,
                        variable_stats=figure_var_out)
                    # plot figure stats/box
                    plot_results_box_snr(figure_file_obj, figure_df_obj, figure_perc_obj, figure_kwargs)

                elif figure_type == 'stats_pearson_box':

                    # compute statistics pearson
                    figure_df_obj, figure_perc_obj = compute_stats_pearson_box(
                        figure_dframe,
                        variable_data=figure_var_in_data, variable_r=figure_var_in_r, variable_p_r=figure_var_in_p_r,
                        variable_stats=figure_var_out)
                    # plot figure stats/box
                    plot_results_box_pearson(figure_file_obj, figure_df_obj, figure_perc_obj, figure_kwargs)

                elif figure_type == 'committed_area':

                    # plot figure committed area
                    plot_committed_area(figure_file_obj, figure_dframe, figure_kwargs)

                else:
                    logging.error(' ===> Figure type "' + figure_type + '" is not expected')
                    raise NotImplemented('Case not implemented yet')

                # info end plot
                logging.info(' ----> Image "' + figure_key + '" ... DONE')

            else:

                # info end plot
                logging.info(' ----> Image "' + figure_key + '" ... SKIPPED. PLOT NOT ACTIVATED')

        # info end method
        logging.info(' ---> Plot datasets ... DONE')

# -------------------------------------------------------------------------------------
