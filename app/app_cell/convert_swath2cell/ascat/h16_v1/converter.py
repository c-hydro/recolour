
from os.path import join, exists
from os import makedirs
import re
from os import walk
from numpy import where, argsort, concatenate, array, asarray, ones, delete
from netCDF4 import Dataset
#import progressbar

def search_cell(path, rexp=r'^\d{4}'):

    cells = []
    for root, dirs, files in walk(path):
        for file in files:
            if re.search(rexp, file):
                cell = re.search(rexp, file).group()
                cell = cell.zfill(4)
                cells.append(cell)

    return sorted(cells)

class Ind2ContTs(object):

    dset_in = None
    dset_out = None

    workspace_file_attrs = {}

    workspace_ts_data = {}
    workspace_ts_attrs = {}

    def __init__(self,
                 filename_in='{:}.nc', filepath_in=None, filename_in_format='{:4}',
                 filename_out='{:}_upd.nc', filepath_out=None, filename_out_format='{:4}',
                 loc_dim_name='locations', obs_dim_name='obs',
                 loc_ids_var='location_id',
                 loc_descr_var='location_description',
                 loc_idx_var='locationIndex',
                 loc_lut_var='row_size',
                 time_var='time', time_units="days since 1858-11-17 00:00:00",
                 lat_var='lat',
                 lon_var='lon',
                 alt_var='alt',
                 **kwargs):

        self.filepath_in = filepath_in

        if filepath_out is None:
            filepath_out = filepath_in
        self.filepath_out = filepath_out

        if not exists(filepath_out):
            makedirs(filepath_out)

        self.filename_in = filename_in
        self.filename_out = filename_out

        self.filename_in_format = filename_in_format
        self.filename_out_format = filename_out_format

        # dimension names
        self.obs_dim_name = obs_dim_name
        self.loc_dim_name = loc_dim_name

        # location names
        self.loc_ids_var = loc_ids_var
        self.loc_descr_var = loc_descr_var
        self.loc_idx_var = loc_idx_var
        self.loc_lut_var = loc_lut_var

        # time, time units and location
        self.time_var = time_var
        self.time_units = time_units
        self.lat_var = lat_var
        self.lon_var = lon_var
        self.alt_var = alt_var

        self.var_not_ts = [loc_ids_var, loc_descr_var, loc_idx_var, loc_lut_var,
                           time_var,
                           lat_var, lon_var, alt_var]

        self.workspace_file_attrs = None
        self.workspace_ts_data = None
        self.workspace_ts_attrs = None

        self.cells = None

    def __read_idx_ts(self):

        # Initialize workspace(s)
        self.workspace_file_attrs = {}
        self.workspace_ts_data = {}
        self.workspace_ts_attrs = {}

        # Get global attributes
        for attr_name in self.dset_in.ncattrs():
            self.workspace_file_attrs[attr_name] = getattr(self.dset_in, attr_name)

        # Get registry data
        gpis_loc = self.dset_in.variables[self.loc_ids_var][:]
        lons_loc = self.dset_in.variables[self.lon_var][:]
        lats_loc = self.dset_in.variables[self.lat_var][:]
        alts_loc = self.dset_in.variables[self.alt_var][:]

        gpis_nan = []
        for i, gpi in enumerate(gpis_loc):
            if not gpi:
                gpis_nan.append(i)

        gpis_nan = asarray(gpis_nan, dtype=int)

        if gpis_nan.size > 0:
            gpis = delete(gpis_loc, gpis_nan)
            lons = delete(lons_loc, gpis_nan)
            lats = delete(lats_loc, gpis_nan)
            alts = delete(alts_loc, gpis_nan)
        else:
            gpis = gpis_loc
            lons = lons_loc
            lats = lats_loc
            alts = alts_loc

        # Get indexed variable(s)
        idx = self.dset_in.variables[self.loc_idx_var][:]
        time = self.dset_in.variables[self.time_var][:]

        for var_name in self.dset_in.variables:

            var_name = str(var_name)
            var_data = self.dset_in.variables[var_name][:]

            var_attrs = {}
            for attr_name in self.dset_in.variables[var_name].ncattrs():
                var_attrs[attr_name] = self.dset_in.variables[var_name].getncattr(attr_name)
            self.workspace_ts_attrs[var_name] = var_attrs

            if var_name not in self.var_not_ts:

                var_data_ts = array([])
                time_ts = array([])
                row_ts = array([])

                for i, gpi in enumerate(gpis):

                    if gpi:
                        gpi_idx = where(idx == i)

                        var_data_idx = var_data[gpi_idx]
                        time_idx = time[gpi_idx]

                        sort_idx = argsort(time_idx)

                        var_data_sort = var_data_idx[sort_idx]
                        time_sort = time_idx[sort_idx]

                        time_ts = concatenate([time_ts, time_sort])
                        row_ts = concatenate([row_ts, asarray(var_data_sort.__len__() * ones(1))])
                        var_data_ts = concatenate([var_data_ts, var_data_sort])

                self.workspace_ts_data[var_name] = var_data_ts

                if not 'time' in self.workspace_ts_data:
                    self.workspace_ts_data['time'] = time_ts
                    self.workspace_ts_attrs['time'] = \
                        {'standard_name': 'time',
                         'long_name': 'time of measurement',
                         'units': self.time_units}

                if not self.loc_lut_var in self.workspace_ts_data:
                    self.workspace_ts_data[self.loc_lut_var] = asarray(row_ts, dtype=int)
                    self.workspace_ts_attrs[self.loc_lut_var] = \
                        {'long_name': 'number of observations at this location',
                         'sample_dimension': idx.__len__()}

                if not self.lon_var in self.workspace_ts_data:
                    self.workspace_ts_data[self.lon_var] = lons
                if not self.lat_var in self.workspace_ts_data:
                    self.workspace_ts_data[self.lat_var] = lats
                if not self.alt_var in self.workspace_ts_data:
                    self.workspace_ts_data[self.alt_var] = alts
                if not self.loc_ids_var in self.workspace_ts_data:
                    self.workspace_ts_data[self.loc_ids_var] = gpis

        self.dset_in.close()

    def __write_contiguous_ts(self):

        # Create variables dimension(s)
        loc_n = self.workspace_ts_data[self.loc_ids_var].__len__()
        obs_n = self.workspace_ts_data[self.time_var].__len__()

        self.dset_out.createDimension(self.loc_dim_name, loc_n)
        self.dset_out.createDimension(self.obs_dim_name, obs_n)

        for attr_name in self.workspace_file_attrs:
            self.dset_out.setncattr(attr_name, self.workspace_file_attrs[attr_name])

        for var_name in self.workspace_ts_data:

            var_data_ts = self.workspace_ts_data[var_name]
            if var_name in self.workspace_ts_attrs:
                var_attrs = self.workspace_ts_attrs[var_name]
            else:
                var_attrs = {}

            var_dtype = var_data_ts.dtype

            if isinstance(var_data_ts[0], (str, bytes)):
                continue

            #if isinstance(var_data_ts[0], unicode):
            #    continue

            if var_data_ts.__len__() == loc_n:
                var_dim = self.loc_dim_name
            elif var_data_ts.__len__() == obs_n:
                var_dim = self.obs_dim_name
            else:
                var_dim = None

            var_obj = self.dset_out.createVariable(
                var_name, var_dtype, dimensions=var_dim, fill_value=None, zlib=None, complevel=None)

            for attr_name in var_attrs:
                var_obj.setncattr(attr_name, var_attrs[attr_name])

            if var_data_ts is not None:
                var_obj[:] = var_data_ts

        self.dset_out.close()

    def conversion(self, cells=None):

        # Get cell
        if cells is None:
            self.cells = search_cell(self.filepath_in)
        else:
            self.cells = cells

        #bar = progressbar.ProgressBar()
        for cell in self.cells:

            print(f' cell {cell} ...' )

            # Open and get dset_in
            filename_in = self.filename_in.format(self.filename_in_format.format(cell))
            self.dset_in = Dataset(join(self.filepath_in, filename_in), mode='r')
            # Open and get dset_out
            filename_out = self.filename_out.format(self.filename_out_format.format(cell))
            self.dset_out = Dataset(join(self.filepath_out, filename_out), mode='w')

            self.__read_idx_ts()
            self.__write_contiguous_ts()

            print(f' cell {cell} ... DONE')

class H16_Ind2ContTS(Ind2ContTs):

    def __init__(self, filename_in_fmt='{:}.nc', filepath_in=None,
                 filename_out_fmt='h16_{:}.nc', filepath_out=None):

        self.filename_in_fmt = filename_in_fmt
        self.filename_out_fmt = filename_out_fmt

        self.filepath_in = filepath_in
        self.filepath_out = filepath_out

        super(H16_Ind2ContTS, self).__init__(
                                             filename_in=self.filename_in_fmt,
                                             filepath_in=self.filepath_in,
                                             filename_out=self.filename_out_fmt,
                                             filepath_out=self.filepath_out)

class H101_Ind2ContTS(Ind2ContTs):

    def __init__(self, filename_in_fmt='{:}.nc', filepath_in=None,
                 filename_out_fmt='h101_{:}.nc', filepath_out=None):

        self.filename_in_fmt = filename_in_fmt
        self.filename_out_fmt = filename_out_fmt

        self.filepath_in = filepath_in
        self.filepath_out = filepath_out

        super(H101_Ind2ContTS, self).__init__(
                                             filename_in=self.filename_in_fmt,
                                             filepath_in=self.filepath_in,
                                             filename_out=self.filename_out_fmt,
                                             filepath_out=self.filepath_out)

class H102_Ind2ContTS(Ind2ContTs):

    def __init__(self, filename_in_fmt='{:}.nc', filepath_in=None,
                 filename_out_fmt='h102_{:}.nc', filepath_out=None):

        self.filename_in_fmt = filename_in_fmt
        self.filename_out_fmt = filename_out_fmt

        self.filepath_in = filepath_in
        self.filepath_out = filepath_out
        
        super(H102_Ind2ContTS, self).__init__(
                                             filename_in=self.filename_in_fmt,
                                             filepath_in=self.filepath_in,
                                             filename_out=self.filename_out_fmt,
                                             filepath_out=self.filepath_out)

class H103_Ind2ContTS(Ind2ContTs):

    def __init__(self, filename_in_fmt='{:}.nc', filepath_in=None,
                 filename_out_fmt='h103_{:}.nc', filepath_out=None):

        self.filename_in_fmt = filename_in_fmt
        self.filename_out_fmt = filename_out_fmt

        self.filepath_in = filepath_in
        self.filepath_out = filepath_out
        
        super(H103_Ind2ContTS, self).__init__(
                                             filename_in=self.filename_in_fmt,
                                             filepath_in=self.filepath_in,
                                             filename_out=self.filename_out_fmt,
                                             filepath_out=self.filepath_out)
