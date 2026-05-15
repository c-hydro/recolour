# -*- coding: utf-8 -*-
# The MIT License (MIT)
#
# Copyright (c) 2016, TU Wien
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Module for a command line interface to convert the HSAF image into a
time series format using the repurpose package
'''

import numpy as np
import os
import sys
import argparse
from datetime import datetime

from orbit_resempler_wrapper import ASCATOrbitResampler

from ascat_custom_2018.h_saf import H16imgChunked, H101imgChunked, H102imgChunked, H103imgChunked

#from ascat.h_saf import H16img, H101img, H102img, H103img
from ascat_custom_2018.timeseries import load_grid
#from pynetcf.time_series import GriddedNcIndexedRaggedTs

from pynetcf_custom_2018.time_series import GriddedNcIndexedRaggedTs

def mkdate(datestring):
    if len(datestring) == 10:
        return datetime.strptime(datestring, '%Y-%m-%d')
    if len(datestring) == 16:
        return datetime.strptime(datestring, '%Y-%m-%dT%H:%M')

def resampler(product, ascat_input_path, output_path, grid_path,
              start_dt, end_dt,
              grid_filename='TUW_WARP5_grid_info_2_3.nc',
              ascat_month_path_str='%Y/%m/%d/%H/',
              writing_mode='w',
              write_n_resampled=1000):
    
    ascat_io_lut = {
        'H16': H16imgChunked,
        'H101': H101imgChunked,
        'H102': H102imgChunked,
        'H103': H103imgChunked,
    }
    
    #ascat_io_lut = {
    #    'H16': H16img,
    #    'H101': H101img,
    #    'H102': H102img,
    #    'H103': H103img,
    #}

    ascat_io_class = ascat_io_lut[product]
    ascat_io = ascat_io_class(ascat_input_path,
                              month_path_str=ascat_month_path_str,
                              chunk_minutes=50)
                              
    #ascat_io = ascat_io_class(ascat_input_path,
    #                          month_path_str=ascat_month_path_str)

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    target_grid = load_grid(os.path.join(grid_path, grid_filename))

    output_ts = GriddedNcIndexedRaggedTs(output_path, grid=target_grid,
                                         mode=writing_mode,
                                         ioclass_kws={'time_units': "days since 1858-11-17 00:00:00"})

    time_stamps = np.array(ascat_io.tstamps_for_daterange(start_dt, end_dt))

    resampler = ASCATOrbitResampler(ascat_io, output_ts,
                                    spatial_res=25000,
                                    dt=15,
                                    wfunc='hamming',
                                    write_orbit_buffer=True)

    resampler.resample(time_stamps, write_n_resampled=write_n_resampled)

    print('ciao')

def parse_args(args):
    """
    Parse command line parameters for conversion from image to timeseries

    :param args: command line parameters as list of strings
    :return: command line parameters as :obj:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(
        description="Convert swath data to time series format.")
    parser.add_argument("product", type=str,
                        help='Ascat satellite HSAF product (H16, H101, H103 or H103)')

    parser.add_argument("ascat_path", type=str,
                        help='Root of local filesystem where the data is stored.')
    parser.add_argument("timeseries_path", type=str,
                        help='Root of local filesystem where the timeseries should be stored.')

    parser.add_argument("grid_path", type=str,
                        help='Path to grid file.')

    parser.add_argument("start", type=mkdate,
                        help="Startdate. Either in format YYYY-MM-DD or YYYY-MM-DDTHH:MM.")
    parser.add_argument("end", type=mkdate,
                        help="Enddate. Either in format YYYY-MM-DD or YYYY-MM-DDTHH:MM.")

    parser.add_argument("--grid_filename", type=str, default='TUW_WARP5_grid_info_2_3.nc',
                        help='Name of grid file.')

    parser.add_argument("--ascat_month_path_str", type=str, default='%Y/%m/%d/%H/',
                        help="If the ASCAT BUFR files are stored in subfolders "
                             "of the ascat_input_path based on dates then "
                             "the rule of the folder names has to be given here. "
                             "The default is '%Y%m' which means that there are folders "
                             "like 201601 and 201602 and so on under ascat_input_path")

    parser.add_argument("--writing_mode", type=str, default='w',
                        help="'a' 'w' ")

    parser.add_argument("--write_n_resampled", type=int, default=1000,
                        help="Write to disk after having resampled this many swaths")

    args = parser.parse_args(args)
    # set defaults that can not be handled by argparse

    print("Converting data from {} to {} into folder {}.".format(args.start.isoformat(),
                                                                 args.end.isoformat(),
                                                                 args.timeseries_path))
    return args


def main(args):
    args = parse_args(args)

    resampler(args.product,
              args.ascat_path,
              args.timeseries_path,
              args.grid_path,
              args.start,
              args.end,
              grid_filename=args.grid_filename,
              ascat_month_path_str=args.ascat_month_path_str,
              writing_mode=args.writing_mode,
              write_n_resampled=args.write_n_resampled
              )

def run():
    main(sys.argv[1:])

if __name__ == '__main__':
    run()
