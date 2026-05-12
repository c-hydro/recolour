from datetime import datetime, timedelta
from ascat.swath import SwathGridFiles


swath_source = "/home/fabio/Desktop/recolour/dset/ascat_swath/h122/"
swath_collection = SwathGridFiles.from_product_id(swath_source, "H122_NRT")

# where to save the files
cell_file_directory = "/home/fabio/Desktop/recolour/exec/ascat_cells/h122/"

# the maximum size of the data buffer before dumping to file (actual maximum memory used will be higher)
max_nbytes = 150 * 1024**2  # 52428800 bytes 50 mbytes

# the date range to use. This should be a tuple of datetime.datetime objects
date_range = (datetime(2026, 5, 3, hour=7), datetime(2026, 5, 5, hour=8))

# Pass a list of cell numbers (integers) here if you only want to stack data for a certain set of cells. This is mainly useful for testing purposes, since even splitting a day's worth of swath data into files for all of its constituent cells is a lengthy process.
cells=None

# mode : "w" for creating new files if any already exist, "a" to append data to existing cell files
mode = "w"

# stack to cell
fmt_kwargs = {
    'dt_delta': timedelta(hours=1)
}

swath_collection.stack_to_cell_files(
    out_dir=cell_file_directory, max_nbytes=max_nbytes,
    date_range=date_range, fmt_kwargs=fmt_kwargs,
    print_progress=True, parallel=False)



