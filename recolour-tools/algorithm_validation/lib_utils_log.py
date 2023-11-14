"""
Library Features:

Name:          lib_utils_log
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to open log file
def open_file_log(group_id, cell_id, file_path, file_update=True, execution_mode='NA',
                  group_fill_zeros=3, cell_fill_zeros=4):

    file_group = str(group_id).zfill(group_fill_zeros)
    file_cell = str(cell_id).zfill(cell_fill_zeros)

    file_tags = {'execution_mode': execution_mode,
                 'group_n': file_group, 'cell_n': file_cell}

    file_path = file_path.format(**file_tags)

    if file_update:
        if os.path.exists(file_path):
            os.remove(file_path)
    file_handle = open(file_path, 'w')

    return file_handle, file_group, file_cell

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to dump message to log
def dump_message_to_file(file_handle, line_arrow=' ===> ', line_msg='', line_sep=' ', line_break=True):
    if line_break:
        file_line = line_sep.join([line_arrow, line_msg, '\n'])
    else:
        file_line = line_sep.join([line_arrow, line_msg])
    file_handle.write(file_line)
    file_handle.flush()

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to close log file
def close_file_log(file_handle):
    file_handle.close()
# ----------------------------------------------------------------------------------------------------------------------
