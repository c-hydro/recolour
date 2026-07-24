"""
Library Features:

Name:           lib_utils_zip
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260715'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import gzip
import logging
import os
import shutil

from typing import Optional

from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to decompress gzip file
def unzip_file_gzip(
        file_path: str,
        file_path_uncompressed: Optional[str] = None,
        output_folder: Optional[str] = None,
        overwrite: bool = False,
        chunk_size: int = 1024 * 1024,
) -> str:
    """
    Decompress a gzip-compressed NetCDF file.

    Parameters
    ----------
    file_path : str
        Path of the input .gz file.

    file_path_uncompressed : str, optional
        Explicit output path.

    output_folder : str, optional
        Output folder used when file_path_uncompressed is not provided.

    overwrite : bool
        Overwrite an existing decompressed file.

    chunk_size : int
        Copy buffer size in bytes.

    Returns
    -------
    str
        Path of the decompressed NetCDF file.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Compressed file '{file_path}' was not found.")

    if not os.path.isfile(file_path):
        raise RuntimeError(f"Compressed path '{file_path}' is not a regular file.")

    if chunk_size <= 0:
        raise ValueError("The decompression chunk size must be greater than zero.")

    file_name = os.path.basename(file_path)

    if file_name.lower().endswith(".gz"):
        uncompressed_name = file_name[:-3]
    else:
        uncompressed_name = f"{file_name}.nc"

    if file_path_uncompressed is None:

        if output_folder is None:
            output_folder = os.path.dirname(file_path)
        os.makedirs(output_folder,exist_ok=True,)

        file_path_uncompressed = os.path.join(output_folder, uncompressed_name,)

    else:
        output_folder = os.path.dirname(file_path_uncompressed)

        if output_folder:
            os.makedirs(output_folder, exist_ok=True,)

    if os.path.exists(file_path_uncompressed) and not overwrite:
        logger.info(f" ===> Uncompressed NetCDF already exists: {file_path_uncompressed}.")
        return file_path_uncompressed

    temporary_path = f"{file_path_uncompressed}.tmp"

    if os.path.exists(temporary_path):
        os.remove(temporary_path)

    logger.info(f" ---------> Decompress NetCDF {file_path} ...")
    try:

        with gzip.open(file_path, mode="rb") as source_handle:
            with open(temporary_path, mode="wb") as destination_handle:
                shutil.copyfileobj(source_handle, destination_handle, length=chunk_size,)

        os.replace(temporary_path, file_path_uncompressed,)

    except Exception as exc:

        if os.path.exists(temporary_path):
            os.remove(temporary_path)

        raise RuntimeError(f"Unable to decompress NetCDF file '{file_path}'.") from exc

    logger.info(f" ---------> Decompress NetCDF {file_path} ... DONE")

    return file_path_uncompressed
# ----------------------------------------------------------------------------------------------------------------------
