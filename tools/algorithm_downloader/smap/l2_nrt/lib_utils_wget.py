"""
Library Features:

Name:          lib_utils_wget
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)

Date:          '20250410'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import subprocess
import os
import requests

from bs4 import BeautifulSoup
from netrc import netrc

from lib_data_io_pickle import read_obj, write_obj
from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)

# add suppress warnings
urllib3_stream = logging.getLogger('urllib3')
urllib3_stream.setLevel(logging.ERROR)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get credentials
def get_credentials(user: str = None, password: str = None, machine: str ='urs.earthdata.nasa.gov'):

    if user is None or password is None:
        try:
            auth = netrc().authenticators(machine)
            if auth:
                login, _, password = auth
                return login, password
            else:
                logger_stream.error(f" ===> No credentials found for {machine} in .netrc")
                raise ValueError(f"Credentials are mandatory. Please provide them.")
        except Exception as e:
            logger_stream.error(f" ===> Error retrieving credentials: {e}")
            raise RuntimeError(f"Credentials are mandatory. Please provide them.")
    elif user is not None and password is not None:
        return user, password
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to create file url
def create_file_url(remote_url: str = 'https://n5eil01u.ecs.nsidc.org/',
                    remote_folder: str = '/SMAP/SPL2SMP_E.006/2025.04.09/') -> str:
    remote_url = ''.join([remote_url, remote_folder])
    return remote_url
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to search for files in a remote folder
def search_file_url(remote_url: str, file_ext: str ='.h5',
                    file_lock: str = "remote_files.pkl", remove_lock: bool = True) -> list:

    if remove_lock:
        if os.path.exists(file_lock):
            os.remove(file_lock)

    # for debugging
    if not os.path.exists(file_lock):

        response = requests.get(remote_url)
        soup = BeautifulSoup(response.text, 'html.parser')

        remote_files = [remote_url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(file_ext)]
        remote_files = list(set(remote_files))  # Remove duplicates

        write_obj(file_lock, remote_files)

    else:

        remote_files = read_obj(file_lock)

    return remote_files
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# define the command to download files
def create_file_wget(user: str, password: str, save_path: str):

    wget_cmd = 'wget --user {user} --password {password} -P {save_path} -nH --cut-dirs 4 --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off'
    wget_cmd = wget_cmd.format(user=user, password=password, save_path=save_path)

    return wget_cmd
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to define file url
def download_file_url(cmd_root, file_url, path_download, machine='urs.earthdata.nasa.gov'):

    # check if the download path exists, create if not
    if not os.path.exists(path_download):
        os.makedirs(path_download)

    # info start
    logger_stream.info(f" ------> Downloaded file from {file_url} to {path_download} ... ")

    # organize file name, path
    file_name = os.path.basename(file_url)
    file_path = os.path.join(path_download, file_name)

    # check if the file already exists
    if not os.path.exists(file_path):

        # create the command
        cmd_file = cmd_root + ' ' + file_url
        # execute the command
        subprocess.run(cmd_file, shell=True, check=True)

        # info end
        logger_stream.info(f" ------> Downloaded file from {file_url} to {path_download} ... DONE")

    else:
        # info end
        logger_stream.info(f" ------> Downloaded file from {file_url} to {path_download} ... SKIPPED. "
                     f"File already exists")

    return file_path
# ----------------------------------------------------------------------------------------------------------------------
