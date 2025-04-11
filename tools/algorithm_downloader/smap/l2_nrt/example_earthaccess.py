

import wget
import os
import requests
from bs4 import BeautifulSoup

import subprocess
import os
import requests
from netrc import netrc
from bs4 import BeautifulSoup


def get_credentials(machine='urs.earthdata.nasa.gov'):
    """Retrieve credentials from .netrc file."""
    auth = netrc().authenticators(machine)
    if auth:
        login, _, password = auth
        return login, password
    else:
        raise ValueError(f"No credentials found for {machine} in .netrc")

def list_files_in_remote_folder(url, file_ext='.h5'):

    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    files = [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(file_ext)]
    files = list(set(files))  # Remove duplicates

    return files


def download_file_with_auth(cmd_root, file_urls, download_path, machine='urs.earthdata.nasa.gov'):

    # Check if the download path exists, create if not
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    # Download each file
    file_list = []
    for file_url in file_urls:

        file_name = os.path.basename(file_url)
        file_list.append(file_name)

        file_path = os.path.join(download_path, file_name)

        if not os.path.exists(file_path):

            cmd_file = cmd_root + ' ' + file_url

            # Execute the command
            subprocess.run(cmd_file, shell=True, check=True)
            print(f"Downloaded file from {file_url} to {download_path}")

        else:
            print(f"File {file_name} already exists in {download_path}, skipping download.")

    return file_list


if __name__ == "__main__":

    bbox = [-125.0, 32.0, -113.0, 42.0]

    remote_address = 'https://n5eil01u.ecs.nsidc.org/'
    remote_folder_root = '/SMAP/SPL2SMP_E.006/'
    remote_folder_date = '2025.04.09/'
    remote_file = 'SMAP_L2_SM_P_E_54415_D_20250409T055838_R19240_001.h5'

    remote_url = ''.join([remote_address, remote_folder_root, remote_folder_date])

    # Path where the files will be saved
    save_path = "/home/fabio/Desktop/Recolour_Workspace/ws/smap_cell_nrt/spl2smp_e/"

    user, password = get_credentials(machine='urs.earthdata.nasa.gov')

    print(f"Remote URL: {remote_url}")

    wget_cmd = 'wget --user {user} --password {password} -P {save_path} -nH --cut-dirs 4 --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off'
    wget_cmd = wget_cmd.format(user=user, password=password, save_path=save_path)

    # List and download files
    file_urls = list_files_in_remote_folder(remote_url)
    file_list = download_file_with_auth(wget_cmd, file_urls, save_path)
