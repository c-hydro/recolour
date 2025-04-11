"""
Library Features:

Name:          lib_utils_cmr
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
               Martina Natali (martina01.natali@edu.unife.it)
Date:          '20231110'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# libraries
import logging
import base64
import itertools
import json
import netrc
import os
import ssl
import sys
import numpy as np
from multiprocessing import Pool, cpu_count

from urllib.parse import urlparse
from urllib.request import urlopen, Request, build_opener, HTTPCookieProcessor
from urllib.error import HTTPError, URLError

from lib_info_args import logger_name

# logger stream
logger_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get credentials
def get_credentials(url, urs_url='https://urs.earthdata.nasa.gov'):

    credentials = None
    try:
        info = netrc.netrc()
        username, account, password = info.authenticators(urlparse(urs_url).hostname)
        errprefix = 'File netrc error: '; print(username, account, password)
    except Exception as e:
        logger_stream.error(' ===> File netrc error: {0}'.format(str(e)))
        raise RuntimeError('Credentials are not available on netrc file')

    while not credentials:
        credentials = '{0}:{1}'.format(username, password)
        credentials = base64.b64encode(credentials.encode('ascii')).decode('ascii')

        if url:
            try:
                req = Request(url)
                req.add_header('Authorization', 'Basic {0}'.format(credentials))
                opener = build_opener(HTTPCookieProcessor())
                opener.open(req)
            except HTTPError:
                logger_stream.error(' ===> ' + errprefix + 'Incorrect username or password [' + url + ']')
                raise RuntimeError('Credentials are not available on netrc file')

    return credentials
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to correct version string
def build_version_query_params(version):
    desired_pad_length = 3
    if len(version) > desired_pad_length:
        logger_stream.error(' ===> Version string too long: "{0}"'.format(version))
        raise RuntimeError('String is not allowed in this format')

    version = str(int(version))  # Strip off any leading zeros
    query_params = ''

    while len(version) <= desired_pad_length:
        padded_version = version.zfill(desired_pad_length)
        query_params += '&version={0}'.format(padded_version)
        desired_pad_length -= 1
    return query_params
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to build query url
# This method has the proper format of the url to query the NSIDC API
# If more options are needed, refer to
# https://nsidc.org/data/user-resources/help-center/table-key-value-pair-kvp-operands-subsetting-reformatting-and-reprojection-services
def build_cmr_query_url(short_name, version, time_start, time_end,
                        bounding_box=None, polygon=None,
                        filename_filter=None, cmr_file_url=None):
    params = '&short_name={0}'.format(short_name)
    params += build_version_query_params(version)
    params += '&temporal[]={0},{1}'.format(time_start, time_end)
    if polygon:
        params += '&polygon={0}'.format(polygon)
    elif bounding_box:
        params += '&bounding_box={0}'.format(bounding_box)
    if filename_filter:
        option = '&options[producer_granule_id][pattern]=true'
        params += '&producer_granule_id[]={0}{1}'.format(filename_filter, option)
    return cmr_file_url + params
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to execute url request
def cmr_request(data_list):

    src_data = data_list[0]
    dest_data = data_list[1]
    credentials = data_list[2]

    file_name = os.path.split(dest_data)[1]

    logging.info(' -----> Downloading ' + src_data + ' :: ' + file_name + ' ... ')

    if not os.path.exists(dest_data):
        try:
            req = Request(src_data)
            if credentials:
                req.add_header('Authorization', 'Basic {0}'.format(credentials))
            opener = build_opener(HTTPCookieProcessor())
            data = opener.open(req).read()
            open(dest_data, 'wb').write(data)

            logging.info(' -----> Downloading ' + src_data + ' :: ' + file_name + ' ... DONE')

        except HTTPError as e:
            logging.info(' -----> Downloading ' + src_data + ' :: ' + file_name + ' ... FAILED')
            logging.error(' ===> HTTP error {0}, {1}'.format(e.code, e.reason))

        except URLError as e:
            logging.info(' -----> Downloading ' + src_data + ' :: ' + file_name + ' ... FAILED')
            logging.error(' ===> URL error: {0}'.format(e.reason))
        except IOError:
            raise
        except KeyboardInterrupt:
            quit()
    else:
        logging.info(' -----> Downloading ' + src_data + ' :: ' + file_name + ' ... SKIPPED. PREVIOUSLY DONE')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to download url(s) in multiprocessing mode
def cmr_download_mp(urls, dests, process_n=4, process_max=None):

    if not urls:
        return

    if process_max is None:
        process_max = cpu_count() - 1
    if process_n > process_max:
        logging.warning(' ===> Maximum of recommended processes must be less then ' + str(process_max))
        logging.warning(' ===> Set number of process from ' + str(process_n) + ' to ' + str(process_max))
        process_n = process_max

    url_count = len(urls)
    logging.info(' ----> Transferring {0} files in multiprocessing mode ... '.format(url_count))
    credentials = None

    dest_file_list = list(dests.values())

    logging.info(' -----> Preparing files ... ')
    request_list = []
    for index, (url, dest_file) in enumerate(zip(urls, dest_file_list), start=1):

        if isinstance(dest_file, list):
            dest_file = dest_file[0]
        path_name, file_name = os.path.split(dest_file)
        os.makedirs(path_name, exist_ok=True)

        if not credentials and urlparse(url).scheme == 'https':
            credentials = get_credentials(url)

        if not os.path.exists(dest_file):
            request_list.append([url, dest_file, credentials])

    logging.info(' -----> Preparing files ... DONE')

    logging.info(' -----> Pooling requests ... ')
    response_list = []
    if request_list:
        with Pool(processes=process_n, maxtasksperchild=1) as process_pool:
            src_response = process_pool.map(cmr_request, request_list, chunksize=1)
            process_pool.close()
            process_pool.join()
            response_list.append(src_response)
        logging.info(' -----> Pooling requests ... DONE')
    else:
        logging.info(' -----> Pooling requests ... SKIPPED. PREVIOUSLY DONE')

    logging.info(' ----> Transferring {0} files in multiprocessing mode ... DONE'.format(url_count))
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to download url(s) in sequential mode
def cmr_download_seq(urls, dests):

    if not urls:
        return

    url_count = len(urls)
    logging.info(' ----> Transferring {0} files in sequential mode ... '.format(url_count))
    credentials = None

    dest_file_list = dests.values()

    for index, (url, dest_file) in enumerate(zip(urls, dest_file_list), start=1):

        if isinstance(dest_file, list):
            dest_file = dest_file[0]

        path_name, file_name = os.path.split(dest_file)
        os.makedirs(path_name, exist_ok=True)

        logging.info(' -----> {0}/{1}: {2} ... '.format(str(index).zfill(len(str(url_count))), url_count, file_name))

        if not os.path.exists(dest_file):
            if not credentials and urlparse(url).scheme == 'https':
                credentials = get_credentials(url)
            cmr_request([url, dest_file, credentials])
        else:
            logging.info(
                ' -----> {0}/{1}: {2} ... SKIPPED. PREVIOUSLY DONE'.format(str(index).zfill(len(str(url_count))),
                                                                           url_count, file_name))

    logging.info(' ----> Transferring {0} files in sequential mode ... DONE'.format(url_count))
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to filter url(s)
def cmr_filter_urls(search_results):
    """Select only the desired data files from CMR response."""
    if 'feed' not in search_results or 'entry' not in search_results['feed']:
        return []

    entries = [e['links']
               for e in search_results['feed']['entry']
               if 'links' in e]
    # Flatten "entries" to a simple list of links
    links = list(itertools.chain(*entries))

    urls = []
    unique_filenames = set()
    for link in links:
        if 'href' not in link:
            # Exclude links with nothing to download
            continue
        if 'inherited' in link and link['inherited'] is True:
            # Why are we excluding these links?
            continue
        if 'rel' in link and 'data#' not in link['rel']:
            # Exclude links which are not classified by CMR as "data" or "metadata"
            continue
        if 'title' in link and 'opendap' in link['title'].lower():
            # Exclude OPeNDAP links--they are responsible for many duplicates
            # This is a hack; when the metadata is updated to properly identify
            # non-datapool links, we should be able to do this in a non-hack way
            continue

        filename = link['href'].split('/')[-1]
        if filename in unique_filenames:
            # Exclude links with duplicate filenames (they would overwrite)
            continue
        unique_filenames.add(filename)

        urls.append(link['href'])

    return urls
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to search url(s)
def cmr_search(short_name, version, time_start, time_end,
               bounding_box='', polygon='', filename_filter='',
               cmr_page_size=None, cmr_url=None, cmr_file_url=None):

    cmr_file_url = (cmr_file_url.format(cmr_url, cmr_page_size))

    """Perform a scrolling CMR query for files matching input criteria."""
    """Take a look at lib_utils_generic also for setting query of orbit, dir, version, etc..."""
    cmr_query_url = build_cmr_query_url(short_name=short_name,
                                        version=version,
                                        time_start=time_start, time_end=time_end,
                                        bounding_box=bounding_box,
                                        polygon=polygon, filename_filter=filename_filter,
                                        cmr_file_url=cmr_file_url)

   # logging.info(' ----> Querying for data:\n\t{0}\n'.format(cmr_query_url))

    logging.info(' ----> Querying for data:\n' + str(cmr_query_url))

    cmr_scroll_id = None
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    try:
        urls = []
        while True:
            req = Request(cmr_query_url)
            if cmr_scroll_id:
                req.add_header('cmr-scroll-id', cmr_scroll_id)

            try:
                response = urlopen(req, context=ctx)
            except HTTPError as e:
                logging.warning(' ===> HTTP error {0}, {1}'.format(e.code, e.reason))
                if e.code == 429:
                    logging.warning(' ===> No data found for the specified time period.')
                    return []
                else:
                    raise RuntimeError(' ===> HTTP error {0}, {1}'.format(e.code, e.reason))

            if not cmr_scroll_id:
                # Python 2 and 3 have different case for the http headers
                headers = {k.lower(): v for k, v in dict(response.info()).items()}
                cmr_scroll_id = headers['cmr-scroll-id']
                hits = int(headers['cmr-hits'])
                if hits > 0:
                    logging.info(' ----> Found {0} matches.'.format(hits))
                else:
                    logging.info(' ----> Found no matches.')
            search_page = response.read()
            search_page = json.loads(search_page.decode('utf-8'))
            url_scroll_results = cmr_filter_urls(search_page)
            if not url_scroll_results:
                break
            if hits > cmr_page_size:
                print('.', end='')
                sys.stdout.flush()
            urls += url_scroll_results

        if hits > cmr_page_size:
            print()
        return urls
    except KeyboardInterrupt:
        quit()
# -------------------------------------------------------------------------------------
