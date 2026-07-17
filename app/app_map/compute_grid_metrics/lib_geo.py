
"""
Library Features:

Name:          lib_geo
Author(s):     Fabio Delogu
Date:          '20260715'
Version:       '1.1.0'

Purpose:
    Validate, load and organize geographical datasets configured in JSON.

Supported types:
    - ascii_grid
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

from copy import deepcopy
from typing import Any, Dict, Mapping

from lib_utils_geo import check_grids
from lib_utils_datasets import validate_datasets, read_datasets
from config_info import LOGGER_NAME, TIME_FMT_CLI

# logging
logger = logging.getLogger(LOGGER_NAME)

# constants
SUPPORTED_GEO_TYPES = {"ascii_grid",}
DEFAULT_TYPE = "ascii_grid"
DEFAULT_BAND = 1
DEFAULT_CRS = "EPSG:4326"
DEFAULT_SCALE_FACTOR = 1.0
DEFAULT_OFFSET = 0.0
DEFAULT_NODATA_VALUES = [-9999.0]
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# geographical data organizer class
class GeoDatasets:

    def __init__(
            self,
            geo_cfg: Mapping[str, Any],
            reference_group: str = "reference",
            reference_key: str = "terrain",
            check_grids: bool = True,
            raise_error: bool = True,):

        # check datasets configuration
        if not isinstance(geo_cfg, Mapping):
            raise TypeError(
                "The geographical configuration must be a dictionary."
            )

        self.geo_cfg_raw = deepcopy(dict(geo_cfg))

        self.reference_group = reference_group
        self.reference_key = reference_key

        self.check_grids = check_grids
        self.raise_error = raise_error

        self.geo_cfg, self.geo_data = {}, {}

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public method to organize geographical datasets
    def organize(self) -> Dict[str, Dict[str, Dict[str, Any]]]:

        # info method start
        logger.info(" ----> Organize geographical datasets ...")

        # validata datasets
        self.geo_cfg = self._validate_datasets()
        # load datasets
        self.geo_data = self._get_datasets()

        # check datasets
        if self.check_grids: self._check_datasets()

        # info method end
        logger.info(" ----> Organize geographical datasets ... DONE")

        return self.geo_data

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to load datasets
    def _get_datasets(self) -> Dict[str, Dict[str, Dict[str, Any]]]:

        # info method start
        logger.info(" -----> Get ... ")

        # iterate over groups
        geo_data = {}
        for group_name, group_cfg in self.geo_cfg.items():

            # info group start
            logger.info(f" ------> Group {group_name} ... ")

            # iterate over datasets
            geo_data[group_name] = {}
            for dataset_key, dataset_cfg in group_cfg.items():

                # info datasets start
                logger.info(f" -------> Datasets {dataset_key} ... ")

                # read datasets
                datasets_obj = read_datasets(dataset_cfg=dataset_cfg, group_name=group_name, dataset_key=dataset_key)
                # store datasets
                geo_data[group_name][dataset_key] = datasets_obj

                # info datasets end
                logger.info(f" -------> Datasets {dataset_key} ... DONE")

            # info group end
            logger.info(f" ------> Group {group_name} ... DONE")

        return geo_data

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to validate datasets
    def _validate_datasets(self) -> Dict[str, Dict[str, Dict[str, Any]]]:

        # info method start
        logger.info(" -----> Validate ... ")

        # check if geo configuration is available
        if not self.geo_cfg_raw:
            raise ValueError("The geographical configuration is empty. Exit")

        # iterate over group
        geo_cfg_validated = {}
        for group_name, group_cfg in self.geo_cfg_raw.items():

            # info group start
            logger.info(f" ------> Group {group_name} ... ")

            # check format of group
            if not isinstance(group_cfg, Mapping):
                raise TypeError(f"Geo group '{group_name}' must be a dictionary.")

            # insert group in the geo configuration
            geo_cfg_validated[group_name] = {}
            if not group_cfg:
                logger.warning("Geo group '%s' is empty.",group_name,)
                continue

            # iterate over datasets
            for dataset_key, dataset_cfg in group_cfg.items():

                # info datasets start
                logger.info(f" -------> Datasets {dataset_key} ... ")

                # validate datasets
                obj_datasets = validate_datasets(
                    dataset_cfg=dataset_cfg, group_name=group_name, dataset_key=dataset_key,)

                # store datasets
                geo_cfg_validated[group_name][dataset_key] = obj_datasets

                # info datasets end
                logger.info(f" -------> Datasets {dataset_key} ... DONE")

            # info group end
            logger.info(f" ------> Group {group_name} ... DONE")

        # info end
        logger.info(" -----> Validate ... DONE")

        return geo_cfg_validated

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to check all grids against the configured reference grid
    def _check_datasets(self) -> bool:

        # info start
        logger.info(" -----> Check ... ")

        # check reference group
        if self.reference_group not in self.geo_data:
            raise KeyError(f"Reference geo group '{self.reference_group}' was not loaded.")
        # get reference data
        reference_group_data = self.geo_data[self.reference_group]

        # check reference data
        if self.reference_key not in reference_group_data:
            raise KeyError(f"Reference geo dataset '{self.reference_group}.{self.reference_key}' was not loaded.")
        reference_data = reference_group_data[self.reference_key]

        # iterate over groups
        all_compatible = True
        for group_name, group_data in self.geo_data.items():

            # info group start
            logger.info(f" ------> Group {group_name} ... ")

            # iterate over datasets
            for dataset_key, dataset_data in group_data.items():

                # info datasets start
                logger.info(f" -------> Datasets {dataset_key} ... ")

                if group_name == self.reference_group and dataset_key == self.reference_key:
                    continue

                is_compatible = check_grids(
                    reference_data=reference_data,
                    dataset_data=dataset_data,
                )

                all_compatible = (
                    all_compatible
                    and is_compatible
                )

                # info datasets end
                logger.info(f" -------> Datasets {dataset_key} ... DONE")

            # info group end
            logger.info(f" ------> Group {group_name} ... DONE")

        # info end
        logger.info(" -----> Check ... DONE")

        return all_compatible

    # ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------