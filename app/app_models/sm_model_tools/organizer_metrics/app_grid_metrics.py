#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RECOLOUR APPS - GRID METRICS - REprocess paCkage for sOiL mOistUre pRoducts

__date__ = '20260715'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_grid_metrics.py -settings_file configuration.json -time_start "YYYY-MM-DD HH:MM" -time_end "YYYY-MM-DD HH:MM"

Version(s):
20260715 (1.0.0) --> First development
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import time
import sys
import argparse
import numpy as np

from lib_utils_logging import get_logger
from lib_utils_io import read_file_json
from lib_utils_time import create_time_period, create_seasons_period, remove_seasons_empty

from lib_geo import GeoDatasets
from lib_data import DynamicDatasets
from lib_results import Results

from config_info import LOGGER_NAME, ALG_NAME, ALG_RELEASE, ALG_VERSION

# set logger
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()
    # get settings
    settings = read_file_json(args.settings_file)

    # get configuration(s)
    log_cfg = settings.get("log", {})
    time_cfg = settings['time']
    reference_cfg = settings["datasets"]["reference"]
    other_cfg = settings["datasets"]["other"]
    metrics_cfg = settings.get("metrics", {})
    weights_cfg = settings.get("weights", {})
    img_cfg = settings.get("img", {})
    results_cfg = settings.get("results", {})

    # get reference time period
    try:
        time_period_root, time_period_ref, time_period_other, time_info = create_time_period(
            time_cfg=time_cfg, time_start_cli=args.time_start, time_end_cli=args.time_end)
    except Exception as exc:
        print(f" ===> ERROR: parsing time: {exc}")
        sys.exit(1)

    # get seasons time periods
    seasons_time_period = create_seasons_period(
        time_period_root, time_period_ref, time_period_other, seasons=metrics_cfg.get('seasons', 'ALL'))

    # get logger
    get_logger(logger, log_cfg, reference_time=time_period_root[-1])
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start message
    logger.info(' ============================================================================ ')
    logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger.info(' ==> START ... ')
    logger.info(' ')

    logger.info(f" ---> Settings file:      {args.settings_file}")

    logger.info(f" ---> Time start:         {time_info['common']['start']}")
    logger.info(f" ---> Time end:           {time_info['common']['end']}")
    logger.info(f" ---> Time steps:         {time_info['common']['steps']}")
    logger.info(f" ---> Time frequency:     {time_info['frequency']}")
    logger.info(f" ---> Seasons:            {metrics_cfg.get('seasons', 'ALL')}")

    # datasets
    logger.info(" ---> Reference:         %s (%s)",
                reference_cfg["name"], reference_cfg["type"],)
    logger.info(" ---> Other:             %s (%s)",
                other_cfg["name"],other_cfg["type"],)

    # metrics
    logger.info(
        " ---> Metrics:           min_obs=%s dtype=%s missing_thr=%.1f%%",
        metrics_cfg.get("min_observations"),
        metrics_cfg.get("dtype"),
        metrics_cfg.get("missing_threshold", 90.0),
    )

    # weights
    logger.info(" ---> Weights:           %s (max=%.2f)",
                weights_cfg.get("method"), weights_cfg.get("maximum_weight"),)

    # outputs
    logger.info(" ---> Images:            %s",img_cfg.get("folder"),)
    logger.info(" ---> Results:           %s",results_cfg.get("folder"),)

    start_time = time.time()

    # remove empty time periods
    seasons_time_period = remove_seasons_empty(seasons_time_period)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # drive geo datasets
    driver_geo = GeoDatasets(geo_cfg=settings['geo'])
    # organize geo datasets
    geo_datasets = driver_geo.organize()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # iterate over seasons
    for season_tag, season_period in seasons_time_period.items():

        # get periods
        common_period = season_period["common"]
        reference_period, other_period = season_period['reference'], season_period['other']

        # info season start
        logger.info(
            f" ----> Season {season_tag:<4}: "
            f"{common_period.iloc[0]:%Y-%m-%d %H:%M} --> "
            f"{common_period.iloc[-1]:%Y-%m-%d %H:%M} "
            f"({len(common_period)} steps) ..."
        )

        # drive dynamic datasets
        driver_data = DynamicDatasets(
            datasets_cfg=settings["datasets"],
            geo=geo_datasets,
            time_tag=season_tag,
            aggregation_type=metrics_cfg.get("aggregation_type", "season"),
            aggregation_by=metrics_cfg.get("aggregation_by", "name"),
            time_period_common=common_period,
            time_period_reference=reference_period, time_period_other=other_period,
            time_frequency=time_info['frequency'],
            reference_group="reference", other_group="other",
            check_grids=True, raise_error=True,
            missing_threshold=metrics_cfg.get("threshold_percentage", 90.0)
        )
        # get time aggregation
        time_tag = driver_data.time_tag

        # organize dynamic datasets
        dynamic_datasets = driver_data.organize()
        # analyze dynamic datasets
        dynamic_analysis = driver_data.analyze_metrics(
            metrics_cfg=settings.get("metrics",{})
        )

        # compute spatial nudging weights
        dynamic_analysis = driver_data.analyze_weights(
            analysis_data=dynamic_analysis,
            weights_cfg=settings.get("weights",{})
        )

        # summarize dynamic datasets
        dynamic_summary = driver_data.summarize(dynamic_analysis)

        # initialize output driver
        driver_results = Results(
            time_tag=time_tag,
            time_reference=time_period_root[-1],
            time_start=time_period_root[0], time_end=time_period_root[-1],
            img_cfg=settings.get("img", {}),
            results_cfg=settings.get("results", {}),
        )

        # create PNG, GeoTIFF and ASCII outputs
        dynamic_results = driver_results.organize(
            analysis_summary=dynamic_summary,
        )

        # info seasons end
        logger.info(
            f" ----> Season {season_tag:<4}: "
            f"{common_period.iloc[0]:%Y-%m-%d %H:%M} --> "
            f"{common_period.iloc[-1]:%Y-%m-%d %H:%M} "
            f"({len(common_period)} steps) ... DONE"
        )
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    alg_time_elapsed = round(time.time() - start_time, 1)

    logger.info(' ')
    logger.info(' ==> ' + ALG_NAME + ' (Version: ' + ALG_VERSION + ' Release_Date: ' + ALG_RELEASE + ')')
    logger.info(' ==> TIME ELAPSED: ' + str(alg_time_elapsed) + ' seconds')
    logger.info(' ==> ... END')
    logger.info(' ==> Bye, Bye')
    logger.info(' ============================================================================ ')
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# cli
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute grid metrics using streaming accumulators."
    )
    parser.add_argument(
        "-settings_file",
        "--settings-file",
        dest="settings_file",
        required=True,
        help="Path to the JSON configuration file.",
    )
    parser.add_argument(
        "-time_start",
        "--time-start",
        dest="time_start",
        default=None,
        help="Optional override, for example '2024-01-01 23:00'.",
    )
    parser.add_argument(
        "-time_end",
        "--time-end",
        dest="time_end",
        default=None,
        help="Optional override, for example '2025-12-31 23:00'.",
    )
    parser.add_argument(
        "--reset-checkpoint",
        action="store_true",
        help="Ignore and overwrite an existing checkpoint.",
    )
    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------
