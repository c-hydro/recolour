"""
Library Features:

Name:          lib_data
Author(s):     Fabio Delogu
Date:          '20260715'
Version:       '1.0.0'

Purpose:
    Validate, load and organize dynamic datasets configured in JSON.

Supported types:
    - netcdf, tiff
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import gc
import logging
import os
import numpy as np
import pandas as pd

from copy import deepcopy
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

from lib_utils_geo import check_grids
from lib_utils_datasets import validate_datasets, read_datasets
from lib_metrics import Metrics

from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)

# constants
SUPPORTED_DATASET_TYPES = {"netcdf_hmc_state", "tiff",}

SEASON_START_MONTH_TAG = {"ALL": "ALL","DJF": "12","MAM": "03","JJA": "06","SON": "09",}
SEASON_END_MONTH_TAG = {"ALL": "ALL","DJF": "02","MAM": "05","JJA": "08","SON": "11",}

DEFAULT_TAG = "ALL"
DEFAULT_AGGREGATION_TYPE = "season" # season, month
DEFAULT_AGGREGATION_BY = "name" # name, start_month, end_month (only for season type)
DEFAULT_FREQUENCY = "D"
DEFAULT_REFERENCE_GROUP = "reference"
DEFAULT_OTHER_GROUP = "other"
DEFAULT_MISSING_THRESHOLD = 90.0
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# class to organize dynamic datasets
class DynamicDatasets:

    def __init__(
            self,
            datasets_cfg: Mapping[str, Any],
            geo: Mapping[str, Any],
            time_period_common: Any, time_period_reference: Any, time_period_other: Any,
            time_frequency: str = DEFAULT_FREQUENCY,
            time_tag: str = DEFAULT_TAG,
            aggregation_type: str = DEFAULT_AGGREGATION_TYPE, aggregation_by: str = DEFAULT_AGGREGATION_BY,
            reference_group: str = DEFAULT_REFERENCE_GROUP,
            other_group: str = DEFAULT_OTHER_GROUP,
            check_grids: bool = True,
            check_grids_once: bool = True,
            raise_error: bool = True,
            skip_missing: bool = True,
            missing_threshold: float = DEFAULT_MISSING_THRESHOLD,
    ):

        # check datasets configuration
        if not isinstance(datasets_cfg, Mapping):
            raise TypeError("The datasets configuration must be a dictionary.")
        if not isinstance(geo, Mapping):
            raise TypeError("The geo object must be a dictionary.")

        self.datasets_cfg_raw = deepcopy(dict(datasets_cfg))
        self.geo = geo

        self.time_frequency = time_frequency
        self.time_period_common = time_period_common
        self.time_period_reference = time_period_reference
        self.time_period_other = time_period_other

        if len(self.time_period_common) == 0: raise ValueError("The time period common is empty.")
        if len(self.time_period_reference) == 0: raise ValueError("The time period reference is empty.")
        if len(self.time_period_other) == 0: raise ValueError("The time period other is empty.")

        self.time_start = self._parse_time(time_value=time_period_common.iloc[0], time_name="time_start",)
        self.time_end = self._parse_time(time_value=time_period_common.iloc[-1],time_name="time_end",)
        self.aggregation_type = str(aggregation_type).strip().lower()
        self.aggregation_by = str(aggregation_by).strip().lower()
        self.time_tag = self._define_time_tag(
            time_value=self.time_end, type_tag=self.aggregation_type,
            requested_tag=time_tag, season_tag_type=self.aggregation_by)

        if self.time_start > self.time_end:
            raise ValueError(
                f"time_start '{self.time_start}' is later than "
                f"time_end '{self.time_end}'."
            )

        self.reference_group = reference_group
        self.other_group = other_group

        self.check_grids = check_grids
        self.check_grids_once = check_grids_once
        self.raise_error = raise_error
        self.skip_missing = skip_missing

        # Validate missing-data threshold.
        try:
            missing_threshold = float(missing_threshold)
        except (TypeError, ValueError) as exc:
            raise TypeError("'missing_threshold' must be numeric.") from exc

        if not 0.0 <= missing_threshold <= 100.0:
            raise ValueError("'missing_threshold' must be between 0 and 100.")

        self.missing_threshold = missing_threshold

        self.datasets_cfg: Dict[str, Dict[str, Any]] = {}

        # The grid compatibility check can normally be performed only once.
        self._grids_checked = False

        # Availability statistics are reset when iterate() starts.
        self._availability_stats: Dict[str, Any] = {
            "time_steps_total": len(self.time_period_common),
            "time_steps_available": 0,
            "time_steps_missing": 0,
            "missing_percentage": 0.0,
            "missing_threshold": self.missing_threshold,
            "missing_files": [],
        }
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public method to organize configurations
    def organize(self) -> Dict[str, Dict[str, Any]]:

        # info method start
        logger.info(" ----> Organize dynamic datasets ...")

        # validata datasets
        self.datasets_cfg = self._validate_datasets()

        # info dynamic datasets
        logger.info(f" -----> Time start: {self.time_start}")
        logger.info(f" -----> Time end: {self.time_end}")
        logger.info(f" -----> Time frequency: {self.time_frequency}")
        logger.info(f" -----> Number of time steps: {len(self.time_period_common)}")

        # info method end
        logger.info(" ----> Organize dynamic datasets ... DONE")

        return self.datasets_cfg

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public generator to load one pair of maps at a time
    def iterate(self) -> Iterator[Dict[str, Any]]:

        # Define datasets if needed.
        if not self.datasets_cfg:
            self.organize()

        logger.info(" ----> Iterate dynamic datasets ...")

        time_steps_total = len(self.time_period_common)

        # Reset availability statistics for this iteration.
        self._availability_stats = {
            "time_steps_total": time_steps_total,
            "time_steps_available": 0,
            "time_steps_missing": 0,
            "missing_percentage": 0.0,
            "missing_threshold": self.missing_threshold,
            "missing_files": [],
        }

        # iterate
        for time_id, (time_step_common, time_step_reference, time_step_other) in enumerate(
                zip(self.time_period_common, self.time_period_reference,self.time_period_other,), start=1,):

            # info start time
            logger.info(
                f" -----> Time step {time_id}/{time_steps_total}: "
                f"common={time_step_common:%Y-%m-%d %H:%M}, "
                f"reference={time_step_reference:%Y-%m-%d %H:%M}, "
                f"other={time_step_other:%Y-%m-%d %H:%M} ... "
            )

            # Resolve dataset paths without requiring files to exist.
            # File availability is managed explicitly below.
            datasets_step_cfg = self._validate_datasets(
                time_common=time_step_common,
                time_reference=time_step_reference, time_other=time_step_other
            )

            datasets_obj = None
            try:

                files_available, missing_files = self._check_time_files(
                    resolved_cfg=datasets_step_cfg,
                )

                if not files_available:

                    self._availability_stats["time_steps_missing"] += 1

                    self._availability_stats["missing_files"].append({
                        "time": time_step_common,
                        "files": list(missing_files),
                    })

                    missing_message = (f"Missing datasets at time '{time_step_common}': "
                                       + ", ".join(missing_files))

                    if self.skip_missing:
                        logger.warning(f"{missing_message}. Time step skipped.")
                        continue

                    raise FileNotFoundError(missing_message)

                self._availability_stats["time_steps_available"] += 1

                # Load the reference and candidate maps only when all
                # files required for the current time step are available.
                datasets_obj = self._get_datasets(resolved_cfg=datasets_step_cfg,time_step=time_step_common,)

                # Grid geometry normally does not change with time.
                if self._must_check_grids():

                    is_compatible = self._check_datasets(datasets_obj=datasets_obj,)
                    if is_compatible:
                        self._grids_checked = True

                yield {
                    "time": time_step_common,
                    "paths": {group_name: group_cfg["file_path"]for group_name, group_cfg in datasets_step_cfg.items()},
                    **datasets_obj,
                }

            finally:

                if datasets_obj is not None:
                    datasets_obj.clear()

                datasets_obj = None
                gc.collect()

                # info end time
                logger.info(
                    f" -----> Time step {time_id}/{time_steps_total}: "
                    f"common={time_step_common:%Y-%m-%d %H:%M}, "
                    f"reference={time_step_reference:%Y-%m-%d %H:%M}, "
                    f"other={time_step_other:%Y-%m-%d %H:%M} ... DONE"
                )

        # --------------------------------------------------------------------------
        # Evaluate the missing-data percentage after checking all requested times.
        time_steps_missing = self._availability_stats["time_steps_missing"]

        if time_steps_total > 0:
            missing_percentage = (100.0 * time_steps_missing/ time_steps_total)
        else:
            missing_percentage = 0.0

        self._availability_stats["missing_percentage"] = missing_percentage

        logger.info(
            " -----> Dataset availability: "
            "requested=%d, available=%d, missing=%d, "
            "missing_percentage=%.2f%%, threshold=%.2f%%",
            time_steps_total,
            self._availability_stats["time_steps_available"],
            time_steps_missing,
            missing_percentage,
            self.missing_threshold,
        )

        # Exit only when the missing percentage is greater than the threshold.
        if missing_percentage > self.missing_threshold:
            raise RuntimeError(
                "Missing-data threshold exceeded: "
                f"{time_steps_missing}/{time_steps_total} time steps "
                f"are missing ({missing_percentage:.2f}%). "
                f"The configured threshold is "
                f"{self.missing_threshold:.2f}%."
            )

        if time_steps_missing > 0:
            logger.warning(
                " -----> Missing datasets are within the accepted threshold: "
                "%.2f%% <= %.2f%%. Analysis continues.",
                missing_percentage,
                self.missing_threshold,
            )

        logger.info(" ----> Iterate dynamic datasets ... DONE")

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public method to analyze metrics
    def analyze_metrics(self, metrics_cfg: Optional[Dict[str, Any]] = None,) -> Dict[str, Any]:

        # organize datasets if needed
        if not self.datasets_cfg:
            self.organize()

        logger.info(" ----> Analyze dynamic datasets ...")

        # normalize metrics configuration
        metrics_cfg = ({} if metrics_cfg is None else dict(metrics_cfg))

        # get metrics options
        min_observations = metrics_cfg.get("min_observations", metrics_cfg.get("min_observations", 10),)
        metrics_dtype = np.dtype(metrics_cfg.get("dtype", "float64",))
        continue_on_error = bool(metrics_cfg.get("continue_on_error", False,))

        # validate options
        if not isinstance(min_observations, (int, np.integer)):
            raise TypeError("Metrics option 'min_observations' must be an integer.")
        if min_observations < 1:
            raise ValueError("Metrics option 'min_observations' must be greater than zero.")

        logger.info(
            " -----> Metrics configuration: "
            "min_observations=%d, dtype=%s, continue_on_error=%s",
            min_observations,
            metrics_dtype,
            continue_on_error,
        )

        # initialize metrics accumulator
        metrics_obj = Metrics(
            grid_shape=None,
            min_observations=min_observations,
            dtype=metrics_dtype,
        )

        # initialize counters
        time_steps_total = len(self.time_period_common)
        time_steps_processed = 0
        time_steps_skipped = 0
        valid_pairs_total = 0

        first_time_processed = None
        last_time_processed = None
        reference_grid = None

        # iterate over dataset pairs
        for time_index, datasets_step in enumerate(self.iterate(), start=1,):

            time_step = datasets_step["time"]

            logger.info(" -----> Processing time step %d/%d: %s",
                        time_index, time_steps_total, time_step,)

            try:

                # get datasets
                reference_data = datasets_step[self.reference_group]
                other_data = datasets_step[self.other_group]

                # get arrays
                reference_values = np.asarray(reference_data["values"])
                other_values = np.asarray(other_data["values"])

                # store grid information once
                if reference_grid is None:
                    reference_grid = {
                        "longitude": reference_data.get("longitude"),
                        "latitude": reference_data.get("latitude"),
                        "transform": reference_data.get("transform"),
                        "crs": reference_data.get("crs"),
                        "width": reference_data.get("width",reference_values.shape[-1],),
                        "height": reference_data.get("height", reference_values.shape[-2],),
                        "shape": reference_values.shape,
                        "metadata": reference_data.get("metadata", {},),
                    }

                # update metrics
                valid_pairs_step = metrics_obj.update(
                    reference_values=reference_values,
                    dataset_values=other_values,
                )

                # the update is counted as processed only when
                # at least one valid pair is available
                if valid_pairs_step > 0:

                    valid_pairs_total += valid_pairs_step
                    time_steps_processed += 1

                    if first_time_processed is None:
                        first_time_processed = time_step

                    last_time_processed = time_step

                    logger.info(" -----> Valid paired pixels: %d",valid_pairs_step,)

                else:

                    time_steps_skipped += 1

                    logger.warning(
                        " -----> Time step skipped because no valid "
                        "paired pixels are available: %s",
                        time_step,
                    )

            except Exception as exc:

                time_steps_skipped += 1
                logger.error(" -----> Time step '%s' failed: %s", time_step,exc,)

                if not continue_on_error: raise

        # make sure at least one valid update was completed
        if metrics_obj.update_count == 0:
            raise RuntimeError("Metrics analysis did not process any valid dataset pair.")

        # compute grid metrics maps
        metrics_grid = metrics_obj.finalize_grid()
        # compute scalar metric maps
        metrics_scalar = metrics_obj.finalize_scalar()

        # get accumulator information
        metrics_info = metrics_obj.info()

        # organize results
        analysis_data = {
            "metrics_grid": metrics_grid, "metrics_scalar": metrics_scalar,
            "grid": reference_grid,
            "time": {
                "start_requested": self.time_period_common.iloc[0],
                "end_requested": self.time_period_common.iloc[-1],
                "first_processed": first_time_processed,
                "last_processed": last_time_processed,
                "frequency": self.time_frequency,
            },
            "datasets": {
                "reference_group": self.reference_group,
                "other_group": self.other_group,
            },
            "processing": {
                "time_steps_total": time_steps_total,
                "time_steps_processed": time_steps_processed,
                "time_steps_skipped": time_steps_skipped,
                "valid_pairs_total": valid_pairs_total,
                "availability": deepcopy(self._availability_stats),
                **metrics_info,
            },
            "configuration": {
                "min_observations": min_observations,
                "dtype": str(metrics_dtype),
                "continue_on_error": continue_on_error,
                "missing_threshold": self.missing_threshold,
            },
        }

        logger.info(
            " -----> Metrics analysis completed: "
            "requested=%d, processed=%d, skipped=%d, valid_pairs=%d",
            time_steps_total,
            time_steps_processed,
            time_steps_skipped,
            valid_pairs_total,
        )

        logger.info(
            " -----> Metrics grid: shape=%s, cells=%d, maximum_observations=%d",
            metrics_info["grid_shape"],
            metrics_info["cells_with_observations"],
            metrics_info["maximum_observations"],
        )

        logger.info(" ----> Analyze dynamic datasets ... DONE")

        return analysis_data

    # ----------------------------------------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------------------------------
    # public method to analyze weights from computed metric maps
    def analyze_weights(self,
                        analysis_data: Dict[str, Any],
                        weights_cfg: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Compute spatial nudging weights using the 2D metric maps produced
        by the ``analyze`` method.

        The final weight is computed as:

            weight =
                maximum_weight
                * correlation_weight
                * error_weight
                * observation_weight

        where:

            correlation_weight
                represents agreement between reference and other datasets;

            error_weight
                penalizes cells with high normalized ubRMSD;

            observation_weight
                penalizes cells with a limited number of valid observations.

        Parameters
        ----------
        analysis_data : dict
            Output dictionary returned by ``analyze()``.

        weights_cfg : dict, optional
            Weight configuration.

        Returns
        -------
        dict
            Updated analysis dictionary containing the ``weights`` section.
        """

        logger.info(" ----> Analyze nudging weights ...")

        # check analysis data
        if not isinstance(analysis_data, dict):
            raise TypeError(
                "'analysis_data' must be a dictionary."
            )

        # get metric maps
        metrics_data = analysis_data.get("metrics_grid")

        if metrics_data is None:
            raise KeyError(
                "Analysis data do not contain the 'metrics' section."
            )

        if not isinstance(metrics_data, dict):
            raise TypeError(
                "Analysis section 'metrics' must be a dictionary."
            )

        # normalize weight configuration
        weights_cfg = (
            {}
            if weights_cfg is None
            else dict(weights_cfg)
        )

        # get general configuration
        method = str(
            weights_cfg.get(
                "method",
                "correlation_error",
            )
        ).strip().lower()

        maximum_weight = float(
            weights_cfg.get(
                "maximum_weight",
                1.0,
            )
        )

        if not 0.0 <= maximum_weight <= 1.0:
            raise ValueError(
                "Weight option 'maximum_weight' must be between zero and one."
            )

        # correlation configuration
        correlation_cfg = weights_cfg.get(
            "correlation",
            {},
        )

        if correlation_cfg is None:
            correlation_cfg = {}

        if not isinstance(correlation_cfg, dict):
            raise TypeError(
                "Weight option 'correlation' must be a dictionary."
            )

        minimum_correlation = float(
            correlation_cfg.get(
                "minimum",
                0.3,
            )
        )

        correlation_power = float(
            correlation_cfg.get(
                "power",
                1.0,
            )
        )

        if minimum_correlation >= 1.0:
            raise ValueError(
                "Weight option 'correlation.minimum' must be lower than one."
            )

        if correlation_power <= 0.0:
            raise ValueError(
                "Weight option 'correlation.power' must be greater than zero."
            )

        # error configuration
        error_cfg = weights_cfg.get(
            "error",
            {},
        )

        if error_cfg is None:
            error_cfg = {}

        if not isinstance(error_cfg, dict):
            raise TypeError(
                "Weight option 'error' must be a dictionary."
            )

        error_power = float(
            error_cfg.get(
                "power",
                2.0,
            )
        )

        minimum_std = float(
            error_cfg.get(
                "minimum_std",
                0.01,
            )
        )

        if error_power <= 0.0:
            raise ValueError(
                "Weight option 'error.power' must be greater than zero."
            )

        if minimum_std <= 0.0:
            raise ValueError(
                "Weight option 'error.minimum_std' must be greater than zero."
            )

        # observation configuration
        observations_cfg = weights_cfg.get(
            "observations",
            {},
        )

        if observations_cfg is None:
            observations_cfg = {}

        if not isinstance(observations_cfg, dict):
            raise TypeError(
                "Weight option 'observations' must be a dictionary."
            )

        minimum_observations = int(
            observations_cfg.get(
                "minimum",
                analysis_data.get(
                    "configuration",
                    {},
                ).get(
                    "min_observations",
                    10,
                ),
            )
        )

        target_observations = int(
            observations_cfg.get(
                "target",
                minimum_observations,
            )
        )

        if minimum_observations < 1:
            raise ValueError(
                "Weight option 'observations.minimum' must be greater than zero."
            )

        if target_observations < minimum_observations:
            raise ValueError(
                "Weight option 'observations.target' must be greater than or "
                "equal to 'observations.minimum'."
            )

        # check supported method
        supported_methods = [
            "correlation_error",
        ]

        if method not in supported_methods:
            raise NotImplementedError(
                f"Weight method '{method}' is not supported. "
                f"Available methods: {supported_methods}."
            )

        # check required metric maps
        required_metrics = [
            "n_obs",
            "correlation",
            "ubrmsd",
            "std_reference",
        ]

        for metric_name in required_metrics:

            if metric_name not in metrics_data:
                raise KeyError(
                    f"Metric '{metric_name}' is required to compute weights."
                )

        # convert metric maps
        values_n_obs = np.asarray(
            metrics_data["n_obs"],
            dtype=np.float64,
        )

        values_correlation = np.asarray(
            metrics_data["correlation"],
            dtype=np.float64,
        )

        values_ubrmsd = np.asarray(
            metrics_data["ubrmsd"],
            dtype=np.float64,
        )

        values_std_reference = np.asarray(
            metrics_data["std_reference"],
            dtype=np.float64,
        )

        # check metric dimensions
        if values_n_obs.ndim != 2:
            raise ValueError(
                f"Metric 'n_obs' must be 2D, but its shape is "
                f"{values_n_obs.shape}."
            )

        grid_shape = values_n_obs.shape

        metric_maps = {
            "correlation": values_correlation,
            "ubrmsd": values_ubrmsd,
            "std_reference": values_std_reference,
        }

        for metric_name, metric_values in metric_maps.items():

            if metric_values.shape != grid_shape:
                raise ValueError(
                    f"Metric '{metric_name}' shape {metric_values.shape} "
                    f"differs from expected grid shape {grid_shape}."
                )

        # initialize output arrays
        weight_correlation = np.zeros(
            grid_shape,
            dtype=np.float64,
        )

        weight_error = np.zeros(
            grid_shape,
            dtype=np.float64,
        )

        weight_observations = np.zeros(
            grid_shape,
            dtype=np.float64,
        )

        normalized_error = np.full(
            grid_shape,
            np.nan,
            dtype=np.float64,
        )

        # define valid metric cells
        valid_mask = (
                np.isfinite(values_n_obs)
                & np.isfinite(values_correlation)
                & np.isfinite(values_ubrmsd)
                & np.isfinite(values_std_reference)
                & (values_n_obs >= minimum_observations)
                & (values_std_reference >= minimum_std)
        )

        # --------------------------------------------------------------------------
        # compute correlation weight
        #
        # correlation <= minimum -> zero
        # correlation == one     -> one

        correlation_range = (
                1.0 - minimum_correlation
        )

        correlation_score = (
                                    values_correlation
                                    - minimum_correlation
                            ) / correlation_range

        correlation_score = np.clip(
            correlation_score,
            0.0,
            1.0,
        )

        weight_correlation[valid_mask] = (
                correlation_score[valid_mask]
                ** correlation_power
        )

        # --------------------------------------------------------------------------
        # compute normalized error
        #
        # normalized_error = ubRMSD / reference standard deviation

        normalized_error[valid_mask] = (
                values_ubrmsd[valid_mask]
                / values_std_reference[valid_mask]
        )

        # error weight:
        #
        #     1 / (1 + normalized_error ** error_power)

        weight_error[valid_mask] = (
                1.0
                / (
                        1.0
                        + normalized_error[valid_mask]
                        ** error_power
                )
        )

        # --------------------------------------------------------------------------
        # compute observation weight

        if target_observations == minimum_observations:

            observation_score = np.where(
                values_n_obs >= minimum_observations,
                1.0,
                0.0,
            )

        else:

            observation_score = (
                                        values_n_obs
                                        - minimum_observations
                                ) / (
                                        target_observations
                                        - minimum_observations
                                )

            observation_score = np.clip(
                observation_score,
                0.0,
                1.0,
            )

        weight_observations[valid_mask] = (
            observation_score[valid_mask]
        )

        # --------------------------------------------------------------------------
        # combine weight components

        values_weight = (
                maximum_weight
                * weight_correlation
                * weight_error
                * weight_observations
        )

        # remove invalid cells
        values_weight[~valid_mask] = 0.0

        # protect final range
        values_weight = np.clip(
            values_weight,
            0.0,
            maximum_weight,
        )

        # compute summary information
        positive_weight_mask = (
                np.isfinite(values_weight)
                & (values_weight > 0.0)
        )

        positive_weights = values_weight[
            positive_weight_mask
        ]

        valid_cells = int(
            np.count_nonzero(valid_mask)
        )

        positive_cells = int(
            np.count_nonzero(positive_weight_mask)
        )

        total_cells = int(
            values_weight.size
        )

        if positive_cells > 0:

            weight_min = float(
                np.min(positive_weights)
            )

            weight_mean = float(
                np.mean(positive_weights)
            )

            weight_median = float(
                np.median(positive_weights)
            )

            weight_max = float(
                np.max(positive_weights)
            )

            weight_std = float(
                np.std(positive_weights)
            )

        else:

            weight_min = np.nan
            weight_mean = np.nan
            weight_median = np.nan
            weight_max = np.nan
            weight_std = np.nan

        # organize weight results
        weights_data = {
            "weight": values_weight,
            "weight_correlation": weight_correlation,
            "weight_error": weight_error,
            "weight_observations": weight_observations,
            "normalized_error": normalized_error,
            "valid_mask": valid_mask,
            "summary": {
                "total_cells": total_cells,
                "valid_cells": valid_cells,
                "positive_cells": positive_cells,
                "zero_weight_cells": total_cells - positive_cells,
                "weight_min": weight_min,
                "weight_mean": weight_mean,
                "weight_median": weight_median,
                "weight_max": weight_max,
                "weight_std": weight_std,
            },
            "configuration": {
                "method": method,
                "maximum_weight": maximum_weight,
                "correlation": {
                    "minimum": minimum_correlation,
                    "power": correlation_power,
                },
                "error": {
                    "power": error_power,
                    "minimum_std": minimum_std,
                },
                "observations": {
                    "minimum": minimum_observations,
                    "target": target_observations,
                },
            },
        }

        # store weights in analysis results
        analysis_data["weights"] = weights_data

        logger.info(
            " -----> Nudging weights: "
            "valid=%d/%d, positive=%d, min=%.4f, mean=%.4f, max=%.4f",
            valid_cells,
            total_cells,
            positive_cells,
            weight_min,
            weight_mean,
            weight_max,
        )

        logger.info(" ----> Analyze nudging weights ... DONE")

        return analysis_data

    # ----------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to summarize analysis results
    @staticmethod
    def summarize(
            analysis_data: Dict[str, Any],
    ) -> Dict[str, Any]:

        # check analysis data
        if not isinstance(analysis_data, dict):
            raise TypeError(
                "'analysis_data' must be a dictionary."
            )

        # --------------------------------------------------------------------------
        # get metric maps
        if "metrics_grid" not in analysis_data:
            raise KeyError(
                "Analysis data do not contain the 'metrics' section."
            )

        metrics_data = analysis_data["metrics_grid"]

        metrics_summary = {
            "n_obs": metrics_data["n_obs"],
            "mean_reference": metrics_data["mean_reference"],
            "mean_dataset": metrics_data["mean_dataset"],
            "bias": metrics_data["bias"],
            "rmse": metrics_data["rmse"],
            "ubrmsd": metrics_data["ubrmsd"],
            "correlation": metrics_data["correlation"],
            "std_reference": metrics_data["std_reference"],
            "std_dataset": metrics_data["std_dataset"],
        }

        # --------------------------------------------------------------------------
        # get weight maps
        weights_data = analysis_data.get(
            "weights",
            None,
        )

        if weights_data is not None:

            weights_summary = {
                "weight": weights_data["weight"],
                "weight_correlation": weights_data[
                    "weight_correlation"
                ],
                "weight_error": weights_data[
                    "weight_error"
                ],
                "weight_observations": weights_data[
                    "weight_observations"
                ],
                "normalized_error": weights_data[
                    "normalized_error"
                ],
                "valid_mask": weights_data[
                    "valid_mask"
                ],
                "summary": weights_data.get(
                    "summary",
                    {},
                ),
                "configuration": weights_data.get(
                    "configuration",
                    {},
                ),
            }

        else:

            weights_summary = None

            logger.warning(
                " -----> Weight analysis is not available. "
                "Run 'analyze_weights()' before 'summarize()' "
                "to include nudging weights."
            )

        # --------------------------------------------------------------------------
        # get grid information
        if "grid" not in analysis_data:
            raise KeyError(
                "Analysis data do not contain the 'grid' section."
            )

        grid_data = analysis_data["grid"]

        grid_summary = {
            "longitude": grid_data["longitude"],
            "latitude": grid_data["latitude"],
            "transform": grid_data["transform"],
            "crs": grid_data["crs"],
            "width": grid_data.get("width"),
            "height": grid_data.get("height"),
            "shape": grid_data.get("shape"),
            "metadata": grid_data.get(
                "metadata",
                {},
            ),
        }

        # --------------------------------------------------------------------------
        # get processing information

        if "processing" not in analysis_data:
            raise KeyError(
                "Analysis data do not contain the 'processing' section."
            )

        processing_info = analysis_data["processing"]

        logger.info(
            " -----> Dynamic analysis completed: "
            "requested=%d, processed=%d, skipped=%d",
            processing_info["time_steps_total"],
            processing_info["time_steps_processed"],
            processing_info["time_steps_skipped"],
        )

        logger.info(
            " -----> Metrics grid: "
            "shape=%s, valid_cells=%d, max_observations=%d",
            processing_info["grid_shape"],
            processing_info["cells_with_observations"],
            processing_info["maximum_observations"],
        )

        # --------------------------------------------------------------------------
        # log weight information
        if weights_summary is not None:
            weight_info = weights_summary["summary"]

            logger.info(
                " -----> Nudging weights: "
                "valid_cells=%d, positive_cells=%d, "
                "zero_weight_cells=%d",
                weight_info.get(
                    "valid_cells",
                    0,
                ),
                weight_info.get(
                    "positive_cells",
                    0,
                ),
                weight_info.get(
                    "zero_weight_cells",
                    0,
                ),
            )

            logger.info(
                " -----> Nudging weight statistics: "
                "min=%.4f, mean=%.4f, median=%.4f, "
                "max=%.4f, std=%.4f",
                weight_info.get(
                    "weight_min",
                    np.nan,
                ),
                weight_info.get(
                    "weight_mean",
                    np.nan,
                ),
                weight_info.get(
                    "weight_median",
                    np.nan,
                ),
                weight_info.get(
                    "weight_max",
                    np.nan,
                ),
                weight_info.get(
                    "weight_std",
                    np.nan,
                ),
            )

        # --------------------------------------------------------------------------
        # organize summary
        analysis_summary = {
            "metrics": metrics_summary,
            "weights": weights_summary,
            "grid": grid_summary,
            "processing": processing_info,
        }

        # optionally preserve time information
        if "time" in analysis_data:
            analysis_summary["time"] = analysis_data["time"]

        # optionally preserve dataset information
        if "datasets" in analysis_data:
            analysis_summary["datasets"] = analysis_data["datasets"]

        # optionally preserve analysis configuration
        if "configuration" in analysis_data:
            analysis_summary["configuration"] = analysis_data[
                "configuration"
            ]

        return analysis_summary

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to determine whether grid checking is required
    def _must_check_grids(self) -> bool:
        if not self.check_grids:
            return False
        if not self.check_grids_once:
            return True
        return not self._grids_checked
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to format one time template
    @staticmethod
    def _format_time_template(template: str, time_step: pd.Timestamp, field_name: str,) -> str:
        """
        Resolve templates such as:

        /data/hmc/{time:%Y/%m/%d}
        hmc.state-grid.{time:%Y%m%d%H%M}.nc.gz
        """

        if not isinstance(template, str):
            raise TypeError(f"Template '{field_name}' must be a string.")

        try:
            formatted_value = template.format(time=time_step.to_pydatetime(),)
        except (KeyError, ValueError, IndexError) as exc:
            raise ValueError(
                f"Unable to resolve time template '{field_name}': '{template}' using time '{time_step}'."
            ) from exc

        return formatted_value
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to check all files for the current timestamp
    @staticmethod
    def _check_time_files(resolved_cfg: Mapping[str, Mapping[str, Any]],) -> Tuple[bool, list[str]]:
        missing_files = []
        for group_name, dataset_cfg in resolved_cfg.items():
            file_path = dataset_cfg["file_path"]
            if not os.path.isfile(file_path):
                missing_files.append(f"{group_name}='{file_path}'")

        return not missing_files, missing_files
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to load maps for one timestamp
    def _get_datasets(self,
                      resolved_cfg: Mapping[str, Mapping[str, Any]], time_step: pd.Timestamp,) \
            -> Dict[str, Optional[Dict[str, Any]]]:

        # info start
        logger.info(f" ------> Get datasets for {time_step} ...")

        # get datasets geo
        datasets_geo = self.geo

        # iterate over datasets
        datasets_data = {}
        for group_name, dataset_cfg in resolved_cfg.items():

            # info datasets start
            logger.info(f" -------> Dataset group {group_name} ...")

            group_geo = {}
            if group_name in datasets_geo:
                group_geo = datasets_geo[group_name]

            # read datasets
            dataset_obj = read_datasets(dataset_cfg=dataset_cfg,
                                        group_name=group_name, dataset_key=group_name, group_geo=group_geo)
            # store datasets
            datasets_data[group_name] = dataset_obj

            # info datasets start
            logger.info(f" -------> Dataset group {group_name} ... DONE")

        # info end
        logger.info(f" ------> Get datasets for {time_step} ... DONE")

        return datasets_data

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to validate all configurations
    def _validate_datasets(self,
                           time_common: pd.Timestamp = None,
                           time_reference: pd.Timestamp = None, time_other: pd.Timestamp = None) \
            -> Dict[str, Dict[str, Any]]:

        # info start
        logger.info(" -----> Validate datasets ...")

        if not self.datasets_cfg_raw:
            raise ValueError("The datasets configuration is empty. Exit")

        if self.reference_group not in self.datasets_cfg_raw:
            raise KeyError(f"Reference dataset group '{self.reference_group}' is not configured.")

        if self.other_group not in self.datasets_cfg_raw:
            raise KeyError(f"Candidate dataset group '{self.other_group}' is not configured.")

        # iterate over datasets
        datasets_cfg_validated = {}
        for group_name, dataset_cfg in self.datasets_cfg_raw.items():

            # info datasets start
            logger.info(f" ------> Dataset group {group_name} ...")

            # define time datasets
            if time_common is not None:
                if group_name == 'reference':
                    time_datasets = time_reference
                elif group_name == 'other':
                    time_datasets = time_other
                else:
                    raise ValueError(f"Dataset group '{group_name}' is not allowed.")
            else:
                time_datasets = None

            if not isinstance(dataset_cfg, Mapping):
                raise TypeError(f"Dataset group '{group_name}' must be a dictionary.")
            if not dataset_cfg:
                raise ValueError(f"Dataset group '{group_name}' is empty.")

            # check time for static or dynamic datasets
            if time_common is None:

                dataset_validated = validate_datasets(
                    dataset_cfg=dataset_cfg,
                    group_name=group_name,
                    dataset_key=group_name,
                    time_step_common=None, time_step_datasets=None,
                    check_file=False,
                )

            else:

                # Resolve the dynamic path, but do not raise if the file is absent.
                # Missing files are counted and evaluated in iterate().
                dataset_validated = validate_datasets(
                    dataset_cfg=dataset_cfg,
                    group_name=group_name,
                    dataset_key=group_name,
                    time_step_common=time_common, time_step_datasets=time_datasets,
                    check_file=False,
                )

            # check dataset type
            dataset_type = dataset_validated["type"]
            if dataset_type not in SUPPORTED_DATASET_TYPES:
                raise ValueError(
                    f"Dataset type '{dataset_type}' configured for "
                    f"'{group_name}' is not supported. "
                    f"Supported types are: "
                    f"{sorted(SUPPORTED_DATASET_TYPES)}."
                )

            # During validation, folder and filename are still templates.
            # File existence must not be checked here.
            datasets_cfg_validated[group_name] = dataset_validated

            # info datasets end
            logger.info(f" ------> Dataset group {group_name} ... DONE")

        # info end
        logger.info(" -----> Validate datasets ... DONE")

        return datasets_cfg_validated

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to compare grids
    def _check_datasets(self,datasets_obj: Mapping[str, Any],) -> bool:

        # info start
        logger.info(" ------> Check dataset grids ...")

        # get reference and other data
        reference_data = datasets_obj.get(self.reference_group, None)
        other_data = datasets_obj.get(self.other_group, None)

        # check reference and other datasets
        if reference_data is None:
            logger.warning(f"Reference dataset {self.reference_group} is missing.")
            return False
        if other_data is None:
            logger.warning(f"Candidate dataset {self.other_group} is missing.")
            return False

        # check grids compatibility
        try:
            is_compatible = check_grids(reference_data=reference_data, dataset_data=other_data,)
        except ValueError:

            # if error is mandatory
            if self.raise_error:
                raise
            # otherwise continue without grid checking
            logger.warning(" ===> Reference and candidate grids are incompatible.",exc_info=True)
            is_compatible = False

        # info end
        logger.info(" ------> Check dataset grids ... DONE")

        return is_compatible

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to parse a time value
    @staticmethod
    def _parse_time(time_value: Any,time_name: str,) -> pd.Timestamp:

        if time_value is None:
            raise ValueError(f"'{time_name}' is not defined.")

        try:
            time_obj = pd.Timestamp(time_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Unable to parse '{time_name}': '{time_value}'.") from exc

        if pd.isna(time_obj):
            raise ValueError(f"'{time_name}' is NaT.")

        return time_obj
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to define time tag
    @staticmethod
    @staticmethod
    def _define_time_tag(
            time_value: Any,
            type_tag: str,
            requested_tag: str = DEFAULT_TAG,
            season_tag_type: str = "name",
    ) -> str:
        """
        Define the tag used to identify output files.

        Parameters
        ----------
        time_value : Any
            Reference time. Used only for monthly analyses.

        type_tag : {"season", "month"}
            Type of temporal aggregation.

        requested_tag : str, default="ALL"
            Season tag (ALL, DJF, MAM, JJA, SON) when
            ``type_tag="season"``.

        season_tag_type : {"name", "start_month", "end_month"}, default="name"
            Representation used for seasonal analyses:

            - "name": season name (e.g. DJF)
            - "start_month": first month of the season (e.g. DJF -> "12")
            - "end_month": last month of the season (e.g. DJF -> "02")

        Returns
        -------
        str
            Tag identifying the requested period.
        """

        season_start_month = {
            "ALL": "ALL_START",
            "DJF": "12",
            "MAM": "03",
            "JJA": "06",
            "SON": "09",
        }

        season_end_month = {
            "ALL": "ALL_END",
            "DJF": "02",
            "MAM": "05",
            "JJA": "08",
            "SON": "11",
        }

        type_tag = str(type_tag).strip().lower()
        requested_tag = str(requested_tag).strip().upper()
        season_tag_type = str(season_tag_type).strip().lower()
        time_value = pd.Timestamp(time_value)

        if type_tag == "season":

            if requested_tag not in season_start_month:
                raise ValueError(
                    f'Unsupported season tag "{requested_tag}". '
                    f"Supported values are: {list(season_start_month.keys())}"
                )

            if season_tag_type == "name":
                return requested_tag

            if season_tag_type == "start_month":
                return season_start_month[requested_tag]

            if season_tag_type == "end_month":
                return season_end_month[requested_tag]

            raise ValueError(
                f'Unsupported season_tag_type "{season_tag_type}". '
                'Supported values are: "name", "start_month", "end_month".'
            )

        if type_tag == "month":
            return time_value.strftime("%m")

        raise ValueError(
            f'Unsupported type_tag "{type_tag}". '
            'Supported values are: "season", "month".'
        )
    # ------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
