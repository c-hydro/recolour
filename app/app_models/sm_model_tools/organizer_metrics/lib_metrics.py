"""
Library Features:

Name:          lib_metrics
Author(s):     Fabio Delogu
Date:          '20260717'
Version:       '1.0.0'

Purpose:
    Incrementally compute pixel-wise metrics between a reference dataset
    and another dataset without storing the full temporal dataset in memory.

Computed metrics:
    - number of valid observations
    - reference mean
    - candidate mean
    - bias
    - RMSE
    - unbiased RMSD
    - Pearson correlation
    - reference standard deviation
    - candidate standard deviation
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

from typing import Any, Dict, Optional, Tuple

import numpy as np

from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to compute gridded metrics incrementally
class Metrics:

    def __init__(
            self,
            grid_shape: Optional[Tuple[int, ...]] = None,
            min_observations: int = 10,
            dtype: Any = np.float64,
    ):

        # check minimum observations
        if min_observations < 1:
            raise ValueError(
                "'min_observations' must be greater than or equal to one."
            )

        # store settings
        self.grid_shape = grid_shape
        self.min_observations = min_observations
        self.dtype = np.dtype(dtype)

        # initialize accumulator status
        self.initialized = False
        self.update_count = 0

        # initialize cumulative arrays
        self.count = None

        self.sum_reference = None
        self.sum_dataset = None

        self.sum_reference_squared = None
        self.sum_dataset_squared = None

        self.sum_cross = None

        self.sum_error = None
        self.sum_error_squared = None

        # initialize arrays if shape is already defined
        if self.grid_shape is not None:
            self._initialize(
                grid_shape=self.grid_shape,
            )

    # ------------------------------------------------------------------------------------------------------------------
    # method to initialize cumulative arrays
    def _initialize(
            self,
            grid_shape: Tuple[int, ...],
    ) -> None:

        # check shape
        if not isinstance(grid_shape, tuple):
            grid_shape = tuple(grid_shape)

        if not grid_shape:
            raise ValueError(
                "Grid shape cannot be empty."
            )

        if any(grid_size <= 0 for grid_size in grid_shape):
            raise ValueError(
                f"Grid shape contains invalid dimensions: {grid_shape}."
            )

        # store grid shape
        self.grid_shape = grid_shape

        # initialize observation counter
        self.count = np.zeros(
            self.grid_shape,
            dtype=np.uint32,
        )

        # initialize cumulative values
        self.sum_reference = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_dataset = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_reference_squared = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_dataset_squared = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_cross = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_error = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        self.sum_error_squared = np.zeros(
            self.grid_shape,
            dtype=self.dtype,
        )

        # update status
        self.initialized = True

    # ------------------------------------------------------------------------------------------------------------------
    # method to validate input arrays
    def _validate_arrays(
            self,
            reference_values: np.ndarray,
            dataset_values: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:

        # convert input data
        reference_values = np.asarray(
            reference_values,
            dtype=self.dtype,
        )

        dataset_values = np.asarray(
            dataset_values,
            dtype=self.dtype,
        )

        # check input dimensions
        if reference_values.ndim < 1:
            raise ValueError(
                "Reference values must have at least one dimension."
            )

        if dataset_values.ndim < 1:
            raise ValueError(
                "Dataset values must have at least one dimension."
            )

        # check input shapes
        if reference_values.shape != dataset_values.shape:
            raise ValueError(
                f"Reference shape {reference_values.shape} differs from "
                f"dataset shape {dataset_values.shape}."
            )

        # initialize accumulator from the first pair
        if not self.initialized:
            self._initialize(
                grid_shape=reference_values.shape,
            )

        # check accumulator shape
        if reference_values.shape != self.grid_shape:
            raise ValueError(
                f"Input shape {reference_values.shape} differs from "
                f"metrics grid shape {self.grid_shape}."
            )

        return reference_values, dataset_values

    # ------------------------------------------------------------------------------------------------------------------
    # public method to update cumulative metrics
    def update(
            self,
            reference_values: np.ndarray,
            dataset_values: np.ndarray,
    ) -> int:
        """
        Update cumulative terms using one reference-candidate grid pair.

        Returns
        -------
        int
            Number of grid cells valid in both input maps.
        """

        # validate input arrays
        reference_values, dataset_values = self._validate_arrays(
            reference_values=reference_values,
            dataset_values=dataset_values,
        )

        # valid data are cells finite in both datasets
        valid_mask = (
            np.isfinite(reference_values)
            & np.isfinite(dataset_values)
        )

        valid_count = int(
            np.count_nonzero(valid_mask)
        )

        # no valid pair is available
        if valid_count == 0:

            logger.warning(
                " ===> Metrics update skipped because no valid "
                "reference-candidate pairs are available."
            )

            return 0

        # select valid values
        reference_valid = reference_values[
            valid_mask
        ]

        dataset_valid = dataset_values[
            valid_mask
        ]

        error_valid = (
            dataset_valid - reference_valid
        )

        # update cumulative terms
        self.count[valid_mask] += 1

        self.sum_reference[valid_mask] += (
            reference_valid
        )

        self.sum_dataset[valid_mask] += (
            dataset_valid
        )

        self.sum_reference_squared[valid_mask] += (
            reference_valid * reference_valid
        )

        self.sum_dataset_squared[valid_mask] += (
            dataset_valid * dataset_valid
        )

        self.sum_cross[valid_mask] += (
            reference_valid * dataset_valid
        )

        self.sum_error[valid_mask] += (
            error_valid
        )

        self.sum_error_squared[valid_mask] += (
            error_valid * error_valid
        )

        # update number of processed grid pairs
        self.update_count += 1

        return valid_count

    # ------------------------------------------------------------------------------------------------------------------
    # public method to finalize metric maps
    def finalize_grid(self) -> Dict[str, np.ndarray]:
        """
        Compute final metric maps from the accumulated terms.
        """

        # check accumulator status
        if not self.initialized:
            raise RuntimeError(
                "Metrics accumulator was not initialized. "
                "No dataset pair was processed."
            )

        if self.update_count == 0:
            raise RuntimeError(
                "Metrics accumulator does not contain valid updates."
            )

        # create output arrays
        count_float = self.count.astype(
            self.dtype,
        )

        safe_count = np.where(
            self.count > 0,
            count_float,
            np.nan,
        )

        enough_observations = (
            self.count >= self.min_observations
        )

        # compute means
        mean_reference = (
            self.sum_reference / safe_count
        )

        mean_dataset = (
            self.sum_dataset / safe_count
        )

        # compute bias
        bias = (
            self.sum_error / safe_count
        )

        # compute mean squared error and RMSE
        mean_squared_error = (
            self.sum_error_squared / safe_count
        )

        rmse = np.sqrt(
            np.maximum(
                mean_squared_error,
                0.0,
            )
        )

        # compute variances
        variance_reference = (
            self.sum_reference_squared / safe_count
            - mean_reference ** 2
        )

        variance_dataset = (
            self.sum_dataset_squared / safe_count
            - mean_dataset ** 2
        )

        variance_reference = np.maximum(
            variance_reference,
            0.0,
        )

        variance_dataset = np.maximum(
            variance_dataset,
            0.0,
        )

        # compute standard deviations
        standard_deviation_reference = np.sqrt(
            variance_reference
        )

        standard_deviation_dataset = np.sqrt(
            variance_dataset
        )

        # compute covariance
        covariance = (
            self.sum_cross / safe_count
            - mean_reference * mean_dataset
        )

        # compute Pearson correlation
        correlation_denominator = (
            standard_deviation_reference
            * standard_deviation_dataset
        )

        correlation = np.full(
            self.grid_shape,
            np.nan,
            dtype=self.dtype,
        )

        valid_correlation = (
            enough_observations
            & np.isfinite(correlation_denominator)
            & (correlation_denominator > 0.0)
        )

        correlation[valid_correlation] = (
            covariance[valid_correlation]
            / correlation_denominator[valid_correlation]
        )

        correlation = np.clip(
            correlation,
            -1.0,
            1.0,
        )

        # compute unbiased RMSD
        error_variance = (
            mean_squared_error
            - bias ** 2
        )

        error_variance = np.maximum(
            error_variance,
            0.0,
        )

        ubrmsd = np.sqrt(
            error_variance
        )

        # organize output metrics
        metrics = {
            "n_obs": self.count.copy(),
            "mean_reference": mean_reference,
            "mean_dataset": mean_dataset,
            "bias": bias,
            "rmse": rmse,
            "ubrmsd": ubrmsd,
            "correlation": correlation,
            "std_reference": standard_deviation_reference,
            "std_dataset": standard_deviation_dataset,
        }

        # mask metrics with insufficient observations
        for metric_name, metric_values in metrics.items():

            if metric_name == "n_obs":
                continue

            metric_values = np.asarray(
                metric_values,
                dtype=self.dtype,
            )

            metric_values[
                ~enough_observations
            ] = np.nan

            metrics[metric_name] = metric_values

        return metrics
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public method to finalize domain-wide scalar metrics
    def finalize_scalar(self) -> Dict[str, Any]:
        """
        Compute scalar metrics using all valid pixel-time pairs accumulated
        over the complete domain and processing period.
        """

        if not self.initialized:
            raise RuntimeError(
                "Metrics accumulator was not initialized."
            )

        if self.update_count == 0:
            raise RuntimeError(
                "Metrics accumulator does not contain valid updates."
            )

        # total paired observations over space and time
        n_obs = int(
            np.sum(
                self.count,
                dtype=np.uint64,
            )
        )

        if n_obs == 0:
            raise RuntimeError(
                "No valid paired observations are available."
            )

        n_obs_float = float(n_obs)

        # aggregate accumulated terms
        sum_reference = float(
            np.sum(
                self.sum_reference,
                dtype=np.float64,
            )
        )

        sum_dataset = float(
            np.sum(
                self.sum_dataset,
                dtype=np.float64,
            )
        )

        sum_reference_squared = float(
            np.sum(
                self.sum_reference_squared,
                dtype=np.float64,
            )
        )

        sum_dataset_squared = float(
            np.sum(
                self.sum_dataset_squared,
                dtype=np.float64,
            )
        )

        sum_cross = float(
            np.sum(
                self.sum_cross,
                dtype=np.float64,
            )
        )

        sum_error = float(
            np.sum(
                self.sum_error,
                dtype=np.float64,
            )
        )

        sum_error_squared = float(
            np.sum(
                self.sum_error_squared,
                dtype=np.float64,
            )
        )

        # means
        mean_reference = (
                sum_reference / n_obs_float
        )

        mean_dataset = (
                sum_dataset / n_obs_float
        )

        # bias and mean squared error
        bias = (
                sum_error / n_obs_float
        )

        mean_squared_error = (
                sum_error_squared / n_obs_float
        )

        rmse = float(
            np.sqrt(
                max(
                    mean_squared_error,
                    0.0,
                )
            )
        )

        # variances
        variance_reference = (
                sum_reference_squared / n_obs_float
                - mean_reference ** 2
        )

        variance_dataset = (
                sum_dataset_squared / n_obs_float
                - mean_dataset ** 2
        )

        variance_reference = max(
            variance_reference,
            0.0,
        )

        variance_dataset = max(
            variance_dataset,
            0.0,
        )

        std_reference = float(
            np.sqrt(
                variance_reference
            )
        )

        std_dataset = float(
            np.sqrt(
                variance_dataset
            )
        )

        # covariance and correlation
        covariance = (
                sum_cross / n_obs_float
                - mean_reference * mean_dataset
        )

        correlation_denominator = (
                std_reference * std_dataset
        )

        if correlation_denominator > 0.0:
            correlation = float(
                np.clip(
                    covariance / correlation_denominator,
                    -1.0,
                    1.0,
                )
            )
        else:
            correlation = np.nan

        # unbiased RMSD
        error_variance = (
                mean_squared_error
                - bias ** 2
        )

        ubrmsd = float(
            np.sqrt(
                max(
                    error_variance,
                    0.0,
                )
            )
        )

        return {
            "n_obs": n_obs,
            "mean_reference": mean_reference,
            "mean_dataset": mean_dataset,
            "bias": bias,
            "rmse": rmse,
            "ubrmsd": ubrmsd,
            "correlation": correlation,
            "std_reference": std_reference,
            "std_dataset": std_dataset,
        }

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # public method to return accumulator information
    def info(self) -> Dict[str, Any]:

        if self.initialized:
            cells_with_observations = int(
                np.count_nonzero(self.count)
            )

            maximum_observations = int(
                np.max(self.count)
            )
        else:
            cells_with_observations = 0
            maximum_observations = 0

        return {
            "initialized": self.initialized,
            "grid_shape": self.grid_shape,
            "update_count": self.update_count,
            "min_observations": self.min_observations,
            "cells_with_observations": cells_with_observations,
            "maximum_observations": maximum_observations,
        }

    # ------------------------------------------------------------------------------------------------------------------
    # public method to clear cumulative arrays
    def clear(self) -> None:

        self.count = None

        self.sum_reference = None
        self.sum_dataset = None

        self.sum_reference_squared = None
        self.sum_dataset_squared = None

        self.sum_cross = None

        self.sum_error = None
        self.sum_error_squared = None

        self.initialized = False
        self.update_count = 0
        self.grid_shape = None
# ----------------------------------------------------------------------------------------------------------------------
