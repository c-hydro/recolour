"""
Library Features:

Name:          lib_weights
Author(s):     Fabio Delogu
Date:          '20260717'
Version:       '1.0.0'

Purpose:
    Compute spatial reliability weights from historical validation metrics
    and apply weighted nudging assimilation.
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging

from typing import Any, Dict, Optional

import numpy as np

from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class to compute and apply nudging weights
class NudgingWeights:

    def __init__(self,weights_cfg: Optional[Dict[str, Any]] = None,):

        # normalize configuration
        self.weights_cfg = (
            {}
            if weights_cfg is None
            else dict(weights_cfg)
        )

        # general configuration
        self.method = self.weights_cfg.get(
            "method",
            "correlation_error",
        )

        self.maximum_weight = float(
            self.weights_cfg.get(
                "maximum_weight",
                1.0,
            )
        )

        if not 0.0 <= self.maximum_weight <= 1.0:
            raise ValueError(
                "'maximum_weight' must be between zero and one."
            )

        # correlation configuration
        correlation_cfg = self.weights_cfg.get(
            "correlation",
            {},
        )

        self.minimum_correlation = float(
            correlation_cfg.get(
                "minimum",
                0.3,
            )
        )

        self.correlation_power = float(
            correlation_cfg.get(
                "power",
                1.0,
            )
        )

        # error configuration
        error_cfg = self.weights_cfg.get(
            "error",
            {},
        )

        self.error_power = float(
            error_cfg.get(
                "power",
                2.0,
            )
        )

        self.minimum_std = float(
            error_cfg.get(
                "minimum_std",
                0.01,
            )
        )

        # observation configuration
        observations_cfg = self.weights_cfg.get(
            "observations",
            {},
        )

        self.minimum_observations = int(
            observations_cfg.get(
                "minimum",
                100,
            )
        )

        self.target_observations = int(
            observations_cfg.get(
                "target",
                self.minimum_observations,
            )
        )

        if self.minimum_observations < 1:
            raise ValueError(
                "'observations.minimum' must be greater than zero."
            )

        if self.target_observations < self.minimum_observations:
            raise ValueError(
                "'observations.target' must be greater than or equal "
                "to 'observations.minimum'."
            )

        # scaling configuration
        scaling_cfg = self.weights_cfg.get(
            "scaling",
            {},
        )

        self.scaling_enabled = bool(
            scaling_cfg.get(
                "enabled",
                True,
            )
        )

        self.scaling_method = scaling_cfg.get(
            "method",
            "mean_std",
        )

        self.scaling_minimum_std = float(
            scaling_cfg.get(
                "minimum_std",
                self.minimum_std,
            )
        )

        # nudging configuration
        nudging_cfg = self.weights_cfg.get(
            "nudging",
            {},
        )

        self.default_alpha = float(
            nudging_cfg.get(
                "alpha",
                1.0,
            )
        )

        self.clip_min = nudging_cfg.get(
            "clip_min",
            None,
        )

        self.clip_max = nudging_cfg.get(
            "clip_max",
            None,
        )

    # ------------------------------------------------------------------------------------------------------------------
    # public method to compute spatial nudging weights
    def compute(
            self,
            metrics_data: Dict[str, np.ndarray],
    ) -> Dict[str, np.ndarray]:
        """
        Compute a spatial reliability weight from historical metric maps.

        Weight definition:

            final_weight =
                maximum_weight
                * correlation_weight
                * error_weight
                * observation_weight

        Returns
        -------
        dict
            Final weight and its individual components.
        """

        # get required metrics
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
        n_obs = np.asarray(
            metrics_data["n_obs"],
            dtype=np.float64,
        )

        correlation = np.asarray(
            metrics_data["correlation"],
            dtype=np.float64,
        )

        ubrmsd = np.asarray(
            metrics_data["ubrmsd"],
            dtype=np.float64,
        )

        std_reference = np.asarray(
            metrics_data["std_reference"],
            dtype=np.float64,
        )

        # check shapes
        metric_shape = n_obs.shape

        for metric_name, metric_values in {
            "correlation": correlation,
            "ubrmsd": ubrmsd,
            "std_reference": std_reference,
        }.items():

            if metric_values.shape != metric_shape:
                raise ValueError(
                    f"Metric '{metric_name}' shape {metric_values.shape} "
                    f"differs from observation shape {metric_shape}."
                )

        # initialize weight components
        correlation_weight = np.zeros(
            metric_shape,
            dtype=np.float64,
        )

        error_weight = np.zeros(
            metric_shape,
            dtype=np.float64,
        )

        observation_weight = np.zeros(
            metric_shape,
            dtype=np.float64,
        )

        # define basic valid mask
        valid_mask = (
            np.isfinite(correlation)
            & np.isfinite(ubrmsd)
            & np.isfinite(std_reference)
            & (n_obs >= self.minimum_observations)
            & (std_reference >= self.minimum_std)
        )

        # ----------------------------------------------------------------------
        # correlation component
        #
        # correlation <= minimum_correlation -> 0
        # correlation == 1                   -> 1

        correlation_range = (
            1.0 - self.minimum_correlation
        )

        if correlation_range <= 0.0:
            raise ValueError(
                "'correlation.minimum' must be lower than one."
            )

        correlation_score = (
            correlation - self.minimum_correlation
        ) / correlation_range

        correlation_score = np.clip(
            correlation_score,
            0.0,
            1.0,
        )

        correlation_weight[valid_mask] = (
            correlation_score[valid_mask]
            ** self.correlation_power
        )

        # ----------------------------------------------------------------------
        # error component
        #
        # Normalize ubRMSD using temporal variability of the reference:
        #
        # normalized_error = ubRMSD / std_reference
        #
        # weight:
        #     1 / (1 + normalized_error ** power)

        normalized_error = np.full(
            metric_shape,
            np.nan,
            dtype=np.float64,
        )

        normalized_error[valid_mask] = (
            ubrmsd[valid_mask]
            / std_reference[valid_mask]
        )

        error_weight[valid_mask] = (
            1.0
            / (
                1.0
                + normalized_error[valid_mask]
                ** self.error_power
            )
        )

        # ----------------------------------------------------------------------
        # observation component
        #
        # Below minimum -> 0
        # At target or above -> 1

        observation_denominator = max(
            self.target_observations
            - self.minimum_observations,
            1,
        )

        observation_score = (
            n_obs - self.minimum_observations
        ) / observation_denominator

        observation_score = np.clip(
            observation_score,
            0.0,
            1.0,
        )

        # when minimum and target are equal, valid cells receive full confidence
        if self.target_observations == self.minimum_observations:
            observation_score[
                n_obs >= self.minimum_observations
            ] = 1.0

        observation_weight[valid_mask] = (
            observation_score[valid_mask]
        )

        # ----------------------------------------------------------------------
        # combine components

        final_weight = (
            self.maximum_weight
            * correlation_weight
            * error_weight
            * observation_weight
        )

        # protect output
        final_weight[~valid_mask] = 0.0

        final_weight = np.clip(
            final_weight,
            0.0,
            self.maximum_weight,
        )

        logger.info(
            " ----> Nudging weights computed: "
            "valid=%d/%d, min=%.4f, mean=%.4f, max=%.4f",
            np.count_nonzero(valid_mask),
            final_weight.size,
            np.min(final_weight),
            np.mean(final_weight),
            np.max(final_weight),
        )

        return {
            "weight": final_weight,
            "weight_correlation": correlation_weight,
            "weight_error": error_weight,
            "weight_observations": observation_weight,
            "normalized_error": normalized_error,
            "valid_mask": valid_mask,
        }

    # ------------------------------------------------------------------------------------------------------------------
    # public method to scale other dataset to reference statistics
    def scale_other(
            self,
            values_other: np.ndarray,
            metrics_data: Dict[str, np.ndarray],
    ) -> np.ndarray:
        """
        Scale the current other-dataset map to historical reference statistics.

        mean/std scaling:

            other_scaled =
                mean_reference
                + (other - mean_dataset)
                * std_reference / std_dataset
        """

        values_other = np.asarray(
            values_other,
            dtype=np.float64,
        )

        if not self.scaling_enabled:
            return values_other.copy()

        if self.scaling_method != "mean_std":
            raise NotImplementedError(
                f"Scaling method '{self.scaling_method}' is not supported."
            )

        required_metrics = [
            "mean_reference",
            "mean_dataset",
            "std_reference",
            "std_dataset",
        ]

        for metric_name in required_metrics:

            if metric_name not in metrics_data:
                raise KeyError(
                    f"Metric '{metric_name}' is required for scaling."
                )

        mean_reference = np.asarray(
            metrics_data["mean_reference"],
            dtype=np.float64,
        )

        mean_dataset = np.asarray(
            metrics_data["mean_dataset"],
            dtype=np.float64,
        )

        std_reference = np.asarray(
            metrics_data["std_reference"],
            dtype=np.float64,
        )

        std_dataset = np.asarray(
            metrics_data["std_dataset"],
            dtype=np.float64,
        )

        if values_other.shape != mean_reference.shape:
            raise ValueError(
                f"Other-data shape {values_other.shape} differs from "
                f"metrics shape {mean_reference.shape}."
            )

        values_scaled = np.full(
            values_other.shape,
            np.nan,
            dtype=np.float64,
        )

        valid_mask = (
            np.isfinite(values_other)
            & np.isfinite(mean_reference)
            & np.isfinite(mean_dataset)
            & np.isfinite(std_reference)
            & np.isfinite(std_dataset)
            & (std_dataset >= self.scaling_minimum_std)
        )

        values_scaled[valid_mask] = (
            mean_reference[valid_mask]
            + (
                values_other[valid_mask]
                - mean_dataset[valid_mask]
            )
            * (
                std_reference[valid_mask]
                / std_dataset[valid_mask]
            )
        )

        return values_scaled

    # ------------------------------------------------------------------------------------------------------------------
    # public method to apply weighted nudging assimilation
    def assimilate(
            self,
            values_reference: np.ndarray,
            values_other: np.ndarray,
            metrics_data: Dict[str, np.ndarray],
            weights_data: Optional[Dict[str, np.ndarray]] = None,
            alpha: Optional[float] = None,
    ) -> Dict[str, np.ndarray]:
        """
        Apply spatially weighted nudging:

            analysis =
                reference
                + alpha * weight * (other_scaled - reference)

        The effective assimilation gain is:

            gain = alpha * weight
        """

        values_reference = np.asarray(
            values_reference,
            dtype=np.float64,
        )

        values_other = np.asarray(
            values_other,
            dtype=np.float64,
        )

        if values_reference.shape != values_other.shape:
            raise ValueError(
                f"Reference shape {values_reference.shape} differs from "
                f"other-data shape {values_other.shape}."
            )

        # compute or retrieve weights
        if weights_data is None:
            weights_data = self.compute(
                metrics_data=metrics_data,
            )

        spatial_weight = np.asarray(
            weights_data["weight"],
            dtype=np.float64,
        )

        if spatial_weight.shape != values_reference.shape:
            raise ValueError(
                f"Weight shape {spatial_weight.shape} differs from "
                f"data shape {values_reference.shape}."
            )

        # define nudging strength
        if alpha is None:
            alpha = self.default_alpha

        alpha = float(alpha)

        if not 0.0 <= alpha <= 1.0:
            raise ValueError(
                "Nudging 'alpha' must be between zero and one."
            )

        # scale other dataset
        values_other_scaled = self.scale_other(
            values_other=values_other,
            metrics_data=metrics_data,
        )

        # initialize result with reference values
        values_analysis = values_reference.copy()

        # effective gain
        effective_gain = (
            alpha * spatial_weight
        )

        # cells where the complete assimilation can be applied
        assimilation_mask = (
            np.isfinite(values_reference)
            & np.isfinite(values_other_scaled)
            & np.isfinite(effective_gain)
            & (effective_gain > 0.0)
        )

        # innovation
        innovation = np.full(
            values_reference.shape,
            np.nan,
            dtype=np.float64,
        )

        innovation[assimilation_mask] = (
            values_other_scaled[assimilation_mask]
            - values_reference[assimilation_mask]
        )

        # weighted nudging
        values_analysis[assimilation_mask] = (
            values_reference[assimilation_mask]
            + effective_gain[assimilation_mask]
            * innovation[assimilation_mask]
        )

        # optional physical clipping
        if self.clip_min is not None:
            values_analysis = np.maximum(
                values_analysis,
                float(self.clip_min),
            )

        if self.clip_max is not None:
            values_analysis = np.minimum(
                values_analysis,
                float(self.clip_max),
            )

        logger.info(
            " ----> Nudging assimilation completed: "
            "assimilated_pixels=%d/%d, alpha=%.4f",
            np.count_nonzero(assimilation_mask),
            values_analysis.size,
            alpha,
        )

        return {
            "analysis": values_analysis,
            "reference": values_reference,
            "other": values_other,
            "other_scaled": values_other_scaled,
            "innovation": innovation,
            "weight": spatial_weight,
            "gain": effective_gain,
            "assimilation_mask": assimilation_mask,
        }
# ----------------------------------------------------------------------------------------------------------------------
