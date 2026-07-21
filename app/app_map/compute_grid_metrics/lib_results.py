"""
Library Features:

Name:          lib_results
Author(s):     Fabio Delogu
Date:          '20260717'
Version:       '1.0.0'

Purpose:
    Plot nudging-weight maps in PNG format and save them as GeoTIFF
    or ESRI ASCII Grid files.
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from typing import Any, Dict, List, Mapping, Optional

import matplotlib.pyplot as plt
import numpy as np
import rasterio

from rasterio.crs import CRS
from rasterio.transform import Affine

from config_info import LOGGER_NAME

# logging
logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# class to view and save nudging weight maps
class Results:

    def __init__(
            self,
            time_tag: str = "ALL",
            img_cfg: Optional[Dict[str, Any]] = None,
            results_cfg: Optional[Dict[str, Any]] = None,
    ):

        # set time tag
        self.time_tag = time_tag

        # normalize configurations
        self.img_cfg = (
            {}
            if img_cfg is None
            else dict(img_cfg)
        )

        self.results_cfg = (
            {}
            if results_cfg is None
            else dict(results_cfg)
        )

        # initialize image configuration
        self._configure_images()

        # initialize result configuration
        self._configure_results()

    # ------------------------------------------------------------------------------------------------------------------
    # method to configure image options
    def _configure_images(self) -> None:

        # general options
        self.img_enabled = bool(
            self.img_cfg.get(
                "enabled",
                True,
            )
        )

        self.img_folder = str(
            self.img_cfg.get(
                "folder",
                "./output/images",
            )
        )

        self.img_filename = str(
            self.img_cfg.get(
                "filename",
                "weights_{name}.png",
            )
        )

        self.img_dpi = int(
            self.img_cfg.get(
                "dpi",
                150,
            )
        )

        self.img_figsize = self.img_cfg.get(
            "figsize",
            [10, 8],
        )

        if (
                not isinstance(self.img_figsize, (list, tuple))
                or len(self.img_figsize) != 2
        ):
            raise ValueError(
                "Image option 'figsize' must contain two values."
            )

        self.img_figsize = (
            float(self.img_figsize[0]),
            float(self.img_figsize[1]),
        )

        self.img_colormap = str(
            self.img_cfg.get(
                "colormap",
                "viridis",
            )
        )

        self.img_show_colorbar = bool(
            self.img_cfg.get(
                "show_colorbar",
                True,
            )
        )

        self.img_show_axis = bool(
            self.img_cfg.get(
                "show_axis",
                True,
            )
        )

        self.img_grid = bool(
            self.img_cfg.get(
                "grid",
                False,
            )
        )

        self.img_transpose = bool(
            self.img_cfg.get(
                "transpose",
                False,
            )
        )

        self.img_flip_vertical = bool(
            self.img_cfg.get(
                "flip_vertical",
                False,
            )
        )

        self.img_vmin = self.img_cfg.get(
            "vmin",
            None,
        )

        self.img_vmax = self.img_cfg.get(
            "vmax",
            None,
        )

        self.img_title = str(
            self.img_cfg.get(
                "title",
                "Nudging weight: {name}",
            )
        )

        self.img_colorbar_label = str(
            self.img_cfg.get(
                "colorbar_label",
                "Weight",
            )
        )

        self.img_layers = self.img_cfg.get(
            "layers",
            [
                "weight",
                "weight_correlation",
                "weight_error",
                "weight_observations",
                "normalized_error",
            ],
        )

        if not isinstance(self.img_layers, list):
            raise TypeError(
                "Image option 'layers' must be a list."
            )

        # layer-specific plot ranges
        self.img_ranges = self.img_cfg.get(
            "ranges",
            {},
        )

        if not isinstance(self.img_ranges, Mapping):
            raise TypeError(
                "Image option 'ranges' must be a dictionary."
            )

    # ------------------------------------------------------------------------------------------------------------------
    # method to configure output result options
    def _configure_results(self) -> None:

        self.results_enabled = bool(
            self.results_cfg.get(
                "enabled",
                True,
            )
        )

        self.results_folder = str(
            self.results_cfg.get(
                "folder",
                "./output/results",
            )
        )

        self.results_layers = self.results_cfg.get(
            "layers",
            [
                "weight",
                "weight_correlation",
                "weight_error",
                "weight_observations",
                "normalized_error",
            ],
        )

        if not isinstance(self.results_layers, list):
            raise TypeError(
                "Results option 'layers' must be a list."
            )

        self.results_nodata = float(
            self.results_cfg.get(
                "nodata",
                -9999.0,
            )
        )

        self.results_dtype = np.dtype(
            self.results_cfg.get(
                "dtype",
                "float32",
            )
        )

        # GeoTIFF options
        tiff_cfg = self.results_cfg.get(
            "tiff",
            {},
        )

        if tiff_cfg is None:
            tiff_cfg = {}

        if not isinstance(tiff_cfg, Mapping):
            raise TypeError(
                "Results option 'tiff' must be a dictionary."
            )

        self.tiff_enabled = bool(
            tiff_cfg.get(
                "enabled",
                True,
            )
        )

        self.tiff_filename = str(
            tiff_cfg.get(
                "filename",
                "weights_{name}.tif",
            )
        )

        self.tiff_compression = tiff_cfg.get(
            "compression",
            "deflate",
        )

        self.tiff_tiled = bool(
            tiff_cfg.get(
                "tiled",
                True,
            )
        )

        # ASCII options
        ascii_cfg = self.results_cfg.get(
            "ascii",
            {},
        )

        if ascii_cfg is None:
            ascii_cfg = {}

        if not isinstance(ascii_cfg, Mapping):
            raise TypeError(
                "Results option 'ascii' must be a dictionary."
            )

        self.ascii_enabled = bool(
            ascii_cfg.get(
                "enabled",
                False,
            )
        )

        self.ascii_filename = str(
            ascii_cfg.get(
                "filename",
                "weights_{name}.asc",
            )
        )

        self.ascii_decimals = int(
            ascii_cfg.get(
                "decimals",
                6,
            )
        )

    # ------------------------------------------------------------------------------------------------------------------
    # method to validate analysis summary
    @staticmethod
    def _validate_analysis(
            analysis_summary: Dict[str, Any],
    ) -> None:

        if not isinstance(analysis_summary, dict):
            raise TypeError(
                "'analysis_summary' must be a dictionary."
            )

        if "weights" not in analysis_summary:
            raise KeyError(
                "Analysis summary does not contain the 'weights' section."
            )

        if analysis_summary["weights"] is None:
            raise ValueError(
                "Analysis summary contains no weight data."
            )

        if "grid" not in analysis_summary:
            raise KeyError(
                "Analysis summary does not contain the 'grid' section."
            )

    # ------------------------------------------------------------------------------------------------------------------
    # method to get selected weight layers
    @staticmethod
    def _get_layers(
            weights_data: Dict[str, Any],
            requested_layers: List[str],
    ) -> Dict[str, np.ndarray]:

        layers_data = {}

        for layer_name in requested_layers:

            if layer_name not in weights_data:

                logger.warning(
                    " -----> Weight layer '%s' is not available and "
                    "will be skipped.",
                    layer_name,
                )

                continue

            layer_values = np.asarray(
                weights_data[layer_name]
            )

            if layer_values.ndim != 2:

                logger.warning(
                    " -----> Weight layer '%s' is not 2D and will "
                    "be skipped. Shape: %s",
                    layer_name,
                    layer_values.shape,
                )

                continue

            layers_data[layer_name] = layer_values

        return layers_data

    # ------------------------------------------------------------------------------------------------------------------
    # method to prepare values for plotting or saving
    def _prepare_values(
            self,
            values: np.ndarray,
    ) -> np.ndarray:

        values = np.asarray(
            values,
            dtype=np.float64,
        ).copy()

        if self.img_transpose:
            values = values.T

        if self.img_flip_vertical:
            values = np.flipud(values)

        return values

    # ------------------------------------------------------------------------------------------------------------------
    # method to get map extent from longitude and latitude
    @staticmethod
    def _get_extent(
            longitude: Optional[np.ndarray],
            latitude: Optional[np.ndarray],
    ) -> Optional[List[float]]:

        if longitude is None or latitude is None:
            return None

        longitude = np.asarray(
            longitude,
            dtype=np.float64,
        )

        latitude = np.asarray(
            latitude,
            dtype=np.float64,
        )

        finite_lon = longitude[
            np.isfinite(longitude)
        ]

        finite_lat = latitude[
            np.isfinite(latitude)
        ]

        if finite_lon.size == 0 or finite_lat.size == 0:
            return None

        return [
            float(np.min(finite_lon)),
            float(np.max(finite_lon)),
            float(np.min(finite_lat)),
            float(np.max(finite_lat)),
        ]

    # ------------------------------------------------------------------------------------------------------------------
    # method to get plot range
    def _get_plot_range(
            self,
            layer_name: str,
    ):

        layer_range = self.img_ranges.get(
            layer_name,
            {},
        )

        if not isinstance(layer_range, Mapping):
            layer_range = {}

        vmin = layer_range.get(
            "vmin",
            self.img_vmin,
        )

        vmax = layer_range.get(
            "vmax",
            self.img_vmax,
        )

        return vmin, vmax

    # ------------------------------------------------------------------------------------------------------------------
    # public method to create PNG images
    def plot(
            self,
            analysis_summary: Dict[str, Any],
    ) -> Dict[str, str]:

        self._validate_analysis(
            analysis_summary=analysis_summary,
        )

        if not self.img_enabled:

            logger.info(
                " ----> Weight PNG generation is disabled."
            )

            return {}

        logger.info(
            " ----> Create weight PNG images ..."
        )

        weights_data = analysis_summary["weights"]
        grid_data = analysis_summary["grid"]

        layers_data = self._get_layers(
            weights_data=weights_data,
            requested_layers=self.img_layers,
        )

        if not layers_data:
            raise RuntimeError(
                "No valid weight layers are available for plotting."
            )

        os.makedirs(
            self.img_folder,
            exist_ok=True,
        )

        longitude = grid_data.get(
            "longitude"
        )

        latitude = grid_data.get(
            "latitude"
        )

        extent = self._get_extent(
            longitude=longitude,
            latitude=latitude,
        )

        season_name = self.time_tag
        output_files = {}
        for layer_name, layer_values in layers_data.items():

            logger.info(
                " -----> Plot weight layer: %s",
                layer_name,
            )

            values_plot = self._prepare_values(
                values=layer_values,
            )

            vmin, vmax = self._get_plot_range(
                layer_name=layer_name,
            )

            figure, axis = plt.subplots(
                figsize=self.img_figsize,
            )

            image = axis.imshow(
                values_plot,
                origin="upper",
                extent=extent,
                cmap=self.img_colormap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )

            title = self.img_title.format(
                name=layer_name,
            )

            axis.set_title(
                title,
            )

            if self.img_show_axis:

                axis.set_xlabel(
                    "Longitude"
                )

                axis.set_ylabel(
                    "Latitude"
                )

            else:

                axis.set_axis_off()

            if self.img_grid and self.img_show_axis:
                axis.grid(
                    True,
                    alpha=0.3,
                )

            if self.img_show_colorbar:

                colorbar = figure.colorbar(
                    image,
                    ax=axis,
                    shrink=0.85,
                )

                colorbar.set_label(
                    self.img_colorbar_label
                )

            figure.tight_layout()

            file_name = self.img_filename.format(
                name=layer_name, season=season_name
            )

            file_path = os.path.join(
                self.img_folder,
                file_name,
            )

            figure.savefig(
                file_path,
                dpi=self.img_dpi,
                bbox_inches="tight",
            )

            plt.close(
                figure
            )

            output_files[layer_name] = file_path

        logger.info(
            " ----> Create weight PNG images ... DONE"
        )

        return output_files

    # ------------------------------------------------------------------------------------------------------------------
    # public method to save weight layers
    def save(
            self,
            analysis_summary: Dict[str, Any],
    ) -> Dict[str, Dict[str, str]]:

        self._validate_analysis(
            analysis_summary=analysis_summary,
        )

        if not self.results_enabled:

            logger.info(
                " ----> Weight result export is disabled."
            )

            return {
                "tiff": {},
                "ascii": {},
            }

        logger.info(
            " ----> Save weight result maps ..."
        )

        weights_data = analysis_summary["weights"]
        grid_data = analysis_summary["grid"]

        layers_data = self._get_layers(
            weights_data=weights_data,
            requested_layers=self.results_layers,
        )

        if not layers_data:
            raise RuntimeError(
                "No valid weight layers are available for saving."
            )

        transform = grid_data.get(
            "transform"
        )

        crs = grid_data.get(
            "crs"
        )

        if transform is None:
            raise ValueError(
                "Grid transform is required to save raster outputs."
            )

        if not isinstance(transform, Affine):
            transform = Affine(
                *transform[:6]
            )

        if crs is not None and not isinstance(crs, CRS):
            crs = CRS.from_user_input(
                crs
            )

        os.makedirs(
            self.results_folder,
            exist_ok=True,
        )

        output_files = {
            "tiff": {},
            "ascii": {},
        }

        season_name = self.time_tag
        for layer_name, layer_values in layers_data.items():

            values_output = np.asarray(
                layer_values,
                dtype=self.results_dtype,
            ).copy()

            invalid_mask = ~np.isfinite(
                values_output
            )

            values_output[
                invalid_mask
            ] = self.results_nodata

            if self.tiff_enabled:

                file_name_tiff = self.tiff_filename.format(
                    name=layer_name, season=season_name
                )

                file_path_tiff = os.path.join(
                    self.results_folder,
                    file_name_tiff,
                )

                self._write_tiff(
                    file_path=file_path_tiff,
                    values=values_output,
                    transform=transform,
                    crs=crs,
                )

                output_files["tiff"][
                    layer_name
                ] = file_path_tiff

            if self.ascii_enabled:

                file_name_ascii = self.ascii_filename.format(
                    name=layer_name, season=season_name
                )

                file_path_ascii = os.path.join(
                    self.results_folder,
                    file_name_ascii,
                )

                self._write_ascii(
                    file_path=file_path_ascii,
                    values=values_output,
                    transform=transform,
                )

                output_files["ascii"][
                    layer_name
                ] = file_path_ascii

        logger.info(
            " ----> Save weight result maps ... DONE"
        )

        return output_files

    # ------------------------------------------------------------------------------------------------------------------
    # method to write GeoTIFF
    def _write_tiff(
            self,
            file_path: str,
            values: np.ndarray,
            transform: Affine,
            crs: Optional[CRS],
    ) -> None:

        height, width = values.shape

        profile = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": self.results_dtype.name,
            "crs": crs,
            "transform": transform,
            "nodata": self.results_nodata,
            "compress": self.tiff_compression,
            "tiled": self.tiff_tiled,
        }

        with rasterio.open(
                file_path,
                mode="w",
                **profile,
        ) as file_handle:

            file_handle.write(
                values,
                1,
            )

        logger.info(
            " -----> GeoTIFF saved: %s",
            file_path,
        )

    # ------------------------------------------------------------------------------------------------------------------
    # method to write ESRI ASCII Grid
    def _write_ascii(
            self,
            file_path: str,
            values: np.ndarray,
            transform: Affine,
    ) -> None:

        nrows, ncols = values.shape

        cell_size_x = float(
            transform.a
        )

        cell_size_y = abs(
            float(transform.e)
        )

        if not np.isclose(
                cell_size_x,
                cell_size_y,
                rtol=1e-6,
                atol=1e-12,
        ):
            raise ValueError(
                "ESRI ASCII output requires square pixels. "
                f"Pixel sizes are x={cell_size_x}, y={cell_size_y}."
            )

        xllcorner = float(
            transform.c
        )

        yllcorner = float(
            transform.f
            + transform.e * nrows
        )

        number_format = (
            f"%.{self.ascii_decimals}f"
        )

        with open(
                file_path,
                mode="w",
                encoding="utf-8",
        ) as file_handle:

            file_handle.write(
                f"ncols         {ncols}\n"
            )

            file_handle.write(
                f"nrows         {nrows}\n"
            )

            file_handle.write(
                f"xllcorner     {xllcorner:.12f}\n"
            )

            file_handle.write(
                f"yllcorner     {yllcorner:.12f}\n"
            )

            file_handle.write(
                f"cellsize      {cell_size_x:.12f}\n"
            )

            file_handle.write(
                f"NODATA_value  {self.results_nodata}\n"
            )

            np.savetxt(
                file_handle,
                values,
                fmt=number_format,
            )

        logger.info(
            " -----> ASCII grid saved: %s",
            file_path,
        )

    # ------------------------------------------------------------------------------------------------------------------
    # public method to create all outputs
    def organize(
            self,
            analysis_summary: Dict[str, Any],
    ) -> Dict[str, Any]:

        logger.info(
            " ----> Organize weight outputs ..."
        )

        image_files = self.plot(
            analysis_summary=analysis_summary,
        )

        result_files = self.save(
            analysis_summary=analysis_summary,
        )

        output_data = {
            "images": image_files,
            "results": result_files,
        }

        logger.info(
            " ----> Organize weight outputs ... DONE"
        )

        return output_data
# ----------------------------------------------------------------------------------------------------------------------
