
"""
Library Features:

Name:          lib_results
Author(s):     Fabio Delogu
Date:          '20260723'
Version:       '1.1.0'

Purpose:
    Plot nudging-weight maps in PNG format and save them as GeoTIFF
    or ESRI ASCII Grid files.

    Destination folders and filenames support these tags:

        {name}
        {season}
        {time:%Y}
        {time:%m}
        {time:%d}
        {time:%H}
        {time:%M}
        {time:%Y%m%d_%H%M}
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os

from typing import Any, Dict, List, Mapping, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
            time_reference: Optional[Any] = None,
            time_start: Optional[Any] = None, time_end: Optional[Any] = None,
            img_cfg: Optional[Dict[str, Any]] = None,
            results_cfg: Optional[Dict[str, Any]] = None,
    ):

        # set time tag
        self.time_tag = str(time_tag)

        # set reference time used for destination folders and filenames
        self.time_reference = self._parse_time(time_value=time_reference, time_name="time_reference",)
        if time_start is not None:
            self.time_start = self._parse_time(time_value=time_start, time_name="time_reference", )
        else:
            self.time_start = self.time_reference
        if time_end is not None:
            self.time_end = self._parse_time(time_value=time_end, time_name="time_reference", )
        else:
            self.time_end = self.time_reference

        # normalize configurations
        self.img_cfg = ({} if img_cfg is None else dict(img_cfg))
        self.results_cfg = ({} if results_cfg is None else dict(results_cfg))

        # initialize image configuration
        self._configure_images()

        # initialize result configuration
        self._configure_results()

    # ------------------------------------------------------------------------------------------------------------------
    # method to parse a time value
    @staticmethod
    def _parse_time(
            time_value: Any,
            time_name: str,
    ) -> pd.Timestamp:

        if time_value is None:
            raise ValueError(
                f"'{time_name}' is not defined."
            )

        try:
            time_obj = pd.Timestamp(
                time_value
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Unable to parse '{time_name}': '{time_value}'."
            ) from exc

        if pd.isna(time_obj):
            raise ValueError(
                f"'{time_name}' is NaT."
            )

        return time_obj

    # ------------------------------------------------------------------------------------------------------------------
    # method to resolve destination templates
    def _format_destination(
            self,
            template: str,
            field_name: str,
            layer_name: Optional[str] = None,
    ) -> str:

        if not isinstance(template, str):
            raise TypeError(
                f"Destination template '{field_name}' must be a string."
            )

        format_values = {
            "time": self.time_reference.to_pydatetime(),
            "time_reference": self.time_reference.to_pydatetime(),
            "time_start": self.time_start.to_pydatetime(),
            "time_end": self.time_end.to_pydatetime(),
            "season": self.time_tag,
            "time_tag": self.time_tag,
            "name": layer_name if layer_name is not None else "",
        }

        try:
            destination = template.format(**format_values)
        except (KeyError, ValueError, IndexError) as exc:
            raise ValueError(
                f"Unable to resolve destination template "
                f"'{field_name}': '{template}' using "
                f"time='{self.time_reference}', "
                f"season='{self.time_tag}' and "
                f"name='{layer_name}'."
            ) from exc

        return destination

    # ------------------------------------------------------------------------------------------------------------------
    # method to configure image options
    def _configure_images(self) -> None:

        # general options
        self.img_enabled = bool(self.img_cfg.get("enabled",True,))

        img_folder_template = str(self.img_cfg.get("folder","./output/images",))

        self.img_folder = self._format_destination(template=img_folder_template,field_name="img.folder",)
        self.img_filename = str(self.img_cfg.get("filename","weights_{season}_{name}_{time:%Y%m%d_%H%M}.png",))
        self.img_dpi = int(self.img_cfg.get("dpi",150,))
        self.img_figsize = self.img_cfg.get("figsize",[10, 8],)

        if (not isinstance(self.img_figsize, (list, tuple)) or len(self.img_figsize) != 2):
            raise ValueError("Image option 'figsize' must contain two values.")

        self.img_figsize = (float(self.img_figsize[0]), float(self.img_figsize[1]),)
        self.img_colormap = str(self.img_cfg.get("colormap","viridis",) )
        self.img_show_colorbar = bool(self.img_cfg.get("show_colorbar",True,))
        self.img_show_axis = bool(self.img_cfg.get("show_axis", True,))
        self.img_grid = bool(self.img_cfg.get("grid", False,))
        self.img_transpose = bool(self.img_cfg.get("transpose", False,))
        self.img_flip_vertical = bool(self.img_cfg.get("flip_vertical",False,))
        self.img_vmin = self.img_cfg.get("vmin",None,)
        self.img_vmax = self.img_cfg.get("vmax",None,)
        self.img_title = str(self.img_cfg.get("title", "Nudging weight: {name}",))
        self.img_colorbar_label = str(self.img_cfg.get("colorbar_label","Weight",))
        self.img_layers = self.img_cfg.get("layers",
            ["weight", "weight_correlation", "weight_error", "weight_observations", "normalized_error",],
        )

        if not isinstance(self.img_layers, list):
            raise TypeError("Image option 'layers' must be a list.")

        self.img_ranges = self.img_cfg.get("ranges",{},)

        if not isinstance(self.img_ranges, Mapping):
            raise TypeError("Image option 'ranges' must be a dictionary.")

    # ------------------------------------------------------------------------------------------------------------------
    # method to configure output result options
    def _configure_results(self) -> None:

        self.results_enabled = bool(
            self.results_cfg.get(
                "enabled",
                True,
            )
        )

        results_folder_template = str(
            self.results_cfg.get(
                "folder",
                "./output/results",
            )
        )

        self.results_folder = self._format_destination(
            template=results_folder_template,
            field_name="results.folder",
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
                "weights_{season}_{name}_{time:%Y%m%d_%H%M}.tif",
            )
        )

        # Compression is optional.
        # Recommended values: null, "none", "lzw", "deflate".
        compression_raw = tiff_cfg.get(
            "compression",
            "lzw",
        )

        if compression_raw is None:
            self.tiff_compression = None
        else:
            compression_value = str(
                compression_raw
            ).strip().lower()

            if compression_value in {
                "",
                "none",
                "null",
                "false",
            }:
                self.tiff_compression = None

            elif compression_value in {
                "lzw",
                "deflate",
            }:
                self.tiff_compression = compression_value

            else:
                raise ValueError(
                    "Unsupported TIFF compression "
                    f"'{compression_raw}'. "
                    "Supported values are: null, none, lzw, deflate."
                )

        # Tiling is disabled by default because these maps are relatively small.
        self.tiff_tiled = bool(
            tiff_cfg.get(
                "tiled",
                False,
            )
        )

        self.tiff_block_size = int(
            tiff_cfg.get(
                "block_size",
                256,
            )
        )

        if self.tiff_block_size <= 0:
            raise ValueError(
                "TIFF option 'block_size' must be greater than zero."
            )

        if self.tiff_block_size % 16 != 0:
            raise ValueError(
                "TIFF option 'block_size' must be a multiple of 16."
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
                "weights_{season}_{name}_{time:%Y%m%d_%H%M}.asc",
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
    def _validate_analysis(analysis_summary: Dict[str, Any],) -> None:

        if not isinstance(analysis_summary, dict):
            raise TypeError("'analysis_summary' must be a dictionary.")
        if "weights" not in analysis_summary:
            raise KeyError("Analysis summary does not contain the 'weights' section.")
        if analysis_summary["weights"] is None:
            raise ValueError("Analysis summary contains no weight data.")
        if "grid" not in analysis_summary:
            raise KeyError("Analysis summary does not contain the 'grid' section.")

    # ------------------------------------------------------------------------------------------------------------------
    # method to get selected weight layers
    @staticmethod
    def _get_layers(weights_data: Dict[str, Any], requested_layers: List[str],) -> Dict[str, np.ndarray]:

        layers_data = {}

        for layer_name in requested_layers:

            if layer_name not in weights_data:
                logger.warning(
                    " -----> Weight layer '%s' is not available and "
                    "will be skipped.",
                    layer_name,
                )
                continue

            layer_values = np.asarray(weights_data[layer_name])
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
    def _prepare_values(self, values: np.ndarray,) -> np.ndarray:

        values = np.asarray(values, dtype=np.float64,).copy()

        if self.img_transpose: values = values.T
        if self.img_flip_vertical: values = np.flipud(values)

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

        longitude = np.asarray(longitude, dtype=np.float64,)
        latitude = np.asarray(latitude, dtype=np.float64,)
        finite_lon = longitude[np.isfinite(longitude)]
        finite_lat = latitude[np.isfinite(latitude)]

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

        layer_range = self.img_ranges.get(layer_name,{},)

        if not isinstance(layer_range, Mapping):
            layer_range = {}

        vmin = layer_range.get("vmin",self.img_vmin,)
        vmax = layer_range.get("vmax",self.img_vmax,)

        return vmin, vmax

    # ------------------------------------------------------------------------------------------------------------------
    # public method to create PNG images
    def plot(self,analysis_summary: Dict[str, Any],) -> Dict[str, str]:

        self._validate_analysis(analysis_summary=analysis_summary,)

        if not self.img_enabled:
            logger.info(" ----> Weight PNG generation is disabled.")
            return {}

        logger.info(" ----> Create weight PNG images ...")

        weights_data = analysis_summary["weights"]
        grid_data = analysis_summary["grid"]

        layers_data = self._get_layers(
            weights_data=weights_data,
            requested_layers=self.img_layers,
        )

        if not layers_data:
            raise RuntimeError("No valid weight layers are available for plotting.")

        os.makedirs(self.img_folder,exist_ok=True,)

        longitude = grid_data.get("longitude")
        latitude = grid_data.get("latitude")
        extent = self._get_extent(longitude=longitude, latitude=latitude,)

        output_files = {}
        for layer_name, layer_values in layers_data.items():

            logger.info(" -----> Plot weight layer: %s",layer_name,)

            values_plot = self._prepare_values(values=layer_values,)
            vmin, vmax = self._get_plot_range(layer_name=layer_name,)

            figure, axis = plt.subplots(figsize=self.img_figsize,)

            image = axis.imshow(
                values_plot,
                origin="upper",
                extent=extent,
                cmap=self.img_colormap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )

            title = self.img_title.format(name=layer_name,
                season=self.time_tag,time=self.time_reference.to_pydatetime(),
            )
            axis.set_title(title,)

            if self.img_show_axis:
                axis.set_xlabel("Longitude")
                axis.set_ylabel("Latitude")
            else:
                axis.set_axis_off()

            if self.img_grid and self.img_show_axis:
                axis.grid(True,alpha=0.3,)

            if self.img_show_colorbar:
                colorbar = figure.colorbar(image, ax=axis,shrink=0.85,)
                colorbar.set_label(self.img_colorbar_label)

            figure.tight_layout()

            file_name = self._format_destination(
                template=self.img_filename,
                field_name="img.filename",
                layer_name=layer_name,
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
                " -----> PNG saved: %s",
                file_path,
            )

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

        logger.info(
            " -----> Collect weight layers ..."
        )

        layers_data = self._get_layers(
            weights_data=weights_data,
            requested_layers=self.results_layers,
        )

        logger.info(
            " -----> Collect weight layers ... DONE: %s",
            list(layers_data.keys()),
        )

        if not layers_data:
            raise RuntimeError(
                "No valid weight layers are available for saving."
            )

        transform = grid_data.get("transform")
        crs = grid_data.get("crs")

        if transform is None:
            raise ValueError("Grid transform is required to save raster outputs.")

        if not isinstance(transform, Affine):
            transform = Affine(*transform[:6])

        if crs is not None and not isinstance(crs, CRS):
            crs = CRS.from_user_input(crs)

        logger.info(
            " -----> Grid transform: %s",
            transform,
        )

        logger.info(
            " -----> Grid CRS: %s",
            crs,
        )

        os.makedirs(
            self.results_folder,
            exist_ok=True,
        )

        output_files = {
            "tiff": {},
            "ascii": {},
        }

        for layer_name, layer_values in layers_data.items():

            logger.info(
                " -----> Prepare weight layer: %s",
                layer_name,
            )

            # Create an independent, native-endian, C-contiguous array.
            values_output = np.array(
                layer_values,
                dtype=np.dtype(self.results_dtype).newbyteorder("="),
                copy=True,
                order="C",
            )

            if values_output.ndim != 2:
                raise ValueError(
                    f"Weight layer '{layer_name}' must be two-dimensional. "
                    f"Received shape: {values_output.shape}"
                )

            invalid_mask = ~np.isfinite(values_output)
            invalid_count = int(np.count_nonzero(invalid_mask))
            values_output[invalid_mask] = self.results_nodata

            # Recreate the array after replacing invalid values.
            values_output = np.ascontiguousarray(values_output)

            logger.info(
                " -----> Layer '%s': shape=%s, dtype=%s, "
                "C-contiguous=%s, invalid=%s, min=%s, max=%s",
                layer_name,
                values_output.shape,
                values_output.dtype,
                values_output.flags["C_CONTIGUOUS"],
                invalid_count,
                float(np.min(values_output)),
                float(np.max(values_output)),
            )

            if self.tiff_enabled:

                file_name_tiff = self._format_destination(
                    template=self.tiff_filename,
                    field_name="results.tiff.filename",
                    layer_name=layer_name,
                )

                file_path_tiff = os.path.join(
                    self.results_folder,
                    file_name_tiff,
                )

                logger.info(
                    " -----> Write TIFF layer '%s': %s",
                    layer_name,
                    file_path_tiff,
                )

                self._write_tiff(
                    file_path=file_path_tiff,
                    values=values_output,
                    transform=transform,
                    crs=crs,
                )

                logger.info(
                    " -----> TIFF saved: %s",
                    file_path_tiff,
                )

                output_files["tiff"][
                    layer_name
                ] = file_path_tiff

            if self.ascii_enabled:

                file_name_ascii = self._format_destination(
                    template=self.ascii_filename,
                    field_name="results.ascii.filename",
                    layer_name=layer_name,
                )

                file_path_ascii = os.path.join(
                    self.results_folder,
                    file_name_ascii,
                )

                logger.info(
                    " -----> Write ASCII layer '%s': %s",
                    layer_name,
                    file_path_ascii,
                )

                self._write_ascii(
                    file_path=file_path_ascii,
                    values=values_output,
                    transform=transform,
                )

                logger.info(
                    " -----> ASCII saved: %s",
                    file_path_ascii,
                )

                output_files["ascii"][
                    layer_name
                ] = file_path_ascii

            del invalid_mask
            del values_output

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

        # Create an independent, native-endian and C-contiguous array.
        values_write = np.array(
            values,
            dtype=np.dtype(self.results_dtype).newbyteorder("="),
            copy=True,
            order="C",
        )

        if values_write.ndim != 2:
            raise ValueError(
                f"GeoTIFF values must be 2D. "
                f"Received shape: {values_write.shape}"
            )

        height, width = values_write.shape

        profile = {
            "driver": "GTiff",
            "height": height,
            "width": width,
            "count": 1,
            "dtype": values_write.dtype.name,
            "crs": crs,
            "transform": transform,
            "nodata": float(self.results_nodata),
        }

        # Add compression only when configured.
        if self.tiff_compression is not None:

            profile["compress"] = self.tiff_compression

            # Floating-point predictor.
            if self.tiff_compression in {
                "lzw",
                "deflate",
            }:
                profile["predictor"] = 3

        # Add tiling only when explicitly enabled.
        if self.tiff_tiled:
            block_size_x = min(
                self.tiff_block_size,
                width,
            )

            block_size_y = min(
                self.tiff_block_size,
                height,
            )

            # Tile dimensions must be multiples of 16.
            block_size_x = max(
                16,
                block_size_x - block_size_x % 16,
            )

            block_size_y = max(
                16,
                block_size_y - block_size_y % 16,
            )

            profile.update({
                "tiled": True,
                "blockxsize": block_size_x,
                "blockysize": block_size_y,
            })

        logger.info(
            " ------> Open GeoTIFF: %s",
            file_path,
        )

        logger.info(
            " ------> GeoTIFF compression: %s",
            self.tiff_compression
            if self.tiff_compression is not None
            else "disabled",
        )

        logger.info(
            " ------> GeoTIFF tiled: %s",
            self.tiff_tiled,
        )

        logger.debug(
            " ------> GeoTIFF profile: %s",
            profile,
        )

        with rasterio.open(
                file_path,
                mode="w",
                **profile,
        ) as file_handle:

            logger.info(
                " ------> GeoTIFF opened: %s",
                file_path,
            )

            file_handle.write(
                values_write,
                indexes=1,
            )

            logger.info(
                " ------> GeoTIFF values written: %s",
                file_path,
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
    def organize(self,analysis_summary: Dict[str, Any],) -> Dict[str, Any]:

        # message organize start
        logger.info(" ----> Organize weight outputs ...")

        # info
        logger.info(" -----> Reference time: %s",self.time_reference,)
        logger.info(" -----> Start time: %s", self.time_start, )
        logger.info(" -----> End time: %s", self.time_end, )
        logger.info(" -----> Image folder: %s",self.img_folder,)
        logger.info(" -----> Results folder: %s",self.results_folder,)

        # organize images
        image_files = self.plot(analysis_summary=analysis_summary,)
        # organize files
        result_files = self.save(analysis_summary=analysis_summary,)

        # resume output info
        output_data = {
            "time_reference": self.time_reference, "time_tag": self.time_tag,
            "time_start": self.time_start, "time_end": self.time_end,
            "images": image_files, "results": result_files,
        }

        # message organize end
        logger.info(" ----> Organize weight outputs ... DONE")

        return output_data
# ----------------------------------------------------------------------------------------------------------------------

