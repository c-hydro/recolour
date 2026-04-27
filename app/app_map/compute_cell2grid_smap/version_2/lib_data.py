"""
Library Features:

Name:           lib_process
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260421'
Version:        '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import os
import numpy as np
import pandas as pd

from lib_utils_basic import discover_source_files, resolve_generic_tags
from lib_utils_io import load_target_grid, write_output_map
from lib_utils_report import collect_report, save_report
from lib_utils_geo import map_points_to_grid_indices
from lib_utils_time import resolve_time_tags, resolve_time_window
from lib_utils_analysis_points import collect_points_to_dataframe, deduplicate_latest_points
from lib_utils_analysis_grid import (interpolate_points_to_grid,
                                     build_mask_by_pixel_extension, build_mask_boundary, build_smooth_map,
                                     apply_mask_filter)

from config_utils import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Helper to process datasets
def process(settings, reference_time):
    """
    Process source point datasets into gridded output maps.

    Workflow
    --------
    1. Read processing settings
    2. Resolve the time window
    3. Discover source files
    4. Collect and deduplicate source points
    5. Interpolate values and observation times to the target grid
    6. Optionally apply:
       - pixel-extension mask
       - boundary mask
       - smoothing
    7. Compute the time-lag map with respect to reference_time

    Returns
    -------
    soil_moisture_map_final : np.ndarray or None
        Final processed soil moisture map.
    time_lag_map : np.ndarray or None
        Time lag map in hours.
    profile : dict or None
        Output raster profile.
    stats : dict
        Summary statistics of the processing.
    time_start : datetime
        Start time of the selected processing window.
    time_end : datetime
        End time of the selected processing window.
    """

    # -------------------------------------------------------------------------------------
    # start message - script
    logger.info(" ----> Process execution ... ")

    # Read settings sections
    parameters = settings.get("parameters", {})
    source = settings.get("source", {})
    destination = settings.get("destination", {})

    variable_name = source.get("variable_name", "surface_soil_moisture")
    dry_run = bool(destination.get("dry_run", False))

    # Target grid / mask settings
    mask_file = parameters["mask_file"]
    mask_band = int(parameters.get("mask_band", 1))

    # Optional processing flags
    apply_mask_extension = int(parameters.get("apply_mask_extension", False))
    extension_pixels = int(parameters.get("extension_pixels", 0))

    apply_mask_boundary = int(parameters.get("apply_mask_boundary", False))
    boundary_pixels = int(parameters.get("boundary_pixels", 0))

    apply_smoothing = int(parameters.get("apply_smoothing", False))
    smoothing_method = parameters.get("smoothing_method", "mean")
    smoothing_sigma = int(parameters.get("smoothing_sigma", 1))

    # Interpolation radius of influence
    roi_km = float(parameters.get("roi_km", 0.0))

    # Fill value used for nodata cells
    fill_value_raw = parameters.get("fill_value", None)
    fill_value = np.nan if fill_value_raw is None else float(fill_value_raw)

    # start message - configuration
    logger.info(" ----> Configuration ... ")
    logger.info(f" -----> Variable: {variable_name}")
    logger.info(f" -----> ROI (km): {roi_km}")
    logger.info(f" -----> Fill value: {fill_value}")
    logger.info(f" -----> Mask extension: {bool(apply_mask_extension)} (pixels={extension_pixels})")
    logger.info(f" -----> Boundary mask: {bool(apply_mask_boundary)} (pixels={boundary_pixels})")
    logger.info(f" -----> Smoothing: {bool(apply_smoothing)} "
                f"(method={smoothing_method}, sigma={smoothing_sigma})")

    # Initialize processing statistics
    stats = {
        "selected_files": 0,
        "raw_points": 0,
        "deduplicated_points": 0,
        "valid_domain_pixels": 0,
        "pixel_extension_mask_pixels": 0,
        "written_files": 0,
        "errors": 0,
        "roi_km": roi_km,
        "apply_mask_extension": apply_mask_extension,
        "extension_pixels": extension_pixels,
        "apply_mask_boundary": apply_mask_boundary,
        "boundary_pixels": boundary_pixels,
        "apply_smoothing": apply_smoothing,
        "smoothing_sigma": smoothing_sigma,
        "smoothing_method": smoothing_method,
    }

    # end message - configuration
    logger.info(" ----> Configuration ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # Resolve time window for file selection
    time_start, time_end = resolve_time_window(settings, reference_time)

    # message - times
    logger.info(" ----> Times ... ")
    logger.info(f" -----> Time reference: {reference_time}")
    logger.info(f" -----> Time start: {time_start}")
    logger.info(f" -----> Time end:   {time_end}")
    logger.info(" ----> Times ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # start message - discover source file(s)
    logger.info(" ----> Discover source files ...")

    # discover source files within the requested time window
    source_files = discover_source_files(settings, time_start, time_end, reference_time)
    stats["selected_files"] = len(source_files)

    # list selected files
    logger.info(f" -----> Selected files: {len(source_files)} ... ")
    for file_path in source_files:
        logger.info(f" ------> {file_path}")
    logger.info(f" -----> Selected files: {len(source_files)} ... DONE")

    # Case 1: dry-run mode
    if dry_run:
        logger.info(" ----> Discover source files ... DONE")
        logger.info(" ----> Process execution ... DRY RUN. EXIT WITHOUT WRITING PRODUCT")
        return None, None, None, stats, time_start, time_end

    # Case 2: no files available
    if not source_files:
        logger.info(" ----> Discover source files ... NO FILES")
        logger.info(" ----> Process execution ... STOP. NOTHING TO DO")
        return None, None, None, stats, time_start, time_end

    # Case 3: invalid interpolation parameter
    if roi_km <= 0:
        raise RuntimeError("Parameter 'roi_km' must be > 0")

    # end message - discover source file(s)
    logger.info(" ----> Discover source files ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # Collect point data from source files
    logger.info(" ----> Collect points from source files ...")
    points_df = collect_points_to_dataframe(settings, source_files)
    stats["raw_points"] = len(points_df)

    # Case 4: files exist but no valid points were extracted
    if points_df.empty:
        logger.info(" ----> Collect points from source files ... NO VALID POINTS")
        logger.info(" ----> Process execution ... STOP. NOTHING TO DO")
        return None, None, None, stats, time_start, time_end

    # Deduplicate points by keeping the latest value
    logger.info(" ----> Deduplicate points ...")
    points_df = deduplicate_latest_points(points_df, value_var=variable_name)
    stats["deduplicated_points"] = len(points_df)
    logger.info(" ----> Deduplicate points ...")

    # points summary info
    logger.info(" ----> Points collection summary ... ")
    logger.info(f" -----> Raw points: {stats['raw_points']}")
    logger.info(f" -----> Deduplicated points: {stats['deduplicated_points']}")
    if stats["deduplicated_points"] < stats["raw_points"]:
        logger.info(f" -----> Removed duplicates: {stats['raw_points'] - stats['deduplicated_points']}")
    logger.info(" ----> Points collection summary ... DONE")
    # -------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - target grid
    logger.info(" ----> Loading target grid ...")

    # Load target grid and domain mask
    grid_lons, grid_lats, domain_mask, profile = load_target_grid(
        mask_file=mask_file,
        mask_band=mask_band
    )
    stats["valid_domain_pixels"] = int(domain_mask.sum())

    # grid info
    logger.info(f" -----> Grid shape: {grid_lons.shape}")
    logger.info(f" -----> Valid domain pixels: {stats['valid_domain_pixels']}")

    # message end - target grid
    logger.info(" ----> Loading target grid ... DONE")
    # -------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - interpolate points
    logger.info(" ----> Interpole points to grid ...")

    # extract point arrays
    point_lons = points_df["lon"].to_numpy(np.float64)
    point_lats = points_df["lat"].to_numpy(np.float64)
    point_values = points_df[variable_name].to_numpy(np.float64)

    point_time_seconds = (
        pd.to_datetime(points_df["time"]).astype("int64").to_numpy(dtype=np.float64) / 1.0e9
    )

    # Interpolate point values to target grid
    logger.info(f" -----> Number of input points: {len(point_lons)}")
    logger.info(f" -----> Grid size: {grid_lons.shape}")

    # interpolate point sm to target grid
    soil_moisture_map_interp = interpolate_points_to_grid(
        src_lons=point_lons,
        src_lats=point_lats,
        src_vals=point_values,
        grid_lons=grid_lons,
        grid_lats=grid_lats,
        domain_mask=domain_mask,
        roi_km=roi_km,
        fill_value=fill_value
    )

    # interpolate point times to target grid
    time_seconds_map_interp = interpolate_points_to_grid(
        src_lons=point_lons,
        src_lats=point_lats,
        src_vals=point_time_seconds,
        grid_lons=grid_lons,
        grid_lats=grid_lats,
        domain_mask=domain_mask,
        roi_km=roi_km,
        fill_value=fill_value
    )

    valid_interp = np.isfinite(soil_moisture_map_interp).sum()
    logger.info(f" -----> Interpolated valid pixels: {valid_interp}")

    # message end - interpolate points
    logger.info(" ----> Interpole points to grid ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - apply pixel-extension mask
    logger.info(" ----> Apply pixel-extension mask ...")
    # optional case: apply pixel-extension mask around source points
    if apply_mask_extension:

        # convert source points to nearest grid indices
        point_rows, point_cols = map_points_to_grid_indices(
            src_lons=point_lons,
            src_lats=point_lats,
            grid_lons=grid_lons,
            grid_lats=grid_lats,
            domain_mask=domain_mask
        )

        # build mask by expanding around each source point
        extension_mask = build_mask_by_pixel_extension(
            rows=point_rows,
            cols=point_cols,
            domain_mask=domain_mask,
            radius_pixels=extension_pixels
        )
        stats["pixel_extension_mask_pixels"] = int(extension_mask.sum())

        # apply mask filters
        soil_moisture_map_filtered = apply_mask_filter(
            data=soil_moisture_map_interp,
            filter_mask=extension_mask,
            domain_mask=domain_mask,
            fill_value=fill_value
        )
        time_seconds_map_filtered = apply_mask_filter(
            data=time_seconds_map_interp,
            filter_mask=extension_mask,
            domain_mask=domain_mask,
            fill_value=fill_value
        )

        # mask info
        logger.info(f" -----> Extension radius (pixels): {extension_pixels}")
        logger.info(f" -----> Pixel-extension mask pixels: {stats['pixel_extension_mask_pixels']}")
        if stats["valid_domain_pixels"] > 0:
            logger.info(
                f" -----> Mask coverage (%): "
                f"{100.0 * stats['pixel_extension_mask_pixels'] / stats['valid_domain_pixels']:.2f}"
            )

        # message end - apply pixel-extension mask
        logger.info(" ----> Apply pixel-extension mask ... DONE")

    else:

        # copy the interpolated maps
        soil_moisture_map_filtered = soil_moisture_map_interp.copy()
        time_seconds_map_filtered = time_seconds_map_interp.copy()

        # message end - apply pixel-extension mask
        logger.info(" ----> Apply pixel-extension mask ... NOT ACTIVE")
    # ---------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------
    # message start - apply boundary mask
    logger.info(" ----> Apply boundary mask ...")
    # Optional case: apply boundary mask
    if apply_mask_boundary:

        # compute boundary mask
        boundary_mask, boundary_map, valid_data, no_data = build_mask_boundary(
            data_map=soil_moisture_map_filtered,
            domain_mask=domain_mask,
            fill_value=fill_value,
            radius=boundary_pixels
        )
        # apply mask to sm
        soil_moisture_map_masked = soil_moisture_map_filtered.copy()
        soil_moisture_map_masked[boundary_mask] = np.float32(fill_value)
        soil_moisture_map_masked[~domain_mask] = np.float32(fill_value)
        # apply mask to times
        time_seconds_map_masked = time_seconds_map_filtered.copy()
        time_seconds_map_masked[boundary_mask] = np.float32(fill_value)
        time_seconds_map_masked[~domain_mask] = np.float32(fill_value)

        # mask info
        removed_pixels = int(np.sum(boundary_mask))
        logger.info(f" -----> Boundary radius (pixels): {boundary_pixels}")
        logger.info(f" -----> Boundary pixels removed: {removed_pixels}")

        # message end - apply boundary mask
        logger.info(" ----> Apply boundary mask ... DONE")

    else:
        # copy filtered maps
        soil_moisture_map_masked = soil_moisture_map_filtered.copy()
        time_seconds_map_masked = time_seconds_map_filtered.copy()
        # message end - apply boundary mask
        logger.info(" ----> Apply boundary mask ... NOT ACTIVE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - apply smoothing
    logger.info(" ----> Apply smoothing ...")
    # optional case: smooth the soil moisture map
    if apply_smoothing:

        # apply smoothing
        soil_moisture_map_smoothed = build_smooth_map(
            data_map=soil_moisture_map_masked,
            domain_mask=domain_mask,
            fill_value=fill_value,
            method=smoothing_method,
            sigma=smoothing_sigma
        )

        # smoothing info
        logger.info(f" -----> Method: {smoothing_method}")
        logger.info(f" -----> Sigma: {smoothing_sigma}")

        # message end - apply smoothing
        logger.info(" ----> Apply smoothing ... DONW")

    else:
        # copy sm masked map
        soil_moisture_map_smoothed = soil_moisture_map_masked.copy()
        # message end - apply smoothing
        logger.info(" ----> Apply smoothing ... NOT ACTIVE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - compute time-lag
    logger.info(" ----> Compute time-lag ... ")

    # build time lag map in hours:
    # time_lag = reference_time - interpolated observation time
    reference_time_seconds = reference_time.timestamp()

    # IMPORTANT: start from masked interpolated times, not from an empty fill-value array
    time_lag_map = time_seconds_map_masked.copy().astype(np.float32)

    if np.isnan(fill_value):
        valid_time_mask = np.isfinite(time_lag_map) & domain_mask
    else:
        valid_time_mask = (
            np.isfinite(time_lag_map) &
            domain_mask &
            (~np.isclose(time_lag_map, fill_value, atol=1e-6))
        )

    time_lag_map[valid_time_mask] = (
        (reference_time_seconds - time_lag_map[valid_time_mask]) / 3600.0
    ).astype(np.float32)

    time_lag_map[~domain_mask] = np.float32(fill_value)

    valid_time_pixels = int(np.sum(valid_time_mask))
    logger.info(f" -----> Valid time-lag pixels: {valid_time_pixels}")

    # message end - compute time-lag
    logger.info(" ----> Compute time-lag ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # end message - script
    logger.info(f" ----> Output grid shape: {soil_moisture_map_smoothed.shape}")
    logger.info(f" ----> Total selected files: {stats['selected_files']}")
    logger.info(f" ----> Total points used: {stats['deduplicated_points']}")
    logger.info(" ----> Process execution ... DONE")

    return (
        soil_moisture_map_smoothed,
        time_lag_map,
        profile,
        stats,
        time_start,
        time_end,
    )
    # ---------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Helper to save processed outputs
def save(soil_moisture_map_out, time_lag_map,
         profile,
         stats, time_start, time_end,
         time_reference, settings,
         alg_start):
    """
    Save processed output maps and processing report.

    Workflow
    --------
    1. Read destination and output settings
    2. Resolve output folder and filename using time tags
    3. Check overwrite policy
    4. Write output raster
    5. Build and save processing report

    Returns
    -------
    bool
        True if the output file was written successfully, False if skipped.
    """

    # ---------------------------------------------------------------------------------
    # start message - script
    logger.info(" ----> Dumping execution ... ")
    logger.info(f" ----> Reference time: {time_reference}")

    # read settings sections
    parameters = settings.get("parameters", {})
    destination = settings.get("destination", {})
    report = settings.get("report", {})

    # fill value used for nodata cells
    fill_value_raw = parameters.get("fill_value", None)
    fill_value = np.nan if fill_value_raw is None else float(fill_value_raw)

    # destination output settings
    variable_name = destination.get("variable_name", "surface_soil_moisture")
    dry_run = bool(destination.get("dry_run", False))
    overwrite_output = bool(destination.get("overwrite", True))

    report_enabled = report.get("enabled", False)
    report_folder_template = report.get("folder", "./report")
    report_filename_template = report.get(
        "filename", "cell2img_h122_report_%Y%m%d_%H%M.json"
    )

    # configuration info
    logger.info(" ----> Configuration ... ")
    logger.info(f" -----> Variable name: {variable_name}")
    logger.info(f" -----> Fill value: {fill_value}")
    logger.info(f" -----> Dry run: {dry_run}")
    logger.info(f" -----> Overwrite existing file: {overwrite_output}")
    logger.info(" ----> Configuration ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # Case 1: nothing to save because output maps are missing
    if soil_moisture_map_out is None or time_lag_map is None or profile is None:
        # end message - script
        logger.warning(" ----> Output data are not available. Nothing to save")
        logger.info(" ----> Dumping execution ... STOP. NO OUTPUT DATA")
        return False
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # Case 2: dry-run mode enabled
    if dry_run:
        # end message - script
        logger.info(" ----> Dumping execution ... DRY-RUN")
        return False
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - destination path
    logger.info(" ----> Define destination path ...")
    # define destination folder and filename using time tags
    destination_folder = resolve_time_tags(
        destination["folder"],
        {
            "time_step": time_reference,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": time_reference
        }
    )

    destination_filename = resolve_generic_tags(
        destination["filename"],
        {
            "time_step": time_reference,
            "time_start": time_start,
            "time_end": time_end,
            "time_run": time_reference
        }
    )

    destination_file = os.path.join(destination_folder, destination_filename)

    # output info
    logger.info(f" -----> Destination folder: {destination_folder}")
    logger.info(f" -----> Destination filename: {destination_filename}")
    logger.info(f" -----> Destination file: {destination_file}")

    # Create destination folder if needed
    if not os.path.exists(destination_folder):
        logger.info(" ----> Destination folder does not exist. Creating folder...")
        os.makedirs(destination_folder, exist_ok=True)
        logger.info(" ----> Destination folder created")

    # message end - destination path
    logger.info(" ----> Define destination path ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - destination path
    logger.info(" ----> Define report path ... ")

    report_folder = resolve_time_tags(
        report_folder_template,
        {
            "time_step": time_reference,
            "time_start": time_reference,
            "time_end": time_reference,
            "time_run": time_reference
        }
    )

    report_filename = time_reference.strftime(report_filename_template)
    report_path = os.path.join(report_folder, report_filename)

    # message end - destination path
    logger.info(" ----> Define report path ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # Case 3: output file exists and overwrite is disabled
    if os.path.exists(destination_file) and not overwrite_output:
        # end message - script
        logger.info(f" ----> Skip existing destination file: {destination_file}")
        logger.info(" ----> Dumping execution ... SKIPPED. DATA PREVIOUSLY COMPUTED")
        return False
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - write datasets
    logger.info(" ----> Write datasets ...")

    # output info
    logger.info(f" -----> Soil moisture map shape: {soil_moisture_map_out.shape}")
    logger.info(f" -----> Time-lag map shape: {time_lag_map.shape}")

    # write datasets
    write_output_map(
        destination_file=destination_file,
        soil_moisture_map=soil_moisture_map_out,
        time_lag_map=time_lag_map,
        profile=profile,
        value_var=variable_name,
        fill_value=fill_value
    )

    # message end - write datasets
    stats["written_files"] = 1
    logger.info(f" -----> Datasets file: {destination_file}")
    logger.info(" ----> Write datasets ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # message start - write report
    logger.info(" ----> Write report ...")

    if not report_enabled:
        logger.info(" ----> Write report ... SKIPPED. NOT ACTIVE")
        return None

    # build processing report
    report_obj = collect_report(
        settings=settings,
        reference_time=time_reference,
        time_start=time_start,
        time_end=time_end,
        destination_file=destination_file,
        stats=stats,
        start_time=alg_start
    )
    save_report(report_path, report_obj)

    # message end - write report
    logger.info(f" -----> Report file: {report_path}")
    logger.info(" ----> Write report ... DONE")
    # ---------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------
    # end message - script
    logger.info(f" ----> Output file written: {destination_file}")
    logger.info(f" ----> Written files: {stats['written_files']}")
    logger.info(" ----> Dumping execution ... DONE ")

    return True
    # ---------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
