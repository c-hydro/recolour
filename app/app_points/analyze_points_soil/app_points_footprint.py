#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create product footprints around soil-moisture stations.

Footprint definitions:
    ASCAT:
        Circular approximation around the station.
        Default radius: 17.5 km.

    SMAP:
        Circular approximation around the station.
        Default radius: 16 km.

    ECMWF:
        Native raster cell containing the station.
        An ECMWF raster path must be provided.

    HMC:
        Native raster cell containing the station.
        The uploaded HMC GeoTIFF is used.

Output:
    GeoPackage containing:
        stations
        footprint_ascat
        footprint_smap
        footprint_ecmwf
        footprint_hmc
        footprints_all
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.transform import xy
from shapely.geometry import Polygon, box


# -------------------------------------------------------------------------------------
# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------------------
# Constants

CRS_GEOGRAPHIC = "EPSG:4326"

# Suitable metric CRS for Liguria
CRS_METRIC = "EPSG:32632"

DEFAULT_ASCAT_RADIUS_M = 17_500.0
DEFAULT_SMAP_RADIUS_M = 16_000.0


# -------------------------------------------------------------------------------------
def create_station_dataframe() -> pd.DataFrame:
    """Create the station registry."""

    return pd.DataFrame(
        [
            {
                "station_id": 3,
                "station_name": "Amborzasco",
                "station_tag": "AMBOR",
                "longitude": 9.45496,
                "latitude": 44.51447,
            },
            {
                "station_id": 4,
                "station_name": "Colle D'Oggia",
                "station_tag": "CODOG",
                "longitude": 7.86667,
                "latitude": 43.98131,
            },
            {
                "station_id": 5,
                "station_name": "Cuccarello",
                "station_tag": "CUCCA",
                "longitude": 9.69908,
                "latitude": 44.34967,
            },
            {
                "station_id": 6,
                "station_name": "Albenga - Isolabella",
                "station_tag": "ISBLL",
                "longitude": 8.17956,
                "latitude": 44.06875,
            },
            {
                "station_id": 7,
                "station_name": "Loco Carchelli",
                "station_tag": "LOCOC",
                "longitude": 9.28421,
                "latitude": 44.55421,
            },
            {
                "station_id": 8,
                "station_name": "Mignanego",
                "station_tag": "MIGNA",
                "longitude": 8.93816,
                "latitude": 44.54028,
            },
            {
                "station_id": 9,
                "station_name": "Ognio",
                "station_tag": "OGNIO",
                "longitude": 9.16994,
                "latitude": 44.44372,
            },
            {
                "station_id": 10,
                "station_name": "Urbe - Vara Sup.",
                "station_tag": "URVAS",
                "longitude": 8.62739,
                "latitude": 44.46953,
            },
            {
                "station_id": 11,
                "station_name": "Valzemola",
                "station_tag": "VALZE",
                "longitude": 8.19105,
                "latitude": 44.36959,
            },
        ]
    )


# -------------------------------------------------------------------------------------
def create_station_points(
        stations_df: pd.DataFrame,
) -> gpd.GeoDataFrame:
    """Convert the station table to geographic point geometries."""

    mandatory_columns = [
        "station_id",
        "station_name",
        "station_tag",
        "longitude",
        "latitude",
    ]

    missing_columns = [
        column_name
        for column_name in mandatory_columns
        if column_name not in stations_df.columns
    ]

    if missing_columns:
        raise KeyError(
            f"Missing station columns: {missing_columns}"
        )

    stations_gdf = gpd.GeoDataFrame(
        stations_df.copy(),
        geometry=gpd.points_from_xy(
            stations_df["longitude"],
            stations_df["latitude"],
        ),
        crs=CRS_GEOGRAPHIC,
    )

    return stations_gdf


# -------------------------------------------------------------------------------------
def create_circular_footprints(
        stations_gdf: gpd.GeoDataFrame,
        product_name: str,
        radius_m: float,
        metric_crs: str = CRS_METRIC,
) -> gpd.GeoDataFrame:
    """
    Create circular station-centred footprints.

    Buffering is performed in a metric CRS to ensure the radius is expressed
    in metres.
    """

    if radius_m <= 0:
        raise ValueError("The footprint radius must be greater than zero.")

    stations_metric = stations_gdf.to_crs(metric_crs)

    footprints_gdf = stations_metric.copy()

    footprints_gdf["geometry"] = (
        footprints_gdf.geometry.buffer(radius_m)
    )

    footprints_gdf["product"] = product_name.lower()
    footprints_gdf["footprint_type"] = "station_buffer"
    footprints_gdf["radius_m"] = float(radius_m)
    footprints_gdf["diameter_m"] = float(radius_m * 2.0)
    footprints_gdf["area_km2"] = (
        footprints_gdf.geometry.area / 1_000_000.0
    )

    footprints_gdf["source_row"] = pd.NA
    footprints_gdf["source_col"] = pd.NA
    footprints_gdf["grid_center_lon"] = np.nan
    footprints_gdf["grid_center_lat"] = np.nan
    footprints_gdf["station_grid_distance_m"] = 0.0
    footprints_gdf["source_file"] = pd.NA

    return footprints_gdf.to_crs(CRS_GEOGRAPHIC)


# -------------------------------------------------------------------------------------
def get_raster_cell_polygon(
        transform,
        row: int,
        col: int,
) -> Polygon:
    """
    Create the polygon corresponding to one raster cell.

    The polygon is returned in the raster CRS.
    """

    x_left, y_top = xy(
        transform,
        row,
        col,
        offset="ul",
    )

    x_right, y_bottom = xy(
        transform,
        row,
        col,
        offset="lr",
    )

    return box(
        min(x_left, x_right),
        min(y_bottom, y_top),
        max(x_left, x_right),
        max(y_bottom, y_top),
    )


# -------------------------------------------------------------------------------------
def geodesic_distance_m(
        longitude_1: float,
        latitude_1: float,
        longitude_2: float,
        latitude_2: float,
) -> float:
    """Calculate geodesic distance between two longitude/latitude points."""

    geod = Geod(ellps="WGS84")

    _, _, distance_m = geod.inv(
        longitude_1,
        latitude_1,
        longitude_2,
        latitude_2,
    )

    return float(distance_m)


# -------------------------------------------------------------------------------------
def create_raster_cell_footprints(
        stations_gdf: gpd.GeoDataFrame,
        raster_file: str | Path,
        product_name: str,
        read_values: bool = True,
) -> gpd.GeoDataFrame:
    """
    Create the native raster cell footprint containing each station.

    This method can be used for ECMWF and HMC, provided that a representative
    raster file is available.

    Parameters
    ----------
    stations_gdf : GeoDataFrame
        Station points in any valid CRS.
    raster_file : str or Path
        Raster used to define the native grid.
    product_name : str
        Dataset name, such as 'ecmwf' or 'hmc'.
    read_values : bool
        If True, also read the raster value at each station cell.

    Returns
    -------
    GeoDataFrame
        One native raster-cell polygon for each station.
    """

    raster_file = Path(raster_file)

    if not raster_file.exists():
        raise FileNotFoundError(
            f"Raster file not found: {raster_file}"
        )

    records = []

    with rasterio.open(raster_file) as raster_obj:

        if raster_obj.crs is None:
            raise RuntimeError(
                f"Raster CRS is not defined: {raster_file}"
            )

        # Transform station points into the raster CRS.
        stations_raster = stations_gdf.to_crs(raster_obj.crs)

        raster_bounds_polygon = box(*raster_obj.bounds)

        for _, station_row in stations_raster.iterrows():

            station_geometry = station_row.geometry

            record = {
                "station_id": station_row["station_id"],
                "station_name": station_row["station_name"],
                "station_tag": station_row["station_tag"],
                "longitude": station_row["longitude"],
                "latitude": station_row["latitude"],
                "product": product_name.lower(),
                "footprint_type": "native_raster_cell",
                "radius_m": np.nan,
                "diameter_m": np.nan,
                "source_file": raster_file.name,
            }

            # Check whether the station is covered by the raster.
            if not raster_bounds_polygon.covers(station_geometry):

                logger.warning(
                    "Station '%s' is outside the %s raster extent.",
                    station_row["station_tag"],
                    product_name.upper(),
                )

                record.update(
                    {
                        "source_row": pd.NA,
                        "source_col": pd.NA,
                        "grid_center_lon": np.nan,
                        "grid_center_lat": np.nan,
                        "station_grid_distance_m": np.nan,
                        "cell_value": np.nan,
                        "cell_valid": False,
                        "area_km2": np.nan,
                        "geometry": None,
                    }
                )

                records.append(record)
                continue

            row_idx, col_idx = raster_obj.index(
                station_geometry.x,
                station_geometry.y,
            )

            cell_polygon = get_raster_cell_polygon(
                transform=raster_obj.transform,
                row=row_idx,
                col=col_idx,
            )

            center_x, center_y = xy(
                raster_obj.transform,
                row_idx,
                col_idx,
                offset="center",
            )

            # Convert the cell centre to geographic coordinates.
            center_gdf = gpd.GeoDataFrame(
                geometry=gpd.points_from_xy(
                    [center_x],
                    [center_y],
                ),
                crs=raster_obj.crs,
            ).to_crs(CRS_GEOGRAPHIC)

            center_lon = center_gdf.geometry.iloc[0].x
            center_lat = center_gdf.geometry.iloc[0].y

            station_distance_m = geodesic_distance_m(
                longitude_1=float(station_row["longitude"]),
                latitude_1=float(station_row["latitude"]),
                longitude_2=float(center_lon),
                latitude_2=float(center_lat),
            )

            cell_value = np.nan
            cell_valid = False

            if read_values:
                raster_value = raster_obj.read(
                    1,
                    window=((row_idx, row_idx + 1), (col_idx, col_idx + 1)),
                )[0, 0]

                cell_value = float(raster_value)

                cell_valid = bool(
                    np.isfinite(raster_value)
                    and (
                        raster_obj.nodata is None
                        or raster_value != raster_obj.nodata
                    )
                )

            record.update(
                {
                    "source_row": int(row_idx),
                    "source_col": int(col_idx),
                    "grid_center_lon": float(center_lon),
                    "grid_center_lat": float(center_lat),
                    "station_grid_distance_m": station_distance_m,
                    "cell_value": cell_value,
                    "cell_valid": cell_valid,
                    "geometry": cell_polygon,
                }
            )

            records.append(record)

        footprints_gdf = gpd.GeoDataFrame(
            records,
            geometry="geometry",
            crs=raster_obj.crs,
        )

    # Calculate cell area in a metric CRS.
    footprints_metric = footprints_gdf.to_crs(CRS_METRIC)

    footprints_gdf["area_km2"] = (
        footprints_metric.geometry.area / 1_000_000.0
    )

    return footprints_gdf.to_crs(CRS_GEOGRAPHIC)


# -------------------------------------------------------------------------------------
def organize_columns(
        footprints_gdf: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """Use the same output columns for every product."""

    output_columns = [
        "station_id",
        "station_name",
        "station_tag",
        "longitude",
        "latitude",
        "product",
        "footprint_type",
        "radius_m",
        "diameter_m",
        "area_km2",
        "source_row",
        "source_col",
        "grid_center_lon",
        "grid_center_lat",
        "station_grid_distance_m",
        "cell_value",
        "cell_valid",
        "source_file",
        "geometry",
    ]

    output_gdf = footprints_gdf.copy()

    for column_name in output_columns:
        if column_name not in output_gdf.columns:
            output_gdf[column_name] = pd.NA

    return output_gdf[output_columns]


# -------------------------------------------------------------------------------------
def save_footprints(
        output_file: str | Path,
        stations_gdf: gpd.GeoDataFrame,
        ascat_gdf: gpd.GeoDataFrame,
        smap_gdf: gpd.GeoDataFrame,
        hmc_gdf: gpd.GeoDataFrame,
        ecmwf_gdf: Optional[gpd.GeoDataFrame] = None,
) -> None:
    """Save stations and footprint layers to a GeoPackage."""

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if output_file.exists():
        output_file.unlink()

    stations_gdf.to_file(
        output_file,
        layer="stations",
        driver="GPKG",
    )

    ascat_gdf.to_file(
        output_file,
        layer="footprint_ascat",
        driver="GPKG",
    )

    smap_gdf.to_file(
        output_file,
        layer="footprint_smap",
        driver="GPKG",
    )

    hmc_gdf.to_file(
        output_file,
        layer="footprint_hmc",
        driver="GPKG",
    )

    footprint_collection = [
        organize_columns(ascat_gdf),
        organize_columns(smap_gdf),
        organize_columns(hmc_gdf),
    ]

    if ecmwf_gdf is not None:

        ecmwf_gdf.to_file(
            output_file,
            layer="footprint_ecmwf",
            driver="GPKG",
        )

        footprint_collection.append(
            organize_columns(ecmwf_gdf)
        )

    footprints_all = gpd.GeoDataFrame(
        pd.concat(
            footprint_collection,
            ignore_index=True,
        ),
        geometry="geometry",
        crs=CRS_GEOGRAPHIC,
    )

    footprints_all.to_file(
        output_file,
        layer="footprints_all",
        driver="GPKG",
    )

    # Non-geometrical summary table.
    summary_df = footprints_all.drop(columns="geometry")

    summary_file = output_file.with_name(
        f"{output_file.stem}_summary.csv"
    )

    summary_df.to_csv(
        summary_file,
        index=False,
    )

    logger.info("Footprints saved to: %s", output_file)
    logger.info("Summary saved to:    %s", summary_file)


# -------------------------------------------------------------------------------------
def create_product_footprints(
        hmc_raster_file: str | Path,
        output_file: str | Path,
        ecmwf_raster_file: Optional[str | Path] = None,
        ascat_radius_m: float = DEFAULT_ASCAT_RADIUS_M,
        smap_radius_m: float = DEFAULT_SMAP_RADIUS_M,
) -> dict[str, gpd.GeoDataFrame]:
    """
    Create all product footprints.

    ECMWF is optional because its actual raster/grid file was not uploaded.
    """

    stations_df = create_station_dataframe()
    stations_gdf = create_station_points(stations_df)

    logger.info("Creating ASCAT station-centred footprints ...")

    ascat_gdf = create_circular_footprints(
        stations_gdf=stations_gdf,
        product_name="ascat",
        radius_m=ascat_radius_m,
    )

    logger.info("Creating SMAP station-centred footprints ...")

    smap_gdf = create_circular_footprints(
        stations_gdf=stations_gdf,
        product_name="smap",
        radius_m=smap_radius_m,
    )

    logger.info("Creating HMC native-cell footprints ...")

    hmc_gdf = create_raster_cell_footprints(
        stations_gdf=stations_gdf,
        raster_file=hmc_raster_file,
        product_name="hmc",
        read_values=True,
    )

    ecmwf_gdf = None

    if ecmwf_raster_file is not None:

        logger.info("Creating ECMWF native-cell footprints ...")

        ecmwf_gdf = create_raster_cell_footprints(
            stations_gdf=stations_gdf,
            raster_file=ecmwf_raster_file,
            product_name="ecmwf",
            read_values=True,
        )

    else:
        logger.warning(
            "ECMWF raster not provided. "
            "The ECMWF footprint layer will not be created."
        )

    save_footprints(
        output_file=output_file,
        stations_gdf=stations_gdf,
        ascat_gdf=ascat_gdf,
        smap_gdf=smap_gdf,
        hmc_gdf=hmc_gdf,
        ecmwf_gdf=ecmwf_gdf,
    )

    return {
        "stations": stations_gdf,
        "ascat": ascat_gdf,
        "smap": smap_gdf,
        "ecmwf": ecmwf_gdf,
        "hmc": hmc_gdf,
    }


# -------------------------------------------------------------------------------------
if __name__ == "__main__":

    hmc_file = "/home/fabio/Desktop/recolour/dset/soil/footprints/src/sm_hmc_202606041200(1).tiff"

    # Set this when an ECMWF RZSM GeoTIFF/NetCDF raster has been prepared.
    ecmwf_file = None
    output_gpkg = "/home/fabio/Desktop/recolour/dset/soil/footprints/src//soil_moisture_product_footprints.gpkg"

    footprint_data = create_product_footprints(
        hmc_raster_file=hmc_file,
        ecmwf_raster_file=ecmwf_file,
        output_file=output_gpkg,
        ascat_radius_m=17_500,
        smap_radius_m=16_000,
    )

    print(footprint_data["hmc"][
        [
            "station_tag",
            "source_row",
            "source_col",
            "grid_center_lon",
            "grid_center_lat",
            "station_grid_distance_m",
            "cell_value",
            "cell_valid",
            "area_km2",
        ]
    ])