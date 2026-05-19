#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

from rasterio.mask import mask
from shapely.geometry import mapping


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        type=str,
        required=True
    )

    return parser.parse_args()


def clean_values(values, nodata):
    values = values.astype(float).ravel()

    if nodata is not None:
        values = values[values != nodata]

    values = values[np.isfinite(values)]

    return values


def extract_nearest_value(raster_path, lon, lat):
    with rasterio.open(raster_path) as src:

        value = list(src.sample([(lon, lat)]))[0][0]

        if src.nodata is not None and value == src.nodata:
            return None

        if not np.isfinite(value):
            return None

        return float(value)


def extract_buffer_stats(raster_path, point_geom,
                         radius_m, buffer_epsg):

    with rasterio.open(raster_path) as src:

        point_gdf = gpd.GeoDataFrame(
            geometry=[point_geom],
            crs="EPSG:4326"
        )

        buffer_gdf = gpd.GeoDataFrame(
            geometry=point_gdf.to_crs(
                f"EPSG:{buffer_epsg}"
            ).buffer(radius_m),
            crs=f"EPSG:{buffer_epsg}"
        ).to_crs(src.crs)

        try:

            array, _ = mask(
                src,
                [mapping(buffer_gdf.geometry.iloc[0])],
                crop=True
            )

        except Exception:

            return {
                "count": 0,
                "mean": None,
                "median": None,
                "std": None,
                "min": None,
                "max": None
            }

        values = clean_values(array[0], src.nodata)

        if values.size == 0:

            return {
                "count": 0,
                "mean": None,
                "median": None,
                "std": None,
                "min": None,
                "max": None
            }

        return {
            "count": int(values.size),
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values))
        }


def compute_porosity_from_bdod(
        bdod,
        particle_density):

    if bdod is None:
        return None

    bulk_density = bdod / 100.0

    porosity = (
            1.0 -
            bulk_density / particle_density
    )

    return float(porosity)


def read_points(config):

    points_cfg = config["points"]

    points_file = (
            Path(points_cfg["folder_name"]) /
            points_cfg["file_name"]
    )

    dataframe = pd.read_csv(
        points_file,
        skipinitialspace=True
    )

    dataframe.columns = (
        dataframe.columns.str.strip()
    )

    if points_cfg.get("valid_column") is not None:

        dataframe = dataframe[
            dataframe[
                points_cfg["valid_column"]
            ] == points_cfg["valid_value"]
            ].copy()

    geodataframe = gpd.GeoDataFrame(
        dataframe,
        geometry=gpd.points_from_xy(
            dataframe[
                points_cfg["longitude_column"]
            ],
            dataframe[
                points_cfg["latitude_column"]
            ]
        ),
        crs="EPSG:4326"
    )

    return geodataframe


def get_raster_path(
        config,
        variable,
        depth,
        statistic):

    raster_folder = Path(
        config["soilgrids"]["folder_name"]
    )

    template = config["soilgrids"]["file_name"]

    file_name = template.format(
        variable=variable,
        depth=depth,
        statistic=statistic
    )

    raster_path = (
            raster_folder / file_name
    )

    return raster_path


def extract_points(config):

    points_cfg = config["points"]
    extraction_cfg = config["extraction"]

    variables = (
        config["soilgrids"]["variables"]
    )

    depths = (
        config["soilgrids"]["depths"]
    )

    statistic = (
        config["soilgrids"]["statistic"]
    )

    use_buffer = extraction_cfg.get(
        "use_buffer", True
    )

    buffer_radius_m = extraction_cfg.get(
        "buffer_radius_m", 500
    )

    buffer_epsg = extraction_cfg.get(
        "buffer_epsg", 32632
    )

    aggregation = extraction_cfg.get(
        "aggregation", "median"
    )

    particle_density = extraction_cfg.get(
        "particle_density_g_cm3",
        2.65
    )

    results = []

    geodataframe = read_points(config)

    for _, row in geodataframe.iterrows():

        station = {
            "id": int(
                row[
                    points_cfg["id_column"]
                ]
            ),

            "station_name": str(
                row[
                    points_cfg["name_column"]
                ]
            ),

            "station_code": str(
                row[
                    points_cfg["code_column"]
                ]
            ),

            "longitude": float(
                row[
                    points_cfg["longitude_column"]
                ]
            ),

            "latitude": float(
                row[
                    points_cfg["latitude_column"]
                ]
            ),

            "depths": {}
        }

        for depth in depths:

            station["depths"][depth] = {}

            for variable in variables:

                raster_path = get_raster_path(
                    config,
                    variable,
                    depth,
                    statistic
                )

                if not raster_path.exists():

                    station["depths"][depth][
                        variable
                    ] = {
                        "error":
                            f"Missing raster: "
                            f"{raster_path}"
                    }

                    continue

                if use_buffer:

                    extracted = (
                        extract_buffer_stats(
                            raster_path,
                            row.geometry,
                            buffer_radius_m,
                            buffer_epsg
                        )
                    )

                else:

                    extracted = {
                        "value":
                            extract_nearest_value(
                                raster_path,
                                row[
                                    points_cfg[
                                        "longitude_column"
                                    ]
                                ],
                                row[
                                    points_cfg[
                                        "latitude_column"
                                    ]
                                ]
                            )
                    }

                station["depths"][depth][
                    variable
                ] = extracted

            bdod_data = (
                station["depths"][depth]
                .get("bdod", {})
            )

            if use_buffer:

                bdod_reference = (
                    bdod_data.get(aggregation)
                )

            else:

                bdod_reference = (
                    bdod_data.get("value")
                )

            station["depths"][depth][
                "porosity"
            ] = {

                "value":
                    compute_porosity_from_bdod(
                        bdod_reference,
                        particle_density
                    ),

                "units": "fraction",

                "source": "bdod",

                "source_statistic":
                    aggregation
                    if use_buffer
                    else "nearest",

                "formula":
                    "1 - ((bdod / 100)"
                    " / particle_density)"
            }

        results.append(station)

    return results


def main():

    args = get_args()

    with open(
            args.settings_file,
            "r") as file_handle:

        config = json.load(file_handle)

    output_dir = Path(
        config["output"]["folder_name"]
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = (
            output_dir /
            config["output"]["json_file_name"]
    )

    results = extract_points(config)

    output = {

        "metadata": {

            "source": "SoilGrids",

            "points_file": str(
                Path(
                    config["points"]["folder_name"]
                ) /
                config["points"]["file_name"]
            ),

            "soilgrids_folder": str(
                Path(
                    config["soilgrids"][
                        "folder_name"
                    ]
                )
            ),

            "soilgrids_template":
                config["soilgrids"][
                    "file_name"
                ],

            "variables":
                config["soilgrids"][
                    "variables"
                ],

            "depths":
                config["soilgrids"][
                    "depths"
                ],

            "statistic":
                config["soilgrids"][
                    "statistic"
                ],

            "extraction":
                config["extraction"]
        },

        "stations": results
    }

    with open(
            output_file,
            "w") as file_handle:

        json.dump(
            output,
            file_handle,
            indent=2
        )

    print(
        f"Saved: {output_file}"
    )


if __name__ == "__main__":
    main()