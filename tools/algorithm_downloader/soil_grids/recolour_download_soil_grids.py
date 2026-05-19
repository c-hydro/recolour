#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from owslib.wcs import WebCoverageService


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-settings_file",
        dest="settings_file",
        type=str,
        required=True
    )
    return parser.parse_args()


def read_ascii_grid_header(file_path):
    header = {}

    with open(file_path, "r") as file_handle:
        for _ in range(6):
            key, value = file_handle.readline().split()[:2]
            header[key.lower()] = float(value)

    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    xll = header["xllcorner"]
    yll = header["yllcorner"]
    cellsize = header["cellsize"]

    bbox = {
        "xmin": xll,
        "ymin": yll,
        "xmax": xll + ncols * cellsize,
        "ymax": yll + nrows * cellsize
    }

    return header, bbox


def apply_bbox_buffer(bbox, buffer_deg):
    return (
        bbox["xmin"] - buffer_deg,
        bbox["ymin"] - buffer_deg,
        bbox["xmax"] + buffer_deg,
        bbox["ymax"] + buffer_deg
    )


def is_geotiff(content):
    return content[:4] in [b"II*\x00", b"MM\x00*"]


def download_soilgrids(config, bbox):
    raster_dir = Path(config["soilgrids"]["folder_name"])
    raster_dir.mkdir(parents=True, exist_ok=True)

    variables = config["soilgrids"]["variables"]
    depths = config["soilgrids"]["depths"]
    statistic = config["soilgrids"]["statistic"]

    width = int(config["soilgrids"].get("width", 1024))
    height = int(config["soilgrids"].get("height", 1024))
    wcs_version = config["soilgrids"].get("wcs_version", "1.0.0")

    output_format = config["soilgrids"].get(
        "format",
        "GEOTIFF_INT16"
    )

    for variable in variables:
        service_url = f"https://maps.isric.org/mapserv?map=/map/{variable}.map"

        print(f"Connecting to {service_url}")

        wcs = WebCoverageService(
            service_url,
            version=wcs_version
        )

        for depth in depths:
            coverage_id = f"{variable}_{depth}_{statistic}"
            output_file = raster_dir / f"{coverage_id}.tif"

            if output_file.exists():
                print(f"Already exists: {output_file}")
                continue

            if coverage_id not in wcs.contents:
                print(f"Coverage not found: {coverage_id}")
                continue

            print(f"Downloading {coverage_id}")

            response = wcs.getCoverage(
                identifier=coverage_id,
                bbox=bbox,
                crs="EPSG:4326",
                format=output_format,
                width=width,
                height=height
            )

            content = response.read()

            if not is_geotiff(content):
                preview = content[:1000].decode(
                    "utf-8",
                    errors="ignore"
                )

                print("\nServer response is not a GeoTIFF.")
                print("Response preview:")
                print(preview)

                raise RuntimeError(
                    f"Downloaded file is not a GeoTIFF: {coverage_id}"
                )

            with open(output_file, "wb") as file_handle:
                file_handle.write(content)

            print(f"Saved: {output_file}")


def main():
    args = get_args()

    with open(args.settings_file, "r") as file_handle:
        config = json.load(file_handle)

    domain_file = (
        Path(config["domain"]["folder_name"]) /
        config["domain"]["file_name"]
    )

    _, bbox_raw = read_ascii_grid_header(domain_file)

    bbox = apply_bbox_buffer(
        bbox_raw,
        config["domain"].get("bbox_buffer_deg", 0.0)
    )

    print("Domain bbox:", bbox_raw)
    print("Download bbox:", bbox)

    download_soilgrids(config, bbox)


if __name__ == "__main__":
    main()