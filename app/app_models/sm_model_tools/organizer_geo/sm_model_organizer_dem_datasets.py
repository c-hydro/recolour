#!/usr/bin/env python3
import argparse
import math
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import rasterio
from rasterio.merge import merge


def read_ascii(path):
    header = {}

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for _ in range(6):
            key, val = f.readline().strip().split()[:2]
            header[key.lower()] = float(val)

        data = np.loadtxt(f)

    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    data = data.reshape((nrows, ncols))

    return header, data


def tile_id(lat, lon, product):
    res_code = "10" if product == "GLO30" else "30"
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"

    return "Copernicus_DSM_COG_%s_%s%02d_00_%s%03d_00_DEM" % (
        res_code, ns, abs(int(lat)), ew, abs(int(lon))
    )


def tile_url(tile, product):
    bucket = "copernicus-dem-30m" if product == "GLO30" else "copernicus-dem-90m"
    return "https://%s.s3.amazonaws.com/%s/%s.tif" % (bucket, tile, tile)


def download_tiles(west, south, east, north, product, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    files = []

    for lat in range(math.floor(south), math.ceil(north)):
        for lon in range(math.floor(west), math.ceil(east)):
            tile = tile_id(lat, lon, product)
            url = tile_url(tile, product)
            out = out_dir / ("%s.tif" % tile)

            if not out.exists():
                print("Downloading %s" % url)
                try:
                    urlretrieve(url, out)
                except Exception as exc:
                    print("WARNING: could not download %s: %s" % (url, exc))
                    if out.exists():
                        out.unlink()
                    continue
            else:
                print("Using existing tile %s" % out)

            files.append(out)

    if not files:
        raise RuntimeError("No Copernicus DEM tiles found.")

    return files


def bilinear_resample(source, src_transform, west, north, cellsize,
                      nrows, ncols, nodata, fill_nodata):
    dst = np.full((nrows, ncols), fill_nodata, dtype=np.float32)
    inv = ~src_transform

    xs = west + (np.arange(ncols) + 0.5) * cellsize
    ys = north - (np.arange(nrows) + 0.5) * cellsize

    for r, y in enumerate(ys):
        cols_f, rows_f = inv * (xs, np.full_like(xs, y))

        c0 = np.floor(cols_f).astype(int)
        r0 = np.floor(rows_f).astype(int)

        dc = cols_f - c0
        dr = rows_f - r0

        valid = (
            (r0 >= 0)
            & (r0 < source.shape[0] - 1)
            & (c0 >= 0)
            & (c0 < source.shape[1] - 1)
        )

        if not np.any(valid):
            continue

        rr = r0[valid]
        cc = c0[valid]

        q11 = source[rr, cc]
        q21 = source[rr, cc + 1]
        q12 = source[rr + 1, cc]
        q22 = source[rr + 1, cc + 1]

        ddx = dc[valid]
        ddy = dr[valid]

        vals = (
            q11 * (1 - ddx) * (1 - ddy)
            + q21 * ddx * (1 - ddy)
            + q12 * (1 - ddx) * ddy
            + q22 * ddx * ddy
        )

        bad = ~np.isfinite(vals)

        if nodata is not None:
            bad |= (
                (q11 == nodata)
                | (q21 == nodata)
                | (q12 == nodata)
                | (q22 == nodata)
            )

        vals[bad] = fill_nodata
        dst[r, np.where(valid)[0]] = vals.astype(np.float32)

    return dst


def write_ascii(path, arr, ncols, nrows, xll, yll, cellsize, nodata):
    arr = np.asarray(arr, dtype=np.float32)
    arr = np.where(np.isfinite(arr), arr, nodata)

    with open(path, "w", encoding="utf-8") as f:
        f.write("ncols        %d\n" % ncols)
        f.write("nrows        %d\n" % nrows)
        f.write("xllcorner    %.12f\n" % xll)
        f.write("yllcorner    %.12f\n" % yll)
        f.write("cellsize     %.12f\n" % cellsize)

        if float(nodata).is_integer():
            f.write("NODATA_value  %d\n" % int(nodata))
        else:
            f.write("NODATA_value  %s\n" % nodata)

        np.savetxt(f, arr, fmt="%.2f")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--template", default="marche.dem.txt")
    parser.add_argument("--output", default="marche_copernicus_glo30_land_only.asc")
    parser.add_argument("--product", choices=["GLO30", "GLO90"], default="GLO30")
    parser.add_argument("--tiles-dir", default="copernicus_dem_tiles")
    parser.add_argument("--resampling", choices=["bilinear"], default="bilinear")
    parser.add_argument("--fill-nodata", type=float, default=0.0)

    parser.add_argument(
        "--mask-sea",
        action="store_true",
        default=True,
        help="Mask Copernicus sea cells using elevation threshold. Default: true.",
    )
    parser.add_argument(
        "--no-mask-sea",
        action="store_false",
        dest="mask_sea",
        help="Do not mask sea cells.",
    )
    parser.add_argument(
        "--sea-level",
        type=float,
        default=0.0,
        help="Cells <= this elevation are set to NoData. Default: 0.0",
    )
    parser.add_argument(
        "--use-template-mask",
        action="store_true",
        default=False,
        help="Also apply original/template DEM mask. Default: false.",
    )

    args = parser.parse_args()

    h, template_dem = read_ascii(args.template)

    ncols = int(h["ncols"])
    nrows = int(h["nrows"])
    west = h["xllcorner"]
    south = h["yllcorner"]
    cellsize = h["cellsize"]
    nodata = h.get("nodata_value", -9999.0)

    east = west + ncols * cellsize
    north = south + nrows * cellsize

    print("Target grid:")
    print("  west=%.12f east=%.12f" % (west, east))
    print("  south=%.12f north=%.12f" % (south, north))
    print("  ncols=%d nrows=%d cellsize=%.12f" % (ncols, nrows, cellsize))
    print("  mask sea: %s, sea level threshold: %.2f" % (args.mask_sea, args.sea_level))
    print("  use template mask: %s" % args.use_template_mask)

    tile_files = download_tiles(
        west=west,
        south=south,
        east=east,
        north=north,
        product=args.product,
        out_dir=Path(args.tiles_dir),
    )

    print("Mosaicking %d tiles" % len(tile_files))

    srcs = [rasterio.open(str(p)) for p in tile_files]

    try:
        mosaic, mosaic_transform = merge(
            srcs,
            bounds=(west, south, east, north),
        )

        source = mosaic[0].astype(np.float32)
        source_nodata = srcs[0].nodata

    finally:
        for src in srcs:
            src.close()

    dst = bilinear_resample(
        source=source,
        src_transform=mosaic_transform,
        west=west,
        north=north,
        cellsize=cellsize,
        nrows=nrows,
        ncols=ncols,
        nodata=source_nodata,
        fill_nodata=args.fill_nodata,
    )

    if args.mask_sea:
        sea_mask = np.isfinite(dst) & (dst <= args.sea_level)
        print("  sea cells masked: %d" % int(np.sum(sea_mask)))
        dst[sea_mask] = nodata

    if args.use_template_mask:
        template_mask = (
            np.isfinite(template_dem)
            & (template_dem != nodata)
            & (template_dem > 0)
        )
        print("  template-mask cells kept: %d" % int(np.sum(template_mask)))
        dst[~template_mask] = nodata

    write_ascii(
        path=args.output,
        arr=dst,
        ncols=ncols,
        nrows=nrows,
        xll=west,
        yll=south,
        cellsize=cellsize,
        nodata=nodata,
    )

    valid = dst[(dst != nodata) & np.isfinite(dst)]

    print("Saved ASCII DEM: %s" % Path(args.output).resolve())
    print("Valid terrain cells: %d" % valid.size)

    if valid.size > 0:
        print(
            "Elevation min/max: %.2f / %.2f"
            % (float(np.nanmin(valid)), float(np.nanmax(valid)))
        )


if __name__ == "__main__":
    main()
