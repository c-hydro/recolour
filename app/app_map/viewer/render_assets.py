#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate OpenLayers PNG overlays and colorbars from a selected metadata JSON file.

The script is intentionally metadata-driven: change product.filename, layer
colorbar limits, units, labels, cmap, opacity or categories in the selected metadata file;
then run launcher.py. Assets are regenerated before the web server starts.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, Normalize

CUSTOM_CMAPS = {
    # Low values blue, high values red. Good if the user wants "blue to red".
    "soil_moisture_blue_red": [
        "#08306b", "#2171b5", "#6baed6", "#c6dbef",
        "#fee0d2", "#fc9272", "#de2d26", "#7f0000",
    ],
    # Low values red/brown, high values blue. Hydrology convention: dry -> wet.
    "soil_moisture_dry_wet": [
        "#7f0000", "#d73027", "#fc8d59", "#fee090",
        "#e0f3f8", "#91bfdb", "#4575b4", "#08306b",
    ],
    "type_categorical": ["#d9d9d9", "#1f78b4", "#33a02c"],
}


def slugify(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "layer"


def get_cmap(name: str, categories: Dict[str, str] | None = None):
    if name in CUSTOM_CMAPS:
        colors = CUSTOM_CMAPS[name]
        if categories:
            keys = sorted([int(k) for k in categories.keys()])
            n = max(keys) + 1 if keys else len(colors)
            colors = (colors + colors[-1:] * max(0, n - len(colors)))[:n]
            return ListedColormap(colors, name=name)
        return LinearSegmentedColormap.from_list(name, colors)
    if name == "categorical":
        return get_cmap("type_categorical", categories)
    try:
        return plt.get_cmap(name)
    except ValueError:
        print(f"WARNING: unknown cmap '{name}', falling back to viridis")
        return plt.get_cmap("viridis")



def apply_layer_scaling(arr: np.ndarray, layer: Dict[str, Any]) -> Tuple[np.ndarray, Dict[str, Any] | None]:
    """Return data after optional metadata-driven scaling.

    Supported layer keys:
    - scale_factor: numeric multiplier applied directly (for example 100 for 0-1 -> 0-100).
    - scale_mode: "fraction_to_percent_if_needed" converts 0-1 style soil-moisture
      values to 0-100 only when valid data are <= 1.0. Already-percent 0-100 HMC
      files are left untouched.
    - scale_mode: "percent_to_fraction_if_needed" is kept for backwards compatibility.
    """
    out = arr.astype("float64", copy=True)
    scale_info: Dict[str, Any] | None = None

    if "scale_factor" in layer:
        factor = float(layer.get("scale_factor", 1.0))
        out *= factor
        scale_info = {"mode": "scale_factor", "factor": factor}

    mode = str(layer.get("scale_mode", "")).lower().strip()
    finite = out[np.isfinite(out)]
    if mode == "fraction_to_percent_if_needed":
        # HMC soil moisture viewer convention: display as 0-100.
        # Files already stored in 0-100 require no transformation.
        if finite.size and float(np.nanmax(finite)) <= 1.0:
            out *= 100.0
            scale_info = {"mode": mode, "factor": 100.0, "applied": True}
        else:
            scale_info = {"mode": mode, "factor": 1.0, "applied": False}
    elif mode == "percent_to_fraction_if_needed":
        # Backwards-compatible mode for metadata that wants 0-1 output.
        if finite.size and float(np.nanmax(finite)) > 1.0:
            out *= 0.01
            scale_info = {"mode": mode, "factor": 0.01, "applied": True}
        else:
            scale_info = {"mode": mode, "factor": 1.0, "applied": False}

    return out, scale_info

def rgba_image(data: np.ndarray, nodata: float | None, vmin: float, vmax: float,
               cmap_name: str, categories: Dict[str, str] | None,
               mask_outside: bool) -> np.ndarray:
    arr = data.astype("float64")
    mask = ~np.isfinite(arr)
    if nodata is not None and np.isfinite(nodata):
        mask |= arr == float(nodata)

    if categories:
        cmap = get_cmap(cmap_name, categories)
        keys = sorted([int(k) for k in categories.keys()])
        min_key, max_key = min(keys), max(keys)
        bounds = np.arange(min_key - 0.5, max_key + 1.5, 1)
        norm = BoundaryNorm(bounds, cmap.N)
        rgba = cmap(norm(arr))
        mask |= (arr < min_key) | (arr > max_key)
    else:
        cmap = get_cmap(cmap_name)
        if vmax <= vmin:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
        rgba = cmap(norm(arr))
        if mask_outside:
            mask |= (arr < vmin) | (arr > vmax)

    rgba[..., 3] = np.where(mask, 0.0, rgba[..., 3])
    return (rgba * 255).astype("uint8")


def save_png(rgba: np.ndarray, out_file: Path) -> None:
    import PIL.Image
    out_file.parent.mkdir(parents=True, exist_ok=True)
    PIL.Image.fromarray(rgba, mode="RGBA").save(out_file)


def save_colorbar(layer: Dict[str, Any], out_file: Path) -> None:
    colorbar = layer["colorbar"]
    vmin = float(colorbar.get("min", 0))
    vmax = float(colorbar.get("max", 1))
    cmap_name = colorbar.get("cmap", "viridis")
    categories = layer.get("categories")
    label = colorbar.get("label", layer.get("name", ""))

    out_file.parent.mkdir(parents=True, exist_ok=True)

    if categories:
        keys = sorted([int(k) for k in categories.keys()])
        cmap = get_cmap(cmap_name, categories)
        fig, ax = plt.subplots(figsize=(5.6, 0.75), dpi=140)
        bounds = np.arange(min(keys) - 0.5, max(keys) + 1.5, 1)
        norm = BoundaryNorm(bounds, cmap.N)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=ax, orientation="horizontal", ticks=keys)
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=8)
    else:
        cmap = get_cmap(cmap_name)
        fig, ax = plt.subplots(figsize=(5.6, 0.75), dpi=140)
        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=ax, orientation="horizontal")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=8)

    fig.tight_layout(pad=0.15)
    fig.savefig(out_file, transparent=False, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def stats(arr: np.ndarray, nodata: float | None) -> Dict[str, Any]:
    a = arr.astype("float64")
    mask = np.isfinite(a)
    if nodata is not None and np.isfinite(nodata):
        mask &= a != float(nodata)
    valid = a[mask]
    if valid.size == 0:
        return {"min": None, "max": None, "mean": None, "valid_pixels": 0}
    return {
        "min": float(np.nanmin(valid)),
        "max": float(np.nanmax(valid)),
        "mean": float(np.nanmean(valid)),
        "valid_pixels": int(valid.size),
    }


def regenerate(metadata_file: Path, write_metadata: bool = True, raster_override: str | None = None) -> Dict[str, Any]:
    with metadata_file.open("r", encoding="utf-8") as f:
        metadata = json.load(f)

    root = metadata_file.parent
    product = metadata.setdefault("product", {})

    if raster_override:
        tif_path = Path(raster_override).expanduser()
        if not tif_path.is_absolute():
            tif_path = (Path.cwd() / tif_path).resolve()
        else:
            tif_path = tif_path.resolve()
        product["filename"] = str(tif_path)
        product["source_path"] = str(tif_path)
        product["source_basename"] = tif_path.name
    else:
        filename = product.get("filename") or metadata.get("source_file")
        if not filename:
            raise RuntimeError("metadata.product.filename is missing")
        tif_path = Path(filename).expanduser()
        if not tif_path.is_absolute():
            tif_path = root / tif_path
        tif_path = tif_path.resolve()
        product["source_path"] = str(tif_path)
        product["source_basename"] = tif_path.name

    if not tif_path.exists():
        raise FileNotFoundError(f"GeoTIFF not found: {tif_path}")

    with rasterio.open(tif_path) as src:
        nodata = src.nodata if src.nodata is not None else metadata.get("nodata")
        metadata["width"] = int(src.width)
        metadata["height"] = int(src.height)
        metadata["bands"] = int(src.count)
        # CRS and extent handling
        # IMPORTANT: do not call rasterio.warp.transform_bounds here. On some
        # lightweight Python/GDAL installations it raises CRSError even for
        # EPSG:4326 rasters. This viewer uses EPSG:4326 ImageStatic overlays,
        # and the distributed ASCAT file is already in EPSG:4326, so native
        # GeoTIFF bounds are the safest default. If you use a raster in another
        # CRS, provide extent_epsg4326 manually in metadata.json before launch.
        src_crs_text = str(src.crs) if src.crs else metadata.get("crs", "unknown")
        metadata["crs"] = src_crs_text
        metadata["nodata"] = nodata

        src_bounds = [float(src.bounds.left), float(src.bounds.bottom),
                      float(src.bounds.right), float(src.bounds.top)]
        meta_crs_text = str(metadata.get("crs", "")).upper()
        product_crs_text = str(metadata.get("product", {}).get("crs", "")).upper()
        crs_text_joined = " ".join([src_crs_text.upper(), meta_crs_text, product_crs_text])

        if "4326" in crs_text_joined or "WGS 84" in crs_text_joined:
            metadata["extent_epsg4326"] = src_bounds
        elif metadata.get("extent_epsg4326"):
            print("WARNING: raster CRS is not EPSG:4326; using existing extent_epsg4326 from metadata.json")
        else:
            print("WARNING: raster CRS is not EPSG:4326 and no extent_epsg4326 is available.")
            print("WARNING: using native bounds; verify georeferencing in the browser.")
            metadata["extent_epsg4326"] = src_bounds

        product_name = metadata.get("product", {}).get("generic_name") or metadata.get("product", {}).get("name", "raster")

        for layer in metadata.get("layers", []):
            band = int(layer["band"])
            if band < 1 or band > src.count:
                print(f"WARNING: skipping band {band}; GeoTIFF has {src.count} band(s)")
                continue

            raw_arr = src.read(band)
            arr, scale_info = apply_layer_scaling(raw_arr, layer)
            if scale_info is not None:
                layer["scale_applied"] = scale_info
            cb = layer.setdefault("colorbar", {})
            vmin = float(cb.get("min", np.nanmin(arr)))
            vmax = float(cb.get("max", np.nanmax(arr)))
            cmap_name = cb.get("cmap", "viridis")
            mask_outside = bool(layer.get("mask_values_outside_colorbar", False))
            categories = layer.get("categories")

            layer_id = layer.get("id") or slugify(layer.get("name", f"band_{band}"))
            layer["id"] = layer_id
            image_name = layer.get("image") or f"{product_name}_band_{band}_{layer_id}.png"
            cb_name = cb.get("image") or f"colorbar_band_{band}_{layer_id}.png"
            layer["image"] = image_name
            cb["image"] = cb_name

            rgba = rgba_image(arr, nodata, vmin, vmax, cmap_name, categories, mask_outside)
            save_png(rgba, root / image_name)
            save_colorbar(layer, root / cb_name)
            layer["statistics"] = stats(arr, nodata)

    metadata["generated_at_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    if write_metadata:
        with metadata_file.open("w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
            f.write("\n")
    return metadata


def main() -> int:
    print("render_assets.py v12.0-categorical-export-legend")
    parser = argparse.ArgumentParser(description="Regenerate layer PNGs and colorbars from a metadata JSON file.")
    parser.add_argument("--metadata", default="metadata.json", help="Path to metadata.json")
    parser.add_argument("--raster", default=None, help="Override metadata product.filename with this GeoTIFF path. Absolute paths are supported.")
    parser.add_argument("--no-write-metadata", action="store_true", help="Do not update metadata statistics/timestamp")
    args = parser.parse_args()
    metadata_file = Path(args.metadata).resolve()
    metadata = regenerate(metadata_file, write_metadata=not args.no_write_metadata, raster_override=args.raster)
    print(f"Regenerated assets from: {metadata['product'].get('source_path') or metadata['product'].get('filename')}")
    print(f"Generated at UTC     : {metadata.get('generated_at_utc')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
