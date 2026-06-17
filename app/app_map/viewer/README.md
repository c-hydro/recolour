# OpenLayers Raster Viewer v12.0

Metadata-driven OpenLayers viewer for georeferenced multi-band rasters.

## Default launch

```bash
python3 launcher.py
```

This uses `metadata.json` and the raster defined in `metadata.product.filename`.

## Launch with another metadata file

```bash
python3 launcher.py --metadata metadata_smap.json
```

The metadata file must be inside the viewer folder so the web app can fetch it.

## Launch with a raster override

Use `--raster` to override `metadata.product.filename` directly from the command line. Absolute paths are supported.

```bash
python3 launcher.py --metadata metadata.json --raster /absolute/path/to/product.tif
```

Relative raster paths are resolved from the shell directory where you run the command.

When `--raster` is used, the launcher passes the file to `render_assets.py`; assets are regenerated and the selected metadata file is updated with:

- `product.filename`
- `product.source_path`
- `product.source_basename`
- raster size, band count, CRS, extent, statistics and generation timestamp

## Useful options

```bash
python3 launcher.py --help
python3 launcher.py --no-browser
python3 launcher.py --reinstall-deps
python3 launcher.py --no-regenerate
```

## Export

Use **Export PNG** in the web app. The PNG export includes the active map view, compact lower-center palette box, upper-left Info box, and compass.


## v12.0 bugfix

- Categorical layers such as band 5 export a compact legend with colored swatches in the PNG, not text only.
- Continuous layers keep the compact colorbar export.

## HMC single-band soil moisture

This package now includes an HMC example configuration:

```bash
python3 launcher.py --metadata metadata_hmc.json
```

The HMC configuration uses:

- `sm_hmc_202606041200.tiff`
- band 1 only
- soil moisture colorbar range `0` to `100`
- units `[%]`
- `scale_mode: fraction_to_percent_if_needed`

`scale_mode: fraction_to_percent_if_needed` multiplies rasters stored in the `0-1` range by `100` before rendering and before statistics are written. If a future HMC file is already stored as `0-100`, no scaling is applied.

To use another HMC GeoTIFF:

```bash
python3 launcher.py --metadata metadata_hmc.json --raster /absolute/path/to/hmc_soil_moisture.tif
```

