# Triple Collocation Algorithm

Triple collocation is a statistical approach that compares three independent estimates of the same variable (for example, satellite retrievals from SMAP and ASCAT, and model output from GLDAS) to quantify the error variance of each data source without requiring a perfect reference. By analyzing the pairwise cross-covariances among the three datasets, it derives unbiased estimates of each product’s random errors, helping to assess reliability and guide calibration.

## Input

- **Cell Grid Sources** (from `datasets.static.source.cell_grid`):
  - `ref`: reference dataset (e.g., HMC) for collocation.
  - `k1`: first other dataset (e.g., ECMWF).
  - `k2`: second other datasets (e.g., SMAP).
- **Geographic Grid** (from `datasets.static.source.geo_grid`):
  - spatial reference grid (GeoTIFF).

## Dynamic Datasets

- **Soil Moisture Time Series** (from `datasets.dynamic.source.soil_moisture`):
  - `ref`: dynamic cell-level SM for reference.
  - `k1`: dynamic cell-level SM for dataset `k1`.
  - `k2`: dynamic cell-level SM for dataset `k2`.
  
  Each time series dataset contains variables/features such as `soil_moisture_ref`, `soil_moisture_k1`, `soil_moisture_k2`, `time`, `lon`, and `lat`.

- **Weights** (from `datasets.dynamic.source.weights`):
  - `ref_k1`: weight array for reference vs k1 comparisons.
  - `ref_k2`: weight array for reference vs k2 comparisons.
  - `k1_k2`: weight array for k1 vs k2 comparisons.

## Ancillary

- **Points** (`datasets.dynamic.ancillary.points`):
  - cell-based ancillary point (raw, analyzed and metrics).
- **Maps** (`datasets.dynamic.ancillary.maps`):
  - cell-based ancillary maps.

## Output

- **Dynamic Output** (from `datasets.dynamic.destination`):
  - `soil_moisture_masked_resampled`: resampled and masked soil moisture map.
  - `soil_moisture_masked_filtered`: filtered soil moisture map.
  - `flags_masked`: datasets flags map.

All saved as GeoTIFF `sm_{datetime_destination}.tiff` to `/home/fabio/Desktop/Recolour_Workspace/recolour-ws/sm_tc_dynamic/v3/{sub_path_destination}`.

