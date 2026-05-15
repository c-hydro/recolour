#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECOLOUR APPS - GRID INFO

__date__ = '20260421'
__version__ = '1.0.0'
__author__ =
    'Fabio Delogu (fabio.delogu@cimafoundation.org)'
__library__ = 'recolour'

General command line:
python app_grid_info.py
python app_grid_info.py -grid_file "TUW_WARP5_grid_info_2_3.nc" -domain_name "Italy"
python app_grid_info.py -grid_file "TUW_WARP5_grid_info_2_3.nc" -domain_name "Spain" -out_folder "/tmp"
python app_grid_info.py -grid_file "TUW_WARP5_grid_info_2_3.nc" -domain_name "France" -resolution "50m"

Notes:
- the script extracts GPIs and cell IDs falling inside a selected domain
- the domain is searched in Natural Earth admin_0 countries
- output files are saved as:
    <domain_name>_cells.txt
    <domain_name>_gpis.txt
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import os
import warnings
import argparse
import numpy as np
import xarray as xr
import geopandas as gpd

from shapely import contains_xy, intersects_xy
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader

# suppress warnings
warnings.filterwarnings("ignore")
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# algorithm information
project_name = 'recolour'
alg_name = 'Grid domain selection by Natural Earth country geometry'
alg_type = 'DataProcessing'
alg_version = '1.0.0'
alg_release = '2026-04-21'
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# default file and settings
default_grid_file = "TUW_WARP5_grid_info_2_3.nc"
default_domain_name = "Italy"
default_resolution = "10m"
default_out_folder = None
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to sanitize domain name for file naming
def sanitize_domain_name(domain_name):
    domain_string = str(domain_name).strip().lower()
    domain_string = domain_string.replace(" ", "_")
    return domain_string
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to load domain geometry from Natural Earth
def load_domain_geometry(domain_name, resolution="10m"):
    """
    Load selected domain geometry from Natural Earth admin_0 countries.
    """

    shp = shpreader.natural_earth(
        resolution=resolution,
        category="cultural",
        name="admin_0_countries"
    )

    gdf = gpd.read_file(shp)

    domain_name_cmp = domain_name.strip().lower()

    field_names = [
        "NAME",
        "NAME_LONG",
        "ADMIN",
        "SOVEREIGNT",
        "FORMAL_EN",
        "BRK_NAME"
    ]

    domain_mask = np.zeros(len(gdf), dtype=bool)

    for field_name in field_names:
        if field_name in gdf.columns:
            field_values = gdf[field_name].fillna("").astype(str).str.strip().str.lower()
            domain_mask = domain_mask | (field_values == domain_name_cmp)

    domain_gdf = gdf[domain_mask]

    if domain_gdf.empty:
        raise RuntimeError(
            f'Domain "{domain_name}" not found in Natural Earth admin_0 countries'
        )

    domain_geometry = unary_union(domain_gdf.geometry.values)

    return domain_geometry
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to read grid file
def read_grid_file(grid_file):
    """
    Read grid file and required variables.
    """

    if not os.path.exists(grid_file):
        raise FileNotFoundError(f'Grid file "{grid_file}" not found')

    grid_dset = xr.open_dataset(grid_file)

    required_vars = ["lon", "lat", "gpi", "cell"]
    missing_vars = [
        var_name for var_name in required_vars
        if (var_name not in grid_dset.variables) and (var_name not in grid_dset.coords)
    ]

    if missing_vars:
        grid_dset.close()
        raise RuntimeError(f'Missing required variables in grid file: {missing_vars}')

    grid_lons = grid_dset["lon"].values.astype(np.float64)
    grid_lats = grid_dset["lat"].values.astype(np.float64)
    grid_gpis = grid_dset["gpi"].values.astype(np.int64)
    grid_cells = grid_dset["cell"].values.astype(np.int64)

    return grid_dset, grid_lons, grid_lats, grid_gpis, grid_cells
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to select points inside domain
def select_grid_domain(grid_lons, grid_lats, grid_gpis, grid_cells, domain_geometry):
    """
    Select grid points inside selected domain geometry.
    """

    min_x, min_y, max_x, max_y = domain_geometry.bounds

    # apply bounding box filter before polygon selection
    bbox_mask = (
        (grid_lons >= min_x) & (grid_lons <= max_x) &
        (grid_lats >= min_y) & (grid_lats <= max_y)
    )

    grid_lons_bbox = grid_lons[bbox_mask]
    grid_lats_bbox = grid_lats[bbox_mask]
    grid_gpis_bbox = grid_gpis[bbox_mask]
    grid_cells_bbox = grid_cells[bbox_mask]

    # apply vectorized point-in-polygon test
    domain_mask = (
        contains_xy(domain_geometry, grid_lons_bbox, grid_lats_bbox) |
        intersects_xy(domain_geometry, grid_lons_bbox, grid_lats_bbox)
    )

    domain_gpis = grid_gpis_bbox[domain_mask]
    domain_cells = np.unique(grid_cells_bbox[domain_mask])

    return domain_gpis, domain_cells, bbox_mask
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to save outputs
def save_domain_data(out_cells, out_gpis, domain_cells, domain_gpis):

    np.savetxt(out_cells, domain_cells, fmt="%d")
    np.savetxt(out_gpis, domain_gpis, fmt="%d")
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# helper to get script arguments
def get_args():

    parser = argparse.ArgumentParser(
        description="Extract grid GPIs and cell IDs inside a selected Natural Earth domain"
    )

    parser.add_argument(
        "-grid_file",
        dest="grid_file",
        default=default_grid_file,
        help=f'Input grid file (default: "{default_grid_file}")'
    )

    parser.add_argument(
        "-domain_name",
        dest="domain_name",
        default=default_domain_name,
        help=f'Domain name in Natural Earth (default: "{default_domain_name}")'
    )

    parser.add_argument(
        "-resolution",
        dest="resolution",
        default=default_resolution,
        choices=["10m", "50m", "110m"],
        help='Natural Earth resolution: "10m", "50m", "110m"'
    )

    parser.add_argument(
        "-out_folder",
        dest="out_folder",
        default=default_out_folder,
        help="Output folder (default: script working folder)"
    )

    parser.add_argument(
        "-out_cells",
        dest="out_cells",
        default=None,
        help="Output txt file for selected cell IDs"
    )

    parser.add_argument(
        "-out_gpis",
        dest="out_gpis",
        default=None,
        help="Output txt file for selected GPI IDs"
    )

    return parser.parse_args()
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# main
def main():

    # ------------------------------------------------------------------------------------------------------------------
    # get args
    args = get_args()

    grid_file = os.path.abspath(args.grid_file)
    domain_name = args.domain_name
    resolution = args.resolution

    if args.out_folder is not None:
        out_folder = os.path.abspath(args.out_folder)
    else:
        out_folder = os.getcwd()

    if not os.path.exists(out_folder):
        os.makedirs(out_folder, exist_ok=True)

    domain_tag = sanitize_domain_name(domain_name)

    if args.out_cells is not None:
        out_cells = os.path.abspath(args.out_cells)
    else:
        out_cells = os.path.join(out_folder, f"{domain_tag}_cells.txt")

    if args.out_gpis is not None:
        out_gpis = os.path.abspath(args.out_gpis)
    else:
        out_gpis = os.path.join(out_folder, f"{domain_tag}_gpis.txt")
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # start info
    print(' ============================================================================ ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> START ... ')
    print(' ')
    print(f' ---> Grid file:          {grid_file}')
    print(f' ---> Domain name:        {domain_name}')
    print(f' ---> Resolution:         {resolution}')
    print(f' ---> Output folder:      {out_folder}')
    print(f' ---> Output cells file:  {out_cells}')
    print(f' ---> Output gpis file:   {out_gpis}')
    print(' ')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # read grid data
    grid_dset, grid_lons, grid_lats, grid_gpis, grid_cells = read_grid_file(grid_file)

    print(f' ---> Grid points:        {grid_lons.size}')

    # load domain geometry
    domain_geometry = load_domain_geometry(domain_name, resolution=resolution)
    min_x, min_y, max_x, max_y = domain_geometry.bounds

    print(f' ---> Domain bbox:        lon[{min_x:.3f}, {max_x:.3f}] lat[{min_y:.3f}, {max_y:.3f}]')

    # select grid points inside domain
    domain_gpis, domain_cells, bbox_mask = select_grid_domain(
        grid_lons=grid_lons,
        grid_lats=grid_lats,
        grid_gpis=grid_gpis,
        grid_cells=grid_cells,
        domain_geometry=domain_geometry
    )

    print(f' ---> Points in bbox:     {int(np.sum(bbox_mask))}')
    print(f' ---> Domain GPIs:        {domain_gpis.size}')
    print(f' ---> Domain cells:       {domain_cells.size}')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # save outputs
    save_domain_data(
        out_cells=out_cells,
        out_gpis=out_gpis,
        domain_cells=domain_cells,
        domain_gpis=domain_gpis
    )

    print(f' ---> Save cells file ... DONE')
    print(f' ---> Save gpis file  ... DONE')
    print(' ')
    print(f' ---> First 20 cell IDs:  {domain_cells[:20]}')
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # close dataset
    grid_dset.close()
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # end message
    print(' ')
    print(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    print(' ==> ... END')
    print(' ==> Bye, Bye')
    print(' ============================================================================ ')
    # ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# call script
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------------------------------------------------