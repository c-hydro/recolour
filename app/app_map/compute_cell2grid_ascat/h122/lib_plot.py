"""
Library Features:

Name:           lib_plot
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20260528'
Version:        '1.1.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to extract points
def extract_points(df, vars_list, lon_name="lon", lat_name="lat"):

    data = {}
    lons = df[lon_name].values
    lats = df[lat_name].values
    for var_name in vars_list:

        if var_name in [lon_name, lat_name]:
            continue

        data[var_name] = {
            "lon": lons,
            "lat": lats,
            "values": df[var_name].values
        }

    return data
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to plot points to map
def plot_points(grid_lons, grid_lats, grid_vals,
                src_lons, src_lats, src_vals,
                domain_mask=None,
                method=None,
                roi_km=None,
                point_styles=None):

    fig, ax = plt.subplots(figsize=(10, 10))

    extent = [
        np.nanmin(grid_lons), np.nanmax(grid_lons),
        np.nanmin(grid_lats), np.nanmax(grid_lats)
    ]

    # Domain mask
    if domain_mask is not None:

        domain_cmap = ListedColormap(["#6FA8DC", "#F5F5F5"])

        ax.imshow(
            domain_mask.astype(int),
            origin="upper",
            extent=extent,
            cmap=domain_cmap,
            vmin=0,
            vmax=1,
            alpha=0.65
        )

        grid_vals_plot = np.ma.masked_where(domain_mask == 0, grid_vals)

        domain_norm = BoundaryNorm([-0.5, 0.5, 1.5], domain_cmap.N)

        sm = ScalarMappable(cmap=domain_cmap, norm=domain_norm)
        sm.set_array([])

    else:
        grid_vals_plot = grid_vals

    # Background map
    ax.imshow(
        grid_vals_plot,
        origin="upper",
        extent=extent,
        alpha=0.65
    )

    # Default plotting configuration
    default_style = {
        "facecolor": "black",
        "edgecolor": "black",
        "size": 16,
        "linewidth": 0.3,
        "marker": "o",
        "label": None
    }

    # Optional overrides only for selected values
    # Example:
    # point_styles = {
    #     0: {"facecolor": "white"},
    #     1: {"facecolor": "red"},
    #     2: {"facecolor": "orange"}
    # }

    src_vals = np.asarray(src_vals)

    # CASE 1: custom styles defined
    if point_styles is not None:

        unique_values = np.unique(src_vals)

        for value in unique_values:

            mask = src_vals == value

            style = default_style.copy()

            # override only provided fields
            if value in point_styles:
                style.update(point_styles[value])

            ax.scatter(
                src_lons[mask],
                src_lats[mask],
                c=style["facecolor"],
                edgecolors=style["edgecolor"],
                s=style["size"],
                linewidths=style["linewidth"],
                marker=style["marker"],
                label=(
                    style["label"]
                    if style["label"] is not None
                    else f"value = {value}"
                )
            )

    # CASE 2: fallback to old behaviour
    else:

        ok_points = src_vals >= 0
        bad_points = src_vals < 0

        ax.scatter(
            src_lons[ok_points],
            src_lats[ok_points],
            c="black",
            s=default_style["size"],
            edgecolors=default_style["edgecolor"],
            linewidths=default_style["linewidth"],
            marker=default_style["marker"],
            label="points >= 0"
        )

        ax.scatter(
            src_lons[bad_points],
            src_lats[bad_points],
            c="yellow",
            s=default_style["size"],
            edgecolors=default_style["edgecolor"],
            linewidths=default_style["linewidth"],
            marker=default_style["marker"],
            label="points < 0"
        )
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")

    title = "points to map"

    if method is not None:
        title += f" | method={method}"

    if roi_km is not None:
        title += f" | roi={roi_km} km"

    ax.set_title(title)

    ax.legend(loc="best")

    plt.show(block=True)

    print()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# helper to plot map
def plot_map(
        grid_lons, grid_lats, domain_mask,
        grid_vals,
        src_lons, src_lats, src_vals,
        fill_value=-9999.0,
        method=None, roi_km=None,
        title=None,
        show=True, plot_points=True,
        save_path=None):

    grid_plot = np.asarray(grid_vals, dtype=np.float32).copy()
    grid_plot[grid_plot == fill_value] = np.nan

    fig, ax = plt.subplots(figsize=(10, 10))

    # domain: 0 = white, 1 = light gray
    domain_cmap = ListedColormap(["white", "lightgray"])

    ax.imshow(
        domain_mask.astype(int),
        origin="upper",
        extent=[
            np.nanmin(grid_lons), np.nanmax(grid_lons),
            np.nanmin(grid_lats), np.nanmax(grid_lats)
        ],
        cmap=domain_cmap,
        vmin=0,
        vmax=1,
        alpha=1.0
    )

    # interpolated map
    im = ax.imshow(
        grid_plot,
        origin="upper",
        extent=[
            np.nanmin(grid_lons), np.nanmax(grid_lons),
            np.nanmin(grid_lats), np.nanmax(grid_lats)
        ],
        alpha=0.75
    )

    plt.colorbar(im, ax=ax, label="interpolated map")

    # optional point plotting
    if plot_points:

        if (
            src_lons is not None and
            src_lats is not None and
            src_vals is not None
        ):

            src_vals = np.asarray(src_vals, dtype=np.float64)

            valid_points = (
                np.isfinite(src_vals) &
                (src_vals != fill_value) &
                (src_vals > 0)
            )

            invalid_points = ~valid_points

            ax.scatter(
                src_lons[valid_points],
                src_lats[valid_points],
                c="black",
                s=14,
                marker="o",
                label="valid source points"
            )

            ax.scatter(
                src_lons[invalid_points],
                src_lats[invalid_points],
                c="red",
                s=14,
                marker="x",
                label="invalid source points"
            )

            ax.legend(loc="best")

    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")

    if title is None:
        title = "points to interpolated map"

    if method is not None:
        title += f" | method={method}"

    if roi_km is not None:
        title += f" | roi={roi_km} km"

    ax.set_title(title)
    ax.legend(loc="best")

    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show(block=True)
    else:
        plt.close(fig)

    print()
# ----------------------------------------------------------------------------------------------------------------------
