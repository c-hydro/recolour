# **Recolour**

**Recolour** is the **REprocess Package for Soil Moisture Products**, designed to process and analyze soil moisture datasets to support hydrogeological risk prevention and reduction within the Flood‑PROOFS modelling system.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [Authors](#authors)
- [License](#license)
- [Changelog](#changelog)
- [References](#references)

## Features

- Ingest and process multiple soil moisture data products
- Reprocessing routines for calibration and quality control
- Integration with Python workflows
- Visualization tools via Jupyter notebooks
- Modular design for extension and research use

## Prerequisites

Before using **Recolour**, ensure you have:

- **Operating System:** Linux (Debian/Ubuntu 64-bit recommended)
- **Python:** Version 3.7 or higher
- **QGIS:** Version 2.18 or higher (for spatial data management)
- **Additional Tools:** Jupyter Notebook, Panoply, CDO, Ncview

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/c-hydro/recolour.git
   cd recolour
   ```
2. **Set up the Python environment:**
   ```bash
   bash setup_recolour_system_conda_python_notebook.sh
   ```
3. **Activate the environment:**
   ```bash
   source activate fp_env_python
   ```
   > After activation, your prompt should display `(fp_env_python)`.

## Directory Structure

### Applications
Core **Recolour** applications organized by processing domain. Each folder contains Python scripts (`.py`) for specific workflows along with their configuration files (`.json`).
```plaintext
recolour/
├── app_cell
│   ├── compute_cell_metrics
│   │   ├── app_cell_metrics.py
│   │   └── *.json
│   ├── compute_cell_rzsm
│   │   ├── app_cell_rzsm.py
│   │   └── *.json
│   ├── compute_cell_scaling
│   │   ├── app_cell_scaling.py
│   │   └── *.json
│   ├── compute_cell_swi
│   │   ├── app_cell_swi.py
│   │   └── *.json
│   └── convert_swath2cell
│       ├── app_swath2cell_ascat.py
│       └── *.json
├── app_map
│   ├── compute_cell2grid_ascat
│   │   └── *.json
│   ├── compute_cell2grid_ecmwf
│   │   └── *.json
│   ├── compute_cell2grid_gldas
│   │   └── *.json
│   ├── compute_cell2grid_metrics
│   │   ├── app_map_cell2grid_metrics.py
│   │   └── *.json
│   ├── compute_cell2grid_smap
│   │   ├── app_map_cell2grid_smap.py
│   │   └── *.json
│   ├── create_grid_ecmwf
│   │   ├── app_map_grid_ecmwf_nrt.py
│   │   └── *.json
│   ├── create_grid_hmc
│   │   ├── app_map_grid_hmc_nrt.py
│   │   └── *.json
│   ├── create_grid_reference
│   │   ├── app_map_grid_reference.py
│   │   └── *.json
│   └── resample_grid_src2ref
│       ├── app_map_grid_src2ref.py
│       └── *.json
├── app_points
│   └── analyze_points_ascat
│       ├── app_points_ascat.py
│       └── *.json
└── app_ts
    ├── analyze_ts
    │   ├── app_sm_ts_analyze.py
    │   └── *.json
    ├── convert_grid_ecmwf2csv
    │   ├── app_sm_grid_ecmwf2csv.py
    │   └── *.json
    ├── convert_grid_hmc2csv
    │   ├── app_sm_grid_hmc2csv.py
    │   └── *.json
    └── view_ts
        ├── app_sm_ts_view.py
        └── *.json
```

### Tools
Utility modules and algorithms supporting **Recolour** workflows. These tools provide functionality for data handling, format conversion, validation, and automated workflows.
```plaintext
tools/
├── algorithm_assimilation
├── algorithm_converter
├── algorithm_downloader
├── algorithm_organizer
├── algorithm_plot_timeseries
├── algorithm_plot_validation
├── algorithm_system
├── algorithm_transfer
├── algorithm_validation_hsaf
├── algorithm_validation_tc
└── algorithm_xml
```

## Usage

**Recolour** applications are executed via their scripts in the `app_cell`, `app_map`, `app_points`, or `app_ts` folders. Use the following pattern:
```bash
python <app_folder>/<workflow_folder>/<script_name>.py <config.json> -time "YYYY-MM-DD HH:MM"
```
- `<app_folder>`: one of `app_cell`, `app_map`, `app_points`, or `app_ts`
- `<workflow_folder>`: subdirectory under the chosen `app_*` folder
- `<script_name>.py`: the application script to run
- `<config.json>`: the configuration file in the same folder as the script
- `-time`: processing date and time in `YYYY-MM-DD HH:MM` format

For example, to compute root-zone soil moisture on June 1, 2025 at 06:00 using the ECMWF DR configuration:
```bash
python app_cell/compute_cell_rzsm/app_cell_rzsm.py app_cell/compute_cell_rzsm/app_cell_rzsm_ecmwf_dr_local.json -time "2025-06-01 06:00"
```

## Examples

- **Scripts:** Check each `app_*` folder for workflow scripts and JSON configs.
- **Notebooks:** Explore the `notebook/` directory for interactive demonstrations.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch: `git checkout -b my-feature`
3. Add changes and tests (if applicable)
4. Submit a pull request against `main`

## Authors

See [AUTHORS.md](AUTHORS.md) for a complete list of contributors.

## License

This project is licensed under the [EUPL‑1.2 License](LICENSE.md).

## Changelog

All notable changes are documented in [CHANGELOG.md](CHANGELOG.md).

## References

1. CIMA Hydrology and Hydraulics [GitHub Repository](https://github.com/c-hydro)
2. Python programming language ([python.org](https://www.python.org/))
3. QGIS project ([qgis.org](https://qgis.org/))
4. Conda environment manager ([conda.io](https://conda.io/))
5. Hydrological Model Continuum codes

