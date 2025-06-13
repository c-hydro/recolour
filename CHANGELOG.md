# Changelog

## Version 2.2.0 (2025-06-03)
- **PROJECT:** operational framework grid  
  - **APPS – CELL:** convert_img2cell  
    - Fixed bugs for grid orientation (adjust longitudes and latitudes 1d arrays) and improved error messages to alert users.  
  - **APPS – CELL:** convert_swath2cell_ascat  
    - Fixed bugs in file searching based on time information.  
  - **APPS – TS:** convert_grid2csv_hmc  
    - Fixed bugs in time period handling for `PREVIOUS_MONTH` and `CURRENT_MONTH` modes.  
- **PROJECT:** operational framework grid  
  - **APPS – POINTS:** analyze_point_ascat  
    - Beta release.

## Version 2.1.0 (2025-05-27)
- **PROJECT:** operational framework grid  
  - **APPS – MAP:** compute_cell2grid_ecmwf  
    - Fixed bugs for undefined datasets (due to time selection).  
  - **APPS – MAP:** merge_grid2ref  
    - Updated code to include additional cases and fixed related bugs.  
  - **APPS – TS:** convert_grid_ecmwf2csv, convert_grid_hmc2csv, convert_grid_smap2csv  
    - Fixed various bugs.  
  - **APPS – TS:** sync_ts, analyze_ts  
    - Fixed bugs.  
- **TOOLS – DOWNLOADER HMC:**  
  - Added tools to download legacy datasets.  
- **TOOLS – ORGANIZER SYSTEM:**  
  - Added tools to organize folder trees and backup crontab.

## Version 2.0.0 (2025-04-22)
- **PROJECT:** operational framework grid  
  - **APPS – MAP:** merge_grid2ref (first release)  
  - **APPS – TS:** convert_time_step_src2csv (extended to operational mode run)  
- **TOOLS – DOWNLOADER SMAP:** smap_downloader_spl2smp_e.py  
  - Refactored all application methods and algorithms.

## Version 1.9.0 (2025-02-27)
- **PROJECT:** soil moisture triple collocation  
  - **APPS – CELL:** app_map_create_tc (bug fixes)  
- **TOOLS – DOWNLOADER SMAP:**  
  - Bug fixes in time information handling for both `spl2smp_e.py` and `spl3smp_e.py`.

## Version 1.8.0 (2024-09-30)
- **PROJECT:** time-series  
  - **APPS – TS:** convert_grid_hmc2csv, convert_time_step_src2csv, sync_ts, analyze_ts, view_ts  
    - Numerous bug fixes (delimiter issues, time selection, unavailable datasets) and feature extensions (heatmaps via seaborn, resampling).

## Version 1.7.0 (2024-07-09)
- **PROJECT:** soil moisture rescaled (obs/mod)  
  - **APPS – MAP:** create_map_smr (first release)  
- **PROJECT:** grid  
  - **APPS – TS:** resample_grid_src2ref (first release)  
- **PROJECT:** time-series  
  - **APPS – TS:** convert_grid_ecmwf2csv, convert_grid_hmc2csv, convert_grid_ascat2csv, convert_grid_smap2csv, convert_grid_gldas2csv  
    - Fixed reference domain orientation bugs.  
- **TOOLS:** transfer, validation, assimilation, xml (first release).

## Version 1.6.0 (2024-05-29)
- **PROJECT:** soil moisture triple collocation  
  - **APPS – CELL:** app_map_create_tc  
    - Coastal smoothing, domain organization, time-selection tolerance fixes, weight method bug fixes.  
- **PROJECT:** validation framework  
  - **APPS – CELL:** app_cell_swi, app_cell_rzsm, app_cell_scaling, app_cell_metrics (first releases)  
  - **TOOLS – VALIDATION HSAF:** app_validation_main, app_validation_publisher (bulk options, bug fixes).

## Version 1.5.0 (2024-04-15)
- **PROJECT:** validation framework  
  - **APPS – CELL:** app_img2cell_gldas, app_img2cell_ecmwf (updates, new image_buffer option)  
  - **TOOLS – VALIDATION HSAF:** app_validation_main, app_validation_publisher (logging, option additions).

## Version 1.4.1 (2024-04-09)
- **PROJECT:** validation framework  
  - **APPS – CELL:** app_img2cell_gldas (georeference bug fix).

## Version 1.4.0 (2024-03-29)
- **PROJECT:** soil moisture rescaled (obs/mod)  
  - **APPS – MAP:** convert_cell2grid_ascat, convert_cell2grid_metrics (first releases)  
- **PROJECT:** validation framework  
  - **APPS – CELL:** app_img2cell_cci, app_img2cell_gldas (updates)  
  - **TOOLS – VALIDATION HSAF/SM:** app_validation_main (bug fixes).

## Version 1.3.0 (2024-02-28)
- **PROJECT:** soil moisture rescaled (obs/mod)  
  - **APPS – MAP:** convert_swath2cell (bug fixes, expanded product support, file management).  
- **PROJECT:** time-series  
  - **APPS – TS:** join_ts, sync_ts, analyze_ts, view_ts (first release).  
- **PROJECT:** utility framework  
  - **TOOLS:** transfer, validation, assimilation, xml (first release).  
- **PROJECT:** viewer framework  
  - **NOTEBOOK:** notebook_recolour_sm_ts.  
- **PROJECT:** validation framework  
  - **TOOLS – VALIDATION SM:** app_validation_main (first release).

## Version 1.2.0 (2023-12-19)
- **PROJECT:** soil moisture triple collocation  
  - **APPS:** create_grid_tc (added temporal periods, resampling, grid fixes, time selection fixes).

## Version 1.1.0 (2023-11-28)
- **PROJECT:** soil moisture triple collocation and time-series  
  - **APPS:** cell, maps, time-series; **TOOLS:** converter, downloader, plot_validation, plot_timeseries, validation, xml; **NOTEBOOKS:** time-series datasets and products.  
  - Refactored project structure, extended methods, fixed operational-mode bugs.

## Version 1.0.0 (2023-11-14)
- **PROJECT:** beta frameworks  
  - **APPS:** maps and time-series; **TOOLS:** validation, grid2ts, swath2ts, plotting, downloader, xml; **NOTEBOOKS:** time-series.

## Version 0.0.0 (2023-06-06)
- First commit and initialization of default settings.
