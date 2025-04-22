=========
Changelog
=========

Version 2.0.0 [2025-04-22]
**************************
- PROJECT: operational framework grid	
	- APPS -- merge_grid2ref
		- first release to megre partial grid to the reference one
	- APPS -- convert_time_step_src2csv
		- extend settings to the operational mode run
	- TOOLS -- DOWNLOADER SMAP: smap_downloader_spl2smp_e.py
		- refactor of all application methods and algorithm

Version 1.9.0 [2025-02-27]
**************************

- PROJECT: operational framework soil moisture triple collocation
	- APPS -- CELL: app_map_create_tc
		- fix bugs
	- TOOLS -- DOWNLOADER SMAP: smap_downloader_spl2smp_e.py
		- fix bugs in managing time information
	- TOOLS -- DOWNLOADER SMAP: smap_downloader_spl3smp_e.py
		- fix bugs in managing time information

Version 1.8.0 [2024-09-30]
**************************
- PROJECT: operational framework time-series
	- APPS -- convert_grid_hmc2csv
		- fix bugs related to the dumping of time-series (delimiter, no_data and formats for matching to the other ts methods)
	- APPS -- convert_time_step_src2csv
		- fix bugs related to time_start/time_end in realtime applications
	- APPS -- sync_ts
		- fix bugs related to unavailable point or datasets
		- extend method to resample time-series with different frequency
	- APPS -- analyze_ts
		- fix bugs related to unavailable point or datasets
		- add options to rescale methods to make the analysis only in the variable range
	- APPS -- view_ts
		- fix bugs related to unavailable point or datasets
		- extend heatmaps features using the seaborn package

Version 1.7.0 [2024-07-09]
**************************

- PROJECT: operational framework soil moisture rescaled (obs/mod)
	- APPS -- MAP: create_map_smr
		- first release based on old soil moisture packages

- PROJECT: operational framework grid	
	- APPS -- TS: resample_grid_src2ref
		- first release

- PROJECT: operational framework time-series	
	- APPS -- TS: convert_grid_ecmwf2csv, convert_grid_hmc2csv, convert_grid_ascat2csv, convert_grid_smap2csv, convert_grid_gldas2csv
		- fix bugs related to the orientation of the reference domain (in tiff format)
	- APPS -- TS: convert_grid_hmc2csv
		- add format to read all hmc datasets (static and dynamics)

Version 1.6.0 [2024-05-29]
**************************
- PROJECT: operational framework soil moisture triple collocation
	- APPS -- CELL: app_map_create_tc
		- update codes to smooth the coasts areas to avoid artefacts
		- update codes to organize domain datasets and metrics using the same resampling and filtering strategy (testing)
		- update codes to correct the time selection based on tolerance period
		. fix bugs in weight method (case: ref found, k1 not found and k2 found)

- PROJECT: operational framework soil moisture rescaled (obs/mod)
	- APPS -- CELL: app_cell_swi
		- first release to convert the ssm to swi time-series
	- APPS -- CELL: app_cell_rzsm
		- first release to convert the rszm layers to rzsm profile time-series
	- APPS -- CELL: app_cell_scaling
		- first release to scale the nrt time-series using a reference time-series dataset
	- APPS -- CELL: app_cell_metrics
		- first release to compute time-series metrics

- PROJECT: validation framework 
	- APPS -- CELL: app_img2cell_gldas
		- fix bugs in geographical orientation
	- APPS -- CELL: app_img2cell_ecmwf
		- fix bugs in geographical orientation
	- APPS -- CELL: app_img2cell_cci
		- fix bugs in geographical orientation
	- TOOLS -- VALIDATION HSAF: app_validation_main
		- add bulk option in the reference dataset
	- TOOLS -- VALIDATION HSAF: app_validation_publisher
		- fix bugs related to the old datasets

Version 1.5.0 [2024-04-15]
**************************
- PROJECT: validation framework 
	- APPS -- CELL: app_img2cell_gldas
		- update codes
	- APPS -- CELL: app_img2cell_ecmwf
		- update codes and add image_buffer option in the settings file (to manage nrt and dr applications)
	- TOOLS -- VALIDATION HSAF: app_validation_main
		- add options and logging features avaialable in the previous versions outside the recolour package
	- TOOLS -- VALIDATION HSAF: app_validation_publisher
		- update codes based on the previous versions (2017-2022) and adapt the scripts to different configurations

Version 1.4.1 [2024-04-09]
**************************
- PROJECT: validation framework 
	- APPS -- CELL: app_img2cell_gldas
		- fix bug related to the georeference information in the reshuffle tool

Version 1.4.0 [2024-03-29]
**************************
- PROJECT: operational framework soil moisture rescaled (obs/mod)
	- APPS -- MAP: convert_cell2grid_ascat
		- first release (product h16 and h103)
	- APPS -- MAP: convert_cell2grid_metrics
		- first release (product ascat and ecmwf)

- PROJECT: validation framework 
	- APPS -- CELL: app_img2cell_cci
		- update codes
	- APPS -- CELL: app_img2cell_gldas
		- update codes
	- TOOLS -- VALIDATION HSAF: app_validation_main
		- update codes and fix bugs 
	- TOOLS -- VALIDATION SM: app_validation_main
		- update codes and fix bugs 
	
Version 1.3.0 [2024-02-28]
**************************
- PROJECT: operational framework soil moisture rescaled (obs/mod)
	- APPS -- MAP: convert_swath2cell
		- fix bugs
		- update code to product h16, h103, h104 and h105
		- update code to manage tmp file (to check the long analysis)

- PROJECT: operational framework time-series	
	- APPS -- TS: join_ts, sync_ts, analyze_ts, view_ts
		- first release and fix bugs

- PROJECT: utility framework
	- TOOLS: transfer, validation, assimilation and xml
		- first release and fix bugs
	
- PROJECT: viewer framework
	- NOTEBOOK: notebook_recolour_sm_ts

- PROJECT: validation framework	
	- TOOLS -- VALIDATION SM: app_validation_main
		- first release and fix bugs

Version 1.2.0 [2023-12-19]
**************************
- PROJECT: operational framework soil moisture triple collocation
	- APPS: create_grid_tc
		- add temporal periods to match available products (reference, k1 and k2)
		- add resampling procedure to remap products k1 and k2 to the reference grid
		- fix artetacts in k1 and k2 products (due to the generic grid reference)
		- fix selection of time for reference, k1 and k2 products

Version 1.1.0 [2023-11-28]
**************************
- PROJECT: operational framework soil moisture triple collocation and time-series
	- APPS: cell, maps and time-series
	- TOOLS: converter, downloader, plot_validation, plot_timeseries, validation, xml
	- NOTEBOOKS: time-series datasets and products

- Refactor project structure and codes
- Extend methods and functions of img2cell, swath2cell, ecmwf2ts, hmc2ts and smap2ts
- Fix bugs (for operational mode)

Version 1.0.0 [2023-11-14]
**************************
- PROJECT: beta frameworks
	- APPS: maps and time-series
	- TOOLS: validation, grid2ts, swath2ts, plotting, downloader, xml
	- NOTEBOOKS: time-series

Version 0.0.0 [2023-06-06]
**************************
- PROJECT: first commit to open the repository and initialize the default settings
