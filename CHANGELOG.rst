=========
Changelog
=========


Version 1.6.0 [2024-05-02]
**************************
- PROJECT: validation framework 
	- APPS -- CELL: app_cell_swi
		- first release to convert the ssm to swi time-series
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
- PROJECT: operational framework
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
- PROJECT: operational framework
	- APPS -- MAP: convert_swath2cell
		- fix bugs
		- update code to product h16, h103, h104 and h105
		- update code to manage tmp file (to check the long analysis)
	
	- APPS -- TS: join_ts, sync_ts, analyze_ts, view_ts
		- first relaase and fix bugs
	
	- TOOLS: transfer, validation, assimilation and xml
		- first release and fix bugs
	
- PROJECT: viewer framework
	- NOTEBOOK: notebook_recolour_sm_ts

- PROJECT: validation framework	
	- TOOLS -- VALIDATION SM: app_validation_main
		- first release and fix bugs

Version 1.2.0 [2023-12-19]
**************************
- PROJECT: operational framework
	- APPS: create_grid_tc
		- add temporal periods to match available products (reference, k1 and k2)
		- add resampling procedure to remap products k1 and k2 to the reference grid
		- fix artetacts in k1 and k2 products (due to the generic grid reference)
		- fix selection of time for reference, k1 and k2 products

Version 1.1.0 [2023-11-28]
**************************
- PROJECT: operational framework
	- APPS: cell, maps and time-series
	- TOOLS: converter, downloader, plot_validation, plot_timeseries, validation, xml
	- NOTEBOOKS: time-series datasets and products

- Refactor project structure and codes
- Extend methods and functions of img2cell, swath2cell, ecmwf2ts, hmc2ts and smap2ts
- Fix bugs (for operational mode)

Version 1.0.0 [2023-11-14]
**************************
- PROJECT: beta framework
	- APPS: maps and time-series
	- TOOLS: validation, grid2ts, swath2ts, plotting, downloader, xml
	- NOTEBOOKS: time-series

Version 0.0.0 [2023-06-06]
**************************
- PROJECT: first commit to open the repository and initialize the default settings
