=========
Changelog
=========

Version 1.2.0 [2023-12-19]
**************************
- PROJECT: operational release
	- APPS: create_grid_tc
		- add temporal periods to match available products (reference, k1 and k2)
		- add resampling procedure to remap products k1 and k2 to the reference grid
		- fix artetacts in k1 and k2 products (due to the generic grid reference)
		- fix selection of time for reference, k1 and k2 products

Version 1.1.0 [2023-11-28]
**************************
- PROJECT: operational release
	- APPS: cell, maps and time-series
	- TOOLS: converter, downloader, plot_validation, plot_timeseries, validation, xml
	- NOTEBOOKS: time-series datasets and products

- Refactor project structure and codes
- Extend methods and functions of img2cell, swath2cell, ecmwf2ts, hmc2ts and smap2ts
- Fix bugs (for operational mode)

Version 1.0.0 [2023-11-14]
**************************
- PROJECT: beta release
	- APPS: maps and time-series
	- TOOLS: validation, grid2ts, swath2ts, plotting, downloader, xml
	- NOTEBOOKS: time-series

Version 0.0.0 [2023-06-06]
**************************
- PROJECT: first commit to open the repository and initialize the default settings
