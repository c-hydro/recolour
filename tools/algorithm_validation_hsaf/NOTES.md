# NOTES for tools/algorithm_validation



Note:
- smap timeseries and data are produced almost exactly matching hmc ts and data, so the hmc libraries should be your reference throught the process
- notes are a list of functions that make relevant things in the process, in the order in which they appear
- if there's no checkbox, there's nothing (apparently) to work on/modify/update


## Tasks
1. [x] use hmc as reference
2. [x] rewrite code to calculate metrics over hmc and ascat only; check lib_utils_metrics


---

## Task 1

## Review of existing code **validation**

-> app_validation_main.py
    - drv_process.DrvProcess()
    - drv_process.setup_process()
    - drv_process.execute_process()

-> drv_process_main.py
- imports: import cpl_process_* (all the libraries)
**def setup_process()**
    - self.coupler_time = CplTime()
    - self.coupler_time.setup_time()
    -> cpl_process_time.py

    - self.coupler_datasets = CplDatasets()
    - dset_interfaces, dset_modes = self.coupler_datasets.setup_datasets_src()
    - self.dset_path_dst = self.coupler_datasets.setup_datasets_dst()
    -> cpl_process_datasets.py
        - [x] add dataset type, variable, name SMAP
        - [x] **def set_dset_reader()**: has hardcoded information, update

    - self.coupler_metrics = CplMetrics()

    - self.coupler_analysis = CplAnalysis()

    - self.coupler_mode = CplMode()

- [x] -> lib_datasets_smap.py
- [x] -> lib_interface_smap.py

---

## Rewriting of existing code

-> lib_interface_smap.py
written based on algorithm_grid2ts/smap/interface and lib_interface_hmc. The first script contains a lot more functions and classes than are needed, so only the last one (smap_ts) is kept as it is the interface to GriddedNcOrthoMultiTs. 

-> lib_datasets_smap.py
written directly based on lib_datasets_hmc: only changed names of functions (almost)

**class SMAP_Dataset(SMAPTs)**
main wrapper around the timeseries reader.
if 'sm_nan' in kwargs:
    self.sm_nan = kwargs.pop('sm_nan')
else:
    self.sm_nan = np.nan
        # default changed from -9999. to nan
        # should check behaviour with both values passed in kwargs

---

## Debugging

**issue #1**

2023-07-26 12:54:01,601 root         INFO       (2) Process Cell(s): "1394"                                                    drv_process_main.py:[216    - exec_process_wrap   ()] 
/home/hsaf/recolour/conda/envs/recolour_libraries/lib/python3.7/site-packages/pygeogrids/nearest_neighbor.py:219: UserWarning: Less than k=1 points found within max_dist=35000. Distance set to 'Inf'.
  warnings.warn(f"Less than k={k} points found within "
/home/hsaf/recolour/conda/envs/recolour_libraries/lib/python3.7/site-packages/pytesmo/validation_framework/validation.py:141: UserWarning: You are using the default temporal matcher. If you are using one of the newer metric calculators (PairwiseIntercomparisonMetrics, TripleCollocationMetrics) you should probably use `make_combined_temporal_matcher` instead. Have a look at the documentation of the metric calculators for more info.
  "You are using the default temporal matcher. If you are using "

- seems like lut_max_dist = 35000 is too small, automatically set to inf
    - changed to 1000000, still no points found
--> add CHECK: print lon, lats inside pygeogrids.nearest_neighbor.calc_lut()
    - unexpected values (from 9 to -99 lon (only see array's extremes) and 0 to 80 lat) so check if wth gldas is the same --> YES

CHECK:  [9.432300567626953 9.54458999633789 9.656878471374512 ...
 -99.01937866210938 99.65821838378906 -99.65821838378906]
CHECK:  [0.0 0.0 0.0 ... 79.90965270996094 79.90965270996094 79.90965270996094]


this is the path to follow in debugging, breakpoint after breakpoint (read from bottom (shallowest) to top (deepest))

Breakpoint reached: nearest_neighbor.py:217
Stack: 
	find_nearest_index, nearest_neighbor.py:217
	calc_lut, grids.py:580
	get_luts, data_manager.py:183
	__init__, data_manager.py:148
	__init__, validation.py:136
	setup_analysis, cpl_process_analysis.py:53
	exec_process_wrap, drv_process_main.py:228
	execute_process, drv_process_main.py:172
	main, app_validation_main.py:110
	<module>, app_validation_main.py:203 


**issue #2**
at some point HMC stops being reference and ascat becomes, dunno why

try with:

calc, validation.py:190 - rename_cols=False
--> cannot do since like this you don't have 'ref' nor 'k1', 'k2' and calculator doesn't recognize columns' names



---


## Task 2

rewrite code to calculate metrics over hmc and ascat only; check lib_utils_metrics.
now using:
lib_utils_metrics.py:
    **class ExtendedMetrics**

Have to implement usage of:
lib_utils_metrics.py:
    **class BasicMetrics**


---


## Task 3

implement reading of hsaf product h14 (ecmwf)







