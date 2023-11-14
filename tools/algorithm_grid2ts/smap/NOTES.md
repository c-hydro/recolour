# NOTES on algorithm_grid2ts general structure

Legend
$$$ = to rewrite

## hmc section
Start: app_grid2ts_hmc

main_wrapper(): wrapper around reshuffle.main()
- gets arguments from command line and from configuration file (from configuration file gets args to pass into reshuffle.main(args))
- sets log
- gets, parses data settings from configuration file .json
- create_file_grid: defined in lib_utils_hmc, reads a pre-existent tif or tiff file and uses it to produce a grid which is compatible with the product; can read any pre-existent map which has the same gridding of your product
- main_runner: calls reshuffle.main()

reshuffle.main(args):
- args must be compatible with the call in reshuffle()
- **hmc.grid.load_grid(land_points=True, bbox=None, grid_path=None):**
  - can pass land_points (default: True) and hmc_land_grid (alias for hms_grids) will retrieve the grid stored in grid_path; hmc_grids (stays in hmc.grid) creates a global 0.25 deg gldas grid, eventually subset to bbox
  - the grid is, as always, a BasicGrid with .to_cell_grid method to provide correct cells of 5. degs
  - since land_points is True by default, it is recommended to pass a land grid as a nc file 
  - optional bbox (min_lat, min_lon, max_lat, max_lon)
  - returns a BasicGrid 
- **reshuffle(args.dataset_grid_root, args.dataset_ts_root, args.start, args.end, args.grid_path, args.parameters, file_name_tmpl, datetime_format, input_grid=input_grid, img_buffer=args.imgbuffer)**
  - for all args.* meaning, see parse_args in the same file:
    - dataset_grid_root: folder where the source data (tifs) are stored
    - start, end = start date, end date
    - grid_path: path to reference grid (nc)
    - parameters: variables to reshuffle into ts
    - --imgbuffer: optional, # images to reshuffle at once
    - --land_points: optional
    - --bbox: optional
  - needs input_grid as reference grid (made with load_grid)
  - $$$ dataset_grid_driver = interface.hmd_ds: file name template is hardcoded here (must think of changing it and providing in configuration file)

$$$
hmc.interface():
classes and methods to read 2D tif images and 1D nc ts before and after reshuffling

**_class hms_ds(MultiTemporalImageBase):_**
"""this is a child class of MultiTemporalImageBase

The constructor of the parent class is called (super().__init__()) and the stack of images is passed using a wrapper HMC_Wrap_Img around the single-image HMC_Base_Img, just to pass arguments.
"""

**_class HMC_Base_Img(ImageBase):_**
"""classic re-writing of custom TUW ImageBase"""
**__init__(self,
filename, mode='r', parameter='soil_moisture', grid_path=None, subgrid=None, array_1D=False)**
- self.fill_values fixed to -9999
- self.grid = hmc_cell_grid(grid_path) (alias to hmc_grids, returns a BasicGrid with cells) or subgrid

**read(self, timestamp=None)**
- WARNING: reshape of image is fixed to (720, 1440) cells, may need to be changed
$$$

(don't feel lost, go back to reshuffle)

drv_reshuffle = Img2Ts(
        input_dataset=dataset_grid_driver,
        outputpath=dataset_ts_root,
        startdate=start_date,
        enddate=end_date,
        input_grid=grid,
        imgbuffer=img_buffer,
        cellsize_lat=5.0,
        cellsize_lon=5.0,
        global_attr=global_attr,
        zlib=True,
        unlim_chunksize=None,
        ts_attributes=ts_attributes,
)

lib_img2ts_hmc
**_class Img2Ts()_**
"""classic TUW img2ts with nc magic and a mess"""

---

## things to check
- [ ] shape of tif files
- [x] do i have all bands? i remeber something about latitude and longitude... YES THE GRID NC FILE IS PRODUCED WITH ALL BANDS WHICH ARE NEEDED
- [x] can i produce a nc grid as needed? YES THE FUNCTION create_file_grid DOES IT FOR ME

---

## smap section
just the list of matrioska methods, un-matriosked
note: if passages above are not present, it may be that the code can run them correctly

app_grid2ts_smap
    main_runner() = smap.reshuffle.main()
-> smap.reshuffle.main()
    input_grid = smap.grid.load_grid()
    dataset_grid_driver = smap.interface.smap_ds()
    drv_reshuffle = Img2Ts()
    drv_reshuffle.calc()


### managing timestamps
reshuffle
    data = dataset_grid_driver.read(start_date)
-> pygeobase.io_base.read(timestamp, **kwargs)
    pygeobase.io_base._assemble_img(timestamp, **kwargs)
    