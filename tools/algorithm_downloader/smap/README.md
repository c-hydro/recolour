# README for downloader_codes

**2023/6/26**
started working on dowloaders: main script is in .py, launched by bash script (.sh) getting information from configuration file .json.
Added .json file for specific product @ 9 km SPL2SMP_E (see more notes in joplin).
- [x] review of bash script
- [x] review of configuration file: missing variables' names
- [x] review of main() called in script .py

**2023/6/30**
- [x] add --time_start / --time_end flags in command-line arguments/configuration file (instead of actual -time (now), days before and after)

**2023/7/13**
- [x] domain geotiff correct with -9999 no-data
- [x] 0 data in cut data tiff --> -9999