To prepare elevation data for the model, initially I made a QGIS session that included many elevation files, ordered such that in combination they produced a "good" elevation field. 

I then saved this session in the 'old' QGIS format. That produces 2 files, one of which is plain text and easy to parse (initially `FILE_TO_SCRAPE_ELEVATION_RASTERS_FROM.txt` -- later updated to `Revised_file_to_scrape_elevation_files_from.qgs`).

The script `extract_elevation_filenames.R` was then used to create a file with the elevation rasters of interest in their preference order (`Elevation_rasters_scraped_from_QGIS_file.csv`). NOTE: THIS WAS RUN IN ANOTHER FOLDER ON MY HOME MACHINE. The path normalisation step won't work if run from this folder.


