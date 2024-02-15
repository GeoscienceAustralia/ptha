#files = read.csv('ordered_files_home_machine/Elevation_rasters_scraped_from_QGIS_file.csv')
files = read.csv('ordered_files_home_machine_rev1/Elevation_rasters_scraped_from_QGIS_file.csv')

writeLines(files$files, 'swals_elevation_files_in_preference_order.txt')
