# Run on NCI after running
#     Rscript create_elevation_preference_list_NCI.R
#
# Check if the md5sum of elevation files on NCI matches those on my home machine.
#
# If they don't match, it is possible that:
#   - One or other files is out of date
#   - They are both fine, but were created with different libraries (e.g. different gdal versions) that cause unimportant differences.
#

home_files = read.csv('ordered_files_home_machine_rev1/Elevation_rasters_scraped_from_QGIS_file.csv')

# NCI files
files = readLines('swals_elevation_files_in_preference_order.txt')
files_md5 = tools::md5sum(files)
names(files_md5) = ""
nci_files = data.frame(files=files, files_md5 = files_md5)

# Check for matches. 
match(nci_files$files_md5, home_files$files_md5)
##  [1]  1 NA  3  4 NA  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
## [76] 76 77 78

## UPDATED FILE LIST 2023/11/09
## > match(nci_files$files_md5, home_files$files_md5)
## [1]  1 NA  3  4 NA  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
## [76] 76 77 78 79 80 81 82 83 84

## UPDATED FILE LIST 2023/11/14
## > match(nci_files$files_md5, home_files$files_md5)
##  [1]  1 NA  3  4 NA  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
## [76] 76 77 78 79 80 81 82

## UPDATED 2023/11/15
## > match(nci_files$files_md5, home_files$files_md5)
##  [1]  1 NA  3  4 NA  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
## [76] 76 77 78 79 80 81 82

## UPDATED 2023/11/28
## > match(nci_files$files_md5, home_files$files_md5)
##  [1]  1 NA  3  4 NA  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
## [76] 76 77 78 79 80 81 82


#
# Above there are two 'non matches', one due to a vrt difference for the 2018
# bathytopo (here the contributing tifs have the same md5sum) and another due
# likely to a minor processing difference for a Gold Coast City DEM.

