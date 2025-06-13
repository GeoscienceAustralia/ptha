#' Make list of breakwalls for swals to use.
#' 
#' Assumes they're all pre-processed into csv files within this directory.
#' Do this by running the make_breakwalls.R files in each subdirectory first.

all_breakwalls = list.files(pattern="\\.csv$", recursive=TRUE)

all_shaps = list.files(pattern="\\.shp$", recursive=TRUE)
# check that the number of csv files is the same as the number of shape files
if(length(all_breakwalls) != length(all_shaps)){
    stop('ERROR: Number of csv files does not match number of shape files.')
}

all_breakwalls = normalizePath(all_breakwalls)
cat(all_breakwalls, file='swals_breakwall_files.txt', sep="\n")
