#' Make list of breakwalls for swals to use.
#' 
#' Assumes they're all pre-processed into csv files within this directory.
#' Do this by running the make_breakwalls.R files in each subdirectory first.


all_breakwalls = normalizePath(list.files(pattern="\\.csv$", recursive=TRUE))
cat(all_breakwalls, file='swals_breakwall_files.txt', sep="\n")
