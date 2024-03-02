#' Make list of breakwalls for swals to use.
#' 
#' Assumes they're all pre-processed into csv files within this directory.
#' Do this by running the make_breakwalls.R files in each subdirectory first.


all_breakwalls = list.files(pattern="\\.csv$", recursive=TRUE)
# DELIBERATELY IGNORE BUNBURY FLOOD-GATE - the code will manually add it later.
all_breakwalls = all_breakwalls[!grepl("bunbury_floodgate", all_breakwalls)]
all_breakwalls = normalizePath(all_breakwalls)
cat(all_breakwalls, file='swals_breakwall_files.txt', sep="\n")
