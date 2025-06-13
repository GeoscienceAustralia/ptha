# check that each shapefile has a corresponding csv file
test_made_all_csv <- function() {
    # get all the shapefiles
    all_shp <- Sys.glob('*/*/*.shp')
    # get all the csv files
    all_csv <- Sys.glob('*/*/*.csv')
    # check that each shapefile has a corresponding csv file
    stopifnot(all_shp %in% gsub('csv', 'shp', all_csv))
}
test_made_all_csv()

# check that each csv file has no NA values
test_no_na_in_csv <- function() {
    all_csv <- Sys.glob('*/*/*.csv')
    for (csv in all_csv) {
        df <- read.csv(csv)
        stopifnot(!any(is.na(df)))
    }
}
test_no_na_in_csv()
