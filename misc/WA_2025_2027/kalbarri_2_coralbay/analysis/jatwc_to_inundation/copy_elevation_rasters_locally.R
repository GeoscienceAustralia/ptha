#
# Make a local copy of the elevation rasters
#

ptha_elevation_rasts_list = list(
    # This file does NOT have the tidal adjustment
    no_tidal_adjustment = Sys.glob('../../swals/OUTPUTS/smallTestCase_kalbarri2coralbay_B_notidaladjustment-test_load_balance-ambient_sea_level_0/RUN_20260529_113644043/elevation0*.tif'),
    # This files DOES have the tidal adjustment
    with_tidal_adjustment = Sys.glob('../../swals/OUTPUTS/smallTestCase_kalbarri2coralbay_B_tidaladjustment-test_load_balance-ambient_sea_level_0/RUN_20260529_113659708/elevation0*.tif')
    )


for(nm in names(ptha_elevation_rasts_list)){

    # Make a local elevation folder and copy the rasters to it
    elev_dir = paste0('./elevation_in_model_', nm)
    dir.create(elev_dir)
    copy_worked = file.copy(
        from=ptha_elevation_rasts_list[[nm]], 
        to=paste0(elev_dir, '/', basename(ptha_elevation_rasts_list[[nm]])))
    if(!all(copy_worked)) stop('Copy of elevation data failed')

    # Make a vrt as well
    mydir = getwd()
    setwd(elev_dir)
    system('gdalbuildvrt -resolution highest all_elevation_combined.vrt elevation0*.tif')
    setwd(mydir)
}
