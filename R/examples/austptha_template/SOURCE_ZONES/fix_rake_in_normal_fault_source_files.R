#
# DO NOT RUN THIS SCRIPT!
#
# Earlier codes had a bug whereby the unit_source_statistics.* files always had
# rake=90 written out, even for normal fault sources which should have rake=-90
#
# These values were also copied to the all_*_earthquake_events_tsunami_*.nc files
#
# The code has been corrected, but this code shows how to make a manual correction

run_me = FALSE

if(run_me){
    library(ncdf4)

    sourcezone_parameters = read.csv('../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv', stringsAsFactors=FALSE)

    sourcezones_to_correct = sourcezone_parameters$sourcename[which(sourcezone_parameters$rake == -90)]

    for(i in 1:length(sourcezones_to_correct)){
        sz = sourcezones_to_correct[i]
        print(sz)

        #
        # Fix the csv
        #
        unit_source_statistics_csv = Sys.glob(paste0(sz, '/TSUNAMI_EVENTS/unit_source_statistics*.csv'))
        if(length(unit_source_statistics_csv) == 1){
            uss = read.csv(unit_source_statistics_csv)

            # Fix the rake
            uss_new = uss
            uss_new$rake = -90.0

            write.csv(uss_new, file=unit_source_statistics_csv, row.names=FALSE)

            print(paste0('..updated csv ', unit_source_statistics_csv))
        }else{
            print('..no csv')
        }

        #
        # Fix the nc file
        #
        unit_source_statistics_nc = Sys.glob(paste0(sz, '/TSUNAMI_EVENTS/unit_source_statistics*.nc'))
        if(length(unit_source_statistics_nc) == 1){
            fid = nc_open(unit_source_statistics_nc, write=TRUE)
            rake = ncvar_get(fid, 'rake')

            new_rake = rep(-90L, length(rake))

            ncvar_put(fid, 'rake', new_rake)

            nc_close(fid)
            print(paste0('..updated nc ', unit_source_statistics_nc))
        }else{
            print('..no nc')
        }

        #
        # Fix the other nc files
        #
        tsunami_events_nc = Sys.glob(paste0(sz, '/TSUNAMI_EVENTS/all_*_earthquake_events_tsunami*.nc'))
        if(length(tsunami_events_nc)>0){

            for(j in 1:length(tsunami_events_nc)){
                fid = nc_open(tsunami_events_nc[j], write=TRUE, readunlim=FALSE)
                rake = ncvar_get(fid, 'us_rake')
                new_rake = rep(-90.0, length(rake))
                ncvar_put(fid, 'us_rake', new_rake)
                nc_close(fid)
                print(paste0('..updated nc ', tsunami_events_nc[j]))
            }

        }else{
            print('..no tsunami events nc')
        }

    }

}else{
    print('By default this script does nothing -- but it illustrates a work-around for a bug (which is now fixed), see comments above')
    print('Thus you should not need to use this script')
}

