#
# This counts the number of unit-sources, netcdf files, and polygons-in-the-unit-source-grid
# If everything has run correctly, the numbers should be identical
#

library(rptha, quietly=TRUE)

all_dir = dirname(dirname(Sys.glob('*/EQ_SOURCE/config.R')))
all_dir = setdiff(all_dir, 'TEMPLATE') # remove TEMPLATE

# Base output directories
output_base = paste0('/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/', all_dir)

# Loop over all sourcezones, and count files, etc
file_list = list()
for(i in 1:length(output_base)){
    file_list[[i]] = list()
    
    nme = all_dir[i]
    ob = output_base[i]

    print(nme)

    file_list[[i]]$nme = nme
    file_list[[i]]$ob = ob

    unit_source_tifs = Sys.glob(paste0(ob, '/EQ_SOURCE/Unit_source_data/*/*.tif'))
    file_list[[i]]$num_unit_source_tifs = length(unit_source_tifs)
    file_list[[i]]$unit_source_tifs = unit_source_tifs

    netcdf_files = Sys.glob(paste0(ob, '/TSUNAMI_UNIT_SOURCES/*/*/*/*.nc'))
    file_list[[i]]$netcdf_files = netcdf_files
    file_list[[i]]$num_netcdf_files = length(netcdf_files)

    if(length(netcdf_files) > 0){
        # Check that 'station' is the unlimited dimension in the netcdf file
        fid = nc_open(netcdf_files[1], readunlim=FALSE)
        file_list[[i]]$station_is_unlim = fid$dim$station$unlim
        nc_close(fid)
    }else{
        file_list[[i]]$station_is_unlim=FALSE
    }

    # Read the unit source grid shapefile
    mypol = try(readOGR(paste0(ob, '/EQ_SOURCE/unit_source_grid/', nme, '.shp'),
        layer=nme, verbose=FALSE))
    file_list[[i]]$mypol = mypol
    file_list[[i]]$num_unit_sources = try(length(mypol))

    if(length(netcdf_files) > 0){
        uss_name = gsub('.tif', '/', basename(unit_source_tifs))
        match_count = sapply(uss_name, f<-function(x) sum(grepl(x, netcdf_files, fixed=TRUE)))
        if(any(match_count != 1)){
            file_list[[i]]$match = 'Match failure'
        }else{
            file_list[[i]]$match = 'Match OK'
        }
    }else{
        file_list[[i]]$match = '... cannot attempt match'
    }

    file_list[[i]]$has_tsunami_max_stage_pdf = (length(Sys.glob(paste0(nme, '/TSUNAMI_UNIT_SOURCE/*.pdf'))) == 1)
    file_list[[i]]$has_unit_source_statistics = (length(Sys.glob(paste0(nme, '/TSUNAMI_EVENTS/unit_source_statistics*.nc'))) == 1)
    file_list[[i]]$has_check_event_png = (length(Sys.glob(paste0(nme, '/TSUNAMI_EVENTS/event_size_scaling_stochastic_*.png'))) == 1)
      
}
names(file_list) = all_dir

# Check number of unit source tifs, number of unit source regions, number of
# netcdf files, number of rows downdip, and whether we have the max-stage pdf,
# the unit_source_statistics netcdf file, and a png file evidencing check of scaling relations
lapply(file_list, 
    f<-function(x) {
        c(x$num_unit_source_tifs, 
          x$num_netcdf_files, 
          x$num_unit_sources, 
          x$station_is_unlim,
          ifelse( class(x$mypol) != 'try-error', length(unique(x$mypol$dwndp_n)), -1),
          x$match,
          x$has_tsunami_max_stage_pdf,
          x$has_unit_source_statistics,
          x$has_check_event_png
         ) }
    )

## EXAMPLE OUTPUT WHEN EVERYTHING IS DONE ##
##
## 
## $alaskaaleutians
## [1] "312"      "312"      "312"      "TRUE"     "4"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $arutrough
## [1] "6"        "6"        "6"        "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $banda_detachment
## [1] "26"       "26"       "26"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $cascadia
## [1] "66"       "66"       "66"       "TRUE"     "3"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $flores
## [1] "80"       "80"       "80"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $hjort
## [1] "26"       "26"       "26"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $izumariana
## [1] "140"      "140"      "140"      "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $kermadectonga
## [1] "207"      "207"      "207"      "TRUE"     "3"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $kurilsjapan
## [1] "244"      "244"      "244"      "TRUE"     "4"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $macquarienorth
## [1] "13"       "13"       "13"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $makran
## [1] "96"       "96"       "96"       "TRUE"     "6"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $manokwari
## [1] "7"        "7"        "7"        "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $manus
## [1] "36"       "36"       "36"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $mexico
## [1] "174"      "174"      "174"      "TRUE"     "3"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $moresby_trough
## [1] "26"       "26"       "26"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $mussau
## [1] "8"        "8"        "8"        "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $newguinea
## [1] "46"       "46"       "46"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $newhebrides
## [1] "68"       "68"       "68"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $north_sulawesi
## [1] "44"       "44"       "44"       "TRUE"     "4"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outer_rise_timor
## [1] "22"       "22"       "22"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outerrise_kermadectonga
## [1] "54"       "54"       "54"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outerrise_puysegur
## [1] "14"       "14"       "14"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outerrisenewhebrides
## [1] "35"       "35"       "35"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outerrisesolomon
## [1] "35"       "35"       "35"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $outerrisesunda
## [1] "122"      "122"      "122"      "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $philippine
## [1] "66"       "66"       "66"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $puysegur
## [1] "34"       "34"       "34"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $ryuku
## [1] "120"      "120"      "120"      "TRUE"     "3"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $sandwich
## [1] "28"       "28"       "28"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $sangihe
## [1] "26"       "26"       "26"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $sangihe_backthrust
## [1] "30"       "30"       "30"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $se_sulawesi
## [1] "11"       "11"       "11"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $seram_thrust
## [1] "24"       "24"       "24"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $seramsouth
## [1] "5"        "5"        "5"        "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $solomon
## [1] "92"       "92"       "92"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $southamerica
## [1] "660"      "660"      "660"      "TRUE"     "4"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $sunda
## [1] "492"      "492"      "492"      "TRUE"     "4"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $tanimbar
## [1] "12"       "12"       "12"       "TRUE"     "2"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $timor
## [1] "66"       "66"       "66"       "TRUE"     "3"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $tolo_thrust
## [1] "6"        "6"        "6"        "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
## $trobriand
## [1] "12"       "12"       "12"       "TRUE"     "1"        "Match OK" "TRUE"    
## [8] "TRUE"     "TRUE"    
## 
