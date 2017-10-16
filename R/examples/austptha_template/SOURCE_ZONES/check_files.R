#
# This counts the number of unit-sources, netcdf files, and polygons-in-the-unit-source-grid
# If everything has run correctly, the numbers should be identical
#

library(rgdal, quietly=TRUE)

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
          ifelse( class(x$mypol) != 'try-error', length(unique(x$mypol$dwndp_n)), -1),
          x$match,
          x$has_tsunami_max_stage_pdf,
          x$has_unit_source_statistics,
          x$has_check_event_png
         ) }
    )


##  $alaskaaleutians
##  [1] "312"      "312"      "312"      "4"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $arutrough
##  [1] "6"        "6"        "6"        "1"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $banda_detachment
##  [1] "26"       "26"       "26"       "2"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $cascadia
##  [1] "66"       "66"       "66"       "3"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $flores
##  [1] "80"       "80"       "80"       "2"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $hjort
##  [1] "26"       "26"       "26"       "2"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $izumariana
##  [1] "140"      "140"      "140"      "2"        "Match OK" "TRUE"     "TRUE"    
##  [8] "TRUE"    
##  
##  $kermadectonga
##  [1] "207"      "207"      "207"      "3"        "Match OK" "TRUE"     "TRUE"    
##  [8] "TRUE"    
##  
##  $kurilsjapan
##  [1] "244"      "244"      "244"      "4"        "Match OK" "TRUE"     "TRUE"    
##  [8] "TRUE"    
##  
##  $macquarienorth
##  [1] "13"       "13"       "13"       "1"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $makran
##  [1] "96"       "96"       "96"       "6"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $manokwari
##  [1] "7"             "14"            "7"             "1"            
##  [5] "Match failure" "FALSE"         "FALSE"         "FALSE"        
##  
##  $manus
##  [1] "36"            "2"             "36"            "1"            
##  [5] "Match failure" "FALSE"         "FALSE"         "FALSE"        
##  
##  $mexico
##  [1] "174"      "174"      "174"      "3"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $moresby_trough
##  [1] "26"                       "0"                       
##  [3] "26"                       "2"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $mussau
##  [1] "8"                        "0"                       
##  [3] "8"                        "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $newguinea
##  [1] "46"                       "0"                       
##  [3] "46"                       "2"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $newhebrides
##  [1] "68"       "68"       "68"       "2"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $north_sulawesi
##  [1] "44"       "44"       "44"       "4"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $outer_rise_timor
##  [1] "22"       "22"       "22"       "1"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $outerrise_kermadectonga
##  [1] "54"                       "0"                       
##  [3] "54"                       "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $outerrise_puysegur
##  [1] "14"                       "0"                       
##  [3] "14"                       "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $outerrisenewhebrides
##  [1] "35"                       "0"                       
##  [3] "35"                       "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $outerrisesolomon
##  [1] "35"                       "0"                       
##  [3] "35"                       "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $outerrisesunda
##  [1] "122"      "122"      "122"      "1"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $philippine
##  [1] "66"                       "0"                       
##  [3] "66"                       "2"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $puysegur
##  [1] "34"       "34"       "34"       "2"        "Match OK" "TRUE"     "TRUE"    
##  [8] "TRUE"    
##  
##  $ryuku
##  [1] "120"      "120"      "120"      "3"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $sandwich
##  [1] "28"       "28"       "28"       "2"        "Match OK" "TRUE"     "FALSE"   
##  [8] "FALSE"   
##  
##  $sangihe
##  [1] "26"       "26"       "26"       "2"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $sangihe_backthrust
##  [1] "30"       "30"       "30"       "2"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $se_sulawesi
##  [1] "11"       "11"       "11"       "1"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $seram_thrust
##  [1] "24"       "24"       "24"       "2"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $seramsouth
##  [1] "5"                        "0"                       
##  [3] "5"                        "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   
##  
##  $solomon
##  [1] "92"       "92"       "92"       "2"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $southamerica
##  [1] "660"      "660"      "660"      "4"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $sunda
##  [1] "492"      "492"      "492"      "4"        "Match OK" "TRUE"     "TRUE"    
##  [8] "FALSE"   
##  
##  $tanimbar
##  [1] "12"       "12"       "12"       "2"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $timor
##  [1] "66"       "66"       "66"       "3"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $tolo_thrust
##  [1] "6"        "6"        "6"        "1"        "Match OK" "FALSE"    "FALSE"   
##  [8] "FALSE"   
##  
##  $trobriand
##  [1] "12"                       "0"                       
##  [3] "12"                       "1"                       
##  [5] "... cannot attempt match" "FALSE"                   
##  [7] "FALSE"                    "FALSE"                   

