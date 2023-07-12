library(sp)

# Get a world-map, but protect against deprecation of maptools
if(!file.exists('wrld_simpl.RDS')){
    # Note the 'maptools' package is being retired.
    library(maptools)
    data(wrld_simpl)
    saveRDS(wrld_simpl, 'wrld_simpl.RDS')
}else{
    wrld_simpl = readRDS('wrld_simpl.RDS')
}
