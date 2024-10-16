#
# Append additional information to the boxes that SWALS can use
# 
library(rptha)

# Names of CSV files defining the nesting levels.
nesting_levels_to_edit = paste0(
    'domains_2_6_30/nesting_level_',
    c('1.csv', '2.csv', '3.csv'))

# For each element of nesting_levels_to_edit, name a shapefile with points
# where domains should be coarsened. The degree of coarsening will be defined
# in the SWALS multidomain_design_mod.f90. Empty string implies no coarsening
points_to_coarsen_shapefiles = c('', 'second_level_domains_to_coarsen/second_level_domains_to_coarsen.shp', '')
points_to_coarsen_layers = c("", "second_level_domains_to_coarsen", "")

stopifnot(length(nesting_levels_to_edit) == length(points_to_coarsen_shapefiles))

for(i in 1:length(nesting_levels_to_edit)){

    bte = read.csv(nesting_levels_to_edit[i])
    points_to_coarsen = points_to_coarsen_shapefiles[i]
    layer_name = points_to_coarsen_layers[i]

    if(points_to_coarsen != ""){
        # Flag that some boxes should be coarsened
        sites_to_keep_fine = coordinates(readOGR(points_to_coarsen, layer=layer_name))

        # Initialise to 'all fine'
        coarsen_box = rep(0, nrow(bte))

        # Maintain high resolution in boxes containing a sites_to_keep_fine point
        for(j in 1:nrow(sites_to_keep_fine)){
            has_a_point = (bte[,1] < sites_to_keep_fine[j,1]) & (bte[,3] > sites_to_keep_fine[j,1]) &
                          (bte[,2] < sites_to_keep_fine[j,2]) & (bte[,4] > sites_to_keep_fine[j,2])
            coarsen_box[has_a_point] = 1
        }
        
    }else{
        coarsen_box = rep(0, nrow(bte))
    }

    bte = cbind(bte, data.frame('coarsen' = coarsen_box))

    output_filename = gsub('.csv', '_edited.csv', nesting_levels_to_edit[i], fixed=TRUE)
    print(paste0("Saving to ", output_filename))
    write.csv(bte, output_filename, row.names=FALSE)
}
