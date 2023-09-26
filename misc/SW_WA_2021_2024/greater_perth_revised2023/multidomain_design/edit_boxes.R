#
# Append additional information to the boxes that SWALS can use
# 
library(rptha)

# Names of CSV files defining the nesting levels.
nesting_levels_to_edit = paste0('domains_2023_08_29_0.5_0.166666666666667_0.0333333333333333/',
    c('first_level_nesting.csv', 'second_level_nesting.csv', 'third_level_nesting.csv'))

# For each element of nesting_levels_to_edit, name a shapefile with points
# where domains should be coarsened. The degree of coarsening will be defined
# in the SWALS multidomain_design_mod.f90. Empty string implies no coarsening
points_to_coarsen_shapefiles = c('', 'second_level_domains_to_coarsen', '')

stopifnot(length(nesting_levels_to_edit) == length(points_to_coarsen_shapefiles))

for(i in 1:length(nesting_levels_to_edit)){

    bte = read.csv(nesting_levels_to_edit[i])
    points_to_coarsen = points_to_coarsen_shapefiles[i]

    if(points_to_coarsen != ""){
        # Flag that some boxes should be coarsened
        sites_to_keep_fine = coordinates(readOGR(points_to_coarsen, layer=gsub('.shp', '', points_to_coarsen, fixed=TRUE)))

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

    output_filename = gsub('level_nesting.csv', 'level_nesting_edited.csv', nesting_levels_to_edit[i], fixed=TRUE)

    write.csv(bte, output_filename, row.names=FALSE)
}
