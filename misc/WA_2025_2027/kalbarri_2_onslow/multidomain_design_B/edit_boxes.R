#'
#' Append additional information to the boxes that SWALS can use:
#' 1. Which boxes to coarsen
#'
library(rptha)

# Kalbarri to onslow + Xmas/Cocos
nesting_levels_to_edit = paste0(
    'kalbarri2onslow_domains_20251218_0.333333333333333_1_3_9_1_3_15/', 
    c(
        'kalbarri2onslow_C_level1.csv',
        'kalbarri2onslow_C_level2.csv',
        'kalbarri2onslow_C_level3.csv',
        'kalbarri2onslow_C_level4.csv',
        'Xmas_cocos_level1.csv',
        'Xmas_cocos_level2.csv',
        'Xmas_cocos_level3.csv'
    )
)
points_to_coarsen_shapefiles = rep("", length(nesting_levels_to_edit))
points_to_coarsen_layers = rep("", length(nesting_levels_to_edit))


#
# Shouldn't need changes below here
#

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
