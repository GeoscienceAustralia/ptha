#
# Append additional information to the boxes that SWALS can use
# 
boxes_to_edit = Sys.glob('domains_010322_0.5_0.166666666666667_0.0333333333333333/*level_nesting.csv')

for(i in 1:length(boxes_to_edit)){

    bte = read.csv(boxes_to_edit[i])


    if(basename(boxes_to_edit[i]) == 'second_level_nesting.csv'){
        # Flag that some boxes should be coarsened
        library(rgdal) 
        sites_to_keep_fine = coordinates(readOGR('keep_highres', layer='keep_highres'))

        # Initialise to 'all coarse'
        coarsen_box = rep(1, nrow(bte))

        # Refine boxes containing a 'keep_highres' point
        for(j in 1:nrow(sites_to_keep_fine)){
            has_a_point = (bte[,1] < sites_to_keep_fine[j,1]) & (bte[,3] > sites_to_keep_fine[j,1]) &
                          (bte[,2] < sites_to_keep_fine[j,2]) & (bte[,4] > sites_to_keep_fine[j,2])
            coarsen_box[has_a_point] = 0
        }
        
    }else{
        coarsen_box = rep(0, nrow(bte))
    }

    bte = cbind(bte, data.frame('coarsen' = coarsen_box))

    output_filename = gsub('level_nesting.csv', 'level_nesting_edited.csv', boxes_to_edit[i], fixed=TRUE)

    write.csv(bte, output_filename, row.names=FALSE)
}
