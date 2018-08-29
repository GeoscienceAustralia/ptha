#
# A few useful functions to check the sourcezone_parameter file
#

#'
#' Logical check on the row weights
#' Wrap in a function to avoid polluting the namespace
#'
check_sourcezone_parameter_row_weights<-function(sourcezone_parameters){

    full_source_names = sourcezone_parameters$sourcename
    segment_names = sourcezone_parameters$segment_name
    unique_source_names = unique(full_source_names)

    for(i in 1:length(unique_source_names)){

        k = which(full_source_names == unique_source_names[i])
        row_weights = sourcezone_parameters$row_weight[k]

        if(length(row_weights) > 1){
            # Segmented sources

            # Ensure the 'unsegmented' case is first
            if(segment_names[k[1]] != ''){
                print(paste0('Problem in input table for ', unique_source_names[i]))
                stop('The unsegmented source should come first, and have a blank entry in the segment_name column')
            }

            # We need the segment weights to be equal
            segment_weights = row_weights[-1]
            test1 = all(segment_weights == segment_weights[1])
            # We need the full-source-zone weight + segment_weight to be unity, unless they are all zero

            # Exit in the zero case
            if(all(row_weights == 0)) next

            # Non-zero segment weights
            test2 = all( abs(segment_weights + row_weights[1] - 1) < 1.0e-12 )
            if(!(test1 & test2)){
                print(paste0('Problem in input table for ', unique_source_names[i]))
                stop('row_weights for all segments must be equal, and equal (1 - the unsegmented weight), unless all are zero')
            }

        }else{
            # Unsegmented sources

            if( !(row_weights %in% c(0.0, 1.0)) ){
                print(paste0('Problem in input table for ', unique_source_names[i]))
                stop('Unsegmented row_weight must equal 1.0, or zero')
            }
        }
    }

}
