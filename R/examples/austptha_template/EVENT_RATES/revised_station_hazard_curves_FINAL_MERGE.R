#
# To do the FULL revised stage-percentile calculations (with all hazard points,
# all slip and rigidity models), the code
# 'revised_station_hazard_curves_FINAL.R' was run over 50 nodes by giving each
# node a 'chunk' of hazard points to work on (about 400 points/node). On each
# node the work was distributed over 16 cores.
#
# Here we re-combine the results to 'look like' it was all done on a single node,
# which makes it easier to process with variants of our old scripts.
#

merge_distributed_results<-function(source_name){

    # Files, each with a list containing chunks of points for this source.
    distributed_stage_pc_files = Sys.glob(paste0(
        './preprocessed_source_rate_revised_stage_exrates_FULL/stage_exceedance_rate_percentiles_',
        source_name, '_p_*.RDS'))

    # We want to sort the files so that the points will be correctly ordered.
    order_files = order(sapply(distributed_stage_pc_files, 
        f<-function(x){
            # With this split, the first point index will be in the second-last
            # entry of the filename
            filename_split = strsplit(x, '_')[[1]]
            my_first_point = as.numeric(filename_split[length(filename_split)-1])
            return(my_first_point)
        })
        )
    distributed_stage_pc_files = distributed_stage_pc_files[order_files] # Now the ordering is correct

    # Read them all
    all_distributed = lapply(distributed_stage_pc_files, f<-function(x) readRDS(x))

    # By combining the lists in this way, the result is 'just like' we would have got by 
    # processing all points at once
    merged_results = do.call(c, all_distributed)

    # Double-check that the ordering is correct (easy because each entry of the
    # list includes a "point_inds_range" variable giving the integer indices of
    # the points)
    point_ranges = unlist(lapply(merged_results, f<-function(x) x$point_inds_range))
    if(any(diff(point_ranges) < 1)) stop('ERROR in point ordering')

    # Save the results to a new folder, so we don't mix them up with the distributed points
    # (which would break the current algorithm, if we applied it to the directory with the new files)
    output_file = paste0(dirname(distributed_stage_pc_files[1]), '_MERGED/stage_exceedance_rate_percentiles_', 
                         source_name, '_p_1_', max(point_ranges), '.RDS')
    dir.create(dirname(output_file), showWarnings=FALSE)

    saveRDS(merged_results, file=output_file)

    return(0)
}

#
# Do the merger for all source-zones
#

# Get the source names -- here we do it by reading key filesnames and stripping out the source-name
all_source_name_RDS_inputs = Sys.glob('./preprocessed_source_rate_revised_stage_percentiles/preprocessed_rate_info_*.RDS')
tmp = basename(all_source_name_RDS_inputs)
tmp = gsub("preprocessed_rate_info_", "", tmp)
all_source_names = gsub(".RDS", "", tmp, fixed=TRUE)

for(source_name in all_source_names){
    merge_distributed_results(source_name)
}

