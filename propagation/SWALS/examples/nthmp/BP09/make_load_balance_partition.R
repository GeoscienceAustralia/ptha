#
# Script to help load balance multidomain
#
# When we split a multidomain in parallel, the individual sub-domains
# don't really take the same amount of time -- because they differ
# in wet/dryness (influencing some expensive "pow" calculations) and they
# differ in how much comms is needed. 
#
# This script figures how to assign sub-domains to images, so that when
# summed, the run-time is fairly similar. This might make the parallel run
# more efficient
#
# FIXME: Consider rearranging 'on-a-node'. First version didn't worry about
# this -- and I found that although the run times were more consistent, not
# so much waiting, and there was a 5% speedup, on the bad side the time in comms
# went up significantly. Maybe that's because we have more "inter-node" comms?
#

most_recent_run_dir = rev(sort(Sys.glob('OUTPUTS/RUN*')))[1]
print(paste0('Generating load_balance_partition.txt file in ', most_recent_run_dir))
print('To be valid, this assumes that job used coarrays without any load balancing.')
print(' Future jobs which use this file should have the same number of images and threads')

md_files = Sys.glob(paste0(most_recent_run_dir, '/multi*.log'))

md_runtime<-function(md_file){
    md_lines = readLines(md_file)

    domain_timer_starts = grep('Timer     ', md_lines)
    total_wallclock_lines = grep('Total WALLCLOCK', md_lines)

    #
    nd = length(domain_timer_starts)
    total_wallclock = sapply(md_lines[total_wallclock_lines[1:nd]], f<-function(x) as.numeric(strsplit(x, ':', fixed=TRUE)[[1]][2]))
 
    return(as.numeric(total_wallclock))
}
md_times = sapply(md_files, md_runtime)

# 
image_i_domains = matrix(NA, nrow=nrow(md_times), ncol=ncol(md_times))

mean_domain_times = apply(md_times, 1, mean)

counter = 0
for(i in rev(order(mean_domain_times))){
    counter = counter+1
    if(counter == 1){
        # Assign the largest domain in the same order
        image_i_domains[i,] = 1:ncol(image_i_domains)
        total_time = md_times[i, image_i_domains[i,]]
    }else{
        # Assign domains to images based on their "total time" so far
        rank_next_domains = rank(md_times[i,], ties='random')
        rank_total_times = rank(-total_time, ties='random')
        image_i_domains[i,] = match(rank_total_times, rank_next_domains)
        total_time = total_time + md_times[i, image_i_domains[i,]]
    }

    if(length(unique(image_i_domains[i,])) != length(image_i_domains[i,])){
        stop('Some domains are missing / double ups')
    }
}

previous_sum_times = colSums(md_times)
new_sum_times = total_time

print(paste0('Previous time range: ', diff(range(previous_sum_times))))
print(paste0('  Previous ', rev(c('min', 'max')), ' times : ', rev(range(previous_sum_times))))
print(paste0('New time range: ', diff(range(new_sum_times))))
print(paste0('  New ', rev(c('min', 'max')), ' times : ', rev(range(new_sum_times))))
print(paste0('Theoretical time reduction: ', max(previous_sum_times) - max(new_sum_times)))

write.table(image_i_domains, paste0(most_recent_run_dir, '/load_balance_partition.txt'), 
            sep=" ", row.names=FALSE, col.names=FALSE)

