#
# Run this from inside a multidomain directory to compute a load_balance file,
# and report on expected time-savings. If the results are good you can them manually
# copy that file to 
#

mydir = getwd()
if(!grepl('RUN_', basename(mydir))){
    stop('This script must be run from inside a multidomain directory')
}

source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')

create_domain_index_groups<-function(){
    # Create domain index groups. Group domain by
    # - min_timestep & is_dispersive

    all_md_logs = Sys.glob('multi*.log')

    domain_metadata = vector(mode='list', length=length(all_md_logs)) # Storage space

    for(i in 1:length(all_md_logs)){
        log_lines = readLines(all_md_logs[i])
        # Find lines that describe domain index, min_timestep, and whether it is dispersive
        i0 = grep("Gauges are setup", log_lines); stopifnot(length(i0) == 1)
        i1 = grep("End setup", log_lines); stopifnot(length(i1) == 1)
        text = log_lines[(i0+1):(i1-1)]

        # Find the beginning of each domain descriptor
        di_start = grep("domain:", text)
        di_end = c(di_start[-1]-1, length(text))
        nd = length(di_start)

        domain_metadata[[i]] = vector(mode='list', length=nd) # Storage space

        # Loop over the domains and extract their variables
        for(j in 1:nd){
            text_string = paste0(text[di_start[j]:di_end[j]], collapse=" ")
            text_lines = scan(text=text_string, sep=",", what='character', quiet=TRUE)

            # Extract the domain index
            k = grep('domain:', text_lines); stopifnot(length(k) == 1)
            tmp = strsplit(text_lines[k], split=" ")[[1]]
            domain_id = tmp[nchar(tmp) > 0][2]
            domain_index = as.numeric(substring(domain_id, 6, nchar(domain_id))) # Strip leading numbers (image)

            # Extract the local index. This will be 1, except on domains that were partitioned
            # by the load balancing (where there will be a unique integer for each piece)
            k = grep('local_index:', text_lines); stopifnot(length(k) == 1)
            tmp = strsplit(text_lines[k], split=" ")[[1]]
            domain_local_index = as.numeric(tmp[nchar(tmp) > 0][2])
            
            # Is the domain dispersive?
            domain_is_dispersive = !any(grepl("non-dispersive", text_lines))

            k = grep('min_timestep: ', text_lines); stopifnot(length(k) == 1)
            tmp = strsplit(text_lines[k], split=" ")[[1]]
            min_ts_string = tmp[nchar(tmp) > 0][2]
            domain_min_ts = as.numeric(min_ts_string)

            domain_metadata[[i]][[j]] = data.frame(index=domain_index, 
                local_index = domain_local_index, 
                is_dispersive = domain_is_dispersive, 
                min_ts=domain_min_ts)
        }
        domain_metadata[[i]] = do.call(rbind, domain_metadata[[i]])
    }

    combined_metadata = do.call(rbind, domain_metadata)
    o1 = order(combined_metadata$index)
    combined_metadata = combined_metadata[o1,]
   
    return(combined_metadata) 
    # The above can contain repeated entries (e.g. a piece of domain 1 on every mpi_rank)
}
domain_metadata = create_domain_index_groups()
if(!any(domain_metadata$is_dispersive)){
    stop('Error: Did not find any dispersive domains -- use "load_balance_script.R" instead of this script')
}
min_dispersive_timestep = min(domain_metadata$min_ts[domain_metadata$is_dispersive])
# Group domains based on dispersion and timestep, combining all NLSW domains
# with timestep <= min_dispersive_timestep
domain_metadata_group_ids = paste0(domain_metadata$is_dispersive, '_', 
    pmax(domain_metadata$min_ts, min_dispersive_timestep))

# Grouping
unique_ids = unique(domain_metadata_group_ids)
domain_index_groups = lapply(unique_ids, function(x){
    # For each unique_id, return a matrix with columns having
    # domain_index, domain_local_index
    k = domain_metadata_group_ids == x
    cbind(domain_metadata$index[k],
          domain_metadata$local_index[k])
    })
names(domain_index_groups) = unique_ids


#browser()

# Assume we want every processor to have part of domain 1. This makes sense
# if domain 1 is a global domain, and the other domains may utilise static_before_time.
x = make_load_balance_partition('.', domain_index_groups=domain_index_groups)

