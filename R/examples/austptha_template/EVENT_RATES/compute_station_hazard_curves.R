#
# Script to compute stage-vs-exceedance-rate curves for all gauges, for all
# source-zones.
#

library(rptha)
library(parallel)

# Get the config variables
config = new.env()
source('config.R', local=config)

#
# INPUTS
#

# By default run all source zones. But if the following variable is FALSE,
# then only run those for which the function check_if_file_is_ok (below)
# returns FALSE
run_all_source_zones = TRUE

# NetCDF files with uniform slip max_stage for every point, and also event
# rates
all_source_uniform_slip_tsunami = config$all_source_uniform_slip_tsunami

# NetCDF files with stochastic slip max_stage for every point, and also event
# rates
all_source_stochastic_slip_tsunami = config$all_source_stochastic_slip_tsunami 
# ... and variable_uniform
all_source_variable_uniform_slip_tsunami =
    config$all_source_variable_uniform_slip_tsunami 

# Apply hazard curve computation / data extraction to chunks of data
# We read in (nevents x point_chunk_size) max stage values at once, and then
# compute rate curves for the point_chunk_size gauges before moving to the next
# chunk. Note these values will be scaled later, so that they apply to the
# source-zone with the most events
point_chunk_size_uniform = config$point_chunk_size_uniform #10000
point_chunk_size_stochastic = config$point_chunk_size_stochastic #1000

# Sequence of stages at which we compute the rate, for every point, for every
# source-zone
stage_seq = config$stage_seq

# Number of cores to use in shared memory parallel
MC_CORES = config$MC_CORES

#
# END INPUTS
#

stage_seq_len = length(config$stage_seq) # Number of points defining rate curve

mycluster = makeForkCluster(nnodes=MC_CORES)

#'
#' Get lon/lat/elevation/id for all gauges
#'
#' @param nc_file netcdf file with gauge lon/lat/elev/gaugeID
#' @return data.frame with all variables
#'
read_lon_lat_elev<-function(nc_file){

    fid = nc_open(nc_file, readunlim=FALSE)
   
    lon = ncvar_get(fid, 'lon') 
    lat = ncvar_get(fid, 'lat') 
    elev = ncvar_get(fid, 'elev')
    gaugeID = ncvar_get(fid, 'gaugeID')

    nc_close(fid)

    output = data.frame(lon=lon, lat=lat, elev=elev, gaugeID = gaugeID)

    return(output)
}

#
# For the output files, check that rows that begin with NA also end with NA.
# An early version of the code had a bug here, which only affected small source-zones.
# Using this function I could correct the bug only in affected files, which saved
# a lot of run time
#
check_if_file_is_ok<-function(uniform_slip_tsunami_file){

    ###   #
    ###   # VERSION WHICH CHECKS FOR UN-WRITTEN RATES
    ###   #

    ###   # Make name for output file
    ###   source_name = basename(dirname(dirname(uniform_slip_tsunami_file)))
    ###   sourcename_dot_nc = paste0(source_name, '.nc')
    ###   nc_file = paste0(
    ###       dirname(uniform_slip_tsunami_file), '/',
    ###       'tsunami_stage_exceedance_rates_', sourcename_dot_nc)

    ###   if(!file.exists(nc_file)) return(FALSE)

    ###   fid = nc_open(nc_file, readunlim=FALSE)
    ###   rts = ncvar_get(fid, 'variable_mu_stochastic_slip_rate')
    ###   k1 = which(!is.finite(rts[1,]))
    ###   k2 = which(!is.finite(rts[100,]))
    ###   nc_close(fid)
    ###   if(setequal(k1, k2)){
    ###       file_ok = TRUE
    ###   }else{
    ###       file_ok = FALSE
    ###   }
    ###   return(file_ok)

    #
    # In this version, the user passes an integer corresponding to two words in 'subsets' below
    # e.g. 2 --> c('k', 'manokwari0')
    # Then if the source-name is alphabetically 'between' these words it is ran, otherwise it is not.
    #

    # These subset ranges were chosen to make the job time < 12 hours in each case.
    subsets = list(
        c('a', 'k'),
        c('k', 'manokwari0'),
        c('manokwari0', 'p'),
        c('p', 'south'), 
        c('south', 'sun'), # Only southamerica
        c('sun', 'z') )

    subset_to_run = as.numeric(commandArgs(trailingOnly=TRUE))
    stopifnot( 
        (length(subset_to_run) == 1) & 
        (subset_to_run > 0) & 
        (subset_to_run <= length(subsets)) )
    word_range = subsets[[subset_to_run]]

    # Get the source-name
    source_name = basename(dirname(dirname(uniform_slip_tsunami_file)))
    source_name = tolower(source_name)

    if(word_range[1] < source_name & word_range[2] > source_name){
        # Run these ones
        skip_source = FALSE 
    }else{
        skip_source = TRUE
    }
    return(skip_source)
 
}

#'
#' Get stage exceedance rates based on a file with tsunami peak-stages and rates
#' for each event
#' 
#' @param tsunami_file netcdf filename with tsunami max_stage / rates, for
#'     every gauge and event
#' @param gauge_points data.frame with gauge coordinates and elevation
#' @param point_chunk_size Number of gauges for which we read max_stage
#'     simultaneously
#' @param stage_seq stages at which we return the integrated exceedance rate
#' @return A list with 3 arrays, each having dim=c(length(stage_seq),
#'   length(gauge_points[,1])). The 3 arrays give the rate values(1), and
#'   upper(2) and lower(3) credible intervals for the rates, for each stage in
#'   stage_seq, for each gauge.
#'
source_zone_stage_exceedance_rates<-function(
    tsunami_file, 
    gauge_points, 
    point_chunk_size, 
    stage_seq){

    # Check to see if all the rates are zero. If they are, then make a quick exit
    fid_rates = nc_open(gsub('_tsunami', '', tsunami_file), readunlim=FALSE)
    event_rate = ncvar_get(fid_rates, 'rate_annual')
    variable_mu_event_rate = ncvar_get(fid_rates, 'variable_mu_rate_annual')
    nc_close(fid_rates)

    if(all(event_rate == 0) & all(variable_mu_event_rate == 0)){
        # Quick exit without using too much memory
        output = source_zone_stage_exceedance_rates_for_zero_rate_sources(
            tsunami_file, 
            gauge_points,
            point_chunk_size,
            stage_seq)
    }else{
        # Standard case
        output = source_zone_stage_exceedance_rates_standard(
            tsunami_file, 
            gauge_points,
            point_chunk_size,
            stage_seq)
    }

    return(output)
}

#'
#' Get stage exceedance rates based on a file with tsunami peak-stages and rates
#' for each event
#'
#' This is the main workhorse function
#' 
#' @param tsunami_file netcdf filename with tsunami max_stage / rates, for
#'     every gauge and event
#' @param gauge_points data.frame with gauge coordinates and elevation
#' @param point_chunk_size Number of gauges for which we read max_stage
#'     simultaneously
#' @param stage_seq stages at which we return the integrated exceedance rate
#' @return A list with 3 arrays, each having dim=c(length(stage_seq),
#'   length(gauge_points[,1])). The 3 arrays give the rate values(1), and
#'   upper(2) and lower(3) credible intervals for the rates, for each stage in
#'   stage_seq, for each gauge.
#'
source_zone_stage_exceedance_rates_standard<-function(
    tsunami_file, 
    gauge_points, 
    point_chunk_size, 
    stage_seq){

    lgp = length(gauge_points[,1])

    # Store the rates at every stage in stage_seq
    output_rates = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_lower_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_upper_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_median = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_16pc = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    output_rates_84pc = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates_lower_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates_upper_ci = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates_median = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates_16pc = matrix( NA, nrow=stage_seq_len, ncol=lgp )
    variable_mu_output_rates_84pc = matrix( NA, nrow=stage_seq_len, ncol=lgp )

    # Process gauges in chunks to reduce memory usage
    point_chunks = splitIndices(lgp, ceiling(lgp/point_chunk_size))

    # Get rate for every event
    fid_rates = nc_open(gsub('_tsunami', '', tsunami_file), readunlim=FALSE)
    event_rate = ncvar_get(fid_rates, 'rate_annual')
    event_rate_lower = ncvar_get(fid_rates, 'rate_annual_lower_ci')
    event_rate_upper = ncvar_get(fid_rates, 'rate_annual_upper_ci')
    event_rate_median = ncvar_get(fid_rates, 'rate_annual_median')
    event_rate_16pc = ncvar_get(fid_rates, 'rate_annual_16pc')
    event_rate_84pc = ncvar_get(fid_rates, 'rate_annual_84pc')

    variable_mu_event_rate = ncvar_get(fid_rates, 'variable_mu_rate_annual')
    variable_mu_event_rate_lower = ncvar_get(fid_rates, 'variable_mu_rate_annual_lower_ci')
    variable_mu_event_rate_upper = ncvar_get(fid_rates, 'variable_mu_rate_annual_upper_ci')
    variable_mu_event_rate_median = ncvar_get(fid_rates, 'variable_mu_rate_annual_median')
    variable_mu_event_rate_16pc = ncvar_get(fid_rates, 'variable_mu_rate_annual_16pc')
    variable_mu_event_rate_84pc = ncvar_get(fid_rates, 'variable_mu_rate_annual_84pc')

    nc_close(fid_rates)

    # File we read tsunami info from
    fid = nc_open(tsunami_file, readunlim=FALSE)

    # Do the calculation for all points, in chunks
    for(j in 1:length(point_chunks)){

        gauge_indices = point_chunks[[j]]
        # Ensure indices are contiguous
        stopifnot(all(range(diff(gauge_indices)) == 1))

        # Get max stage for every event, for this chunk of gauges
        peak_stages = ncvar_get(fid, 'max_stage', 
            start=c( 1, gauge_indices[1]), 
            count=c(-1, length(gauge_indices)))

        # Function for the inner loop which computes rates for the k'th gauge
        # in this chunk. Will run in parallel.
        rates_kth_gauge<-function(peak_stage_kth_gauge, 
            stage_seq = stage_seq, 
            event_rate=event_rate, 
            event_rate_lower=event_rate_lower, 
            event_rate_upper=event_rate_upper,
            event_rate_median=event_rate_median,
            event_rate_16pc=event_rate_16pc,
            event_rate_84pc=event_rate_84pc,
            variable_mu_event_rate=variable_mu_event_rate, 
            variable_mu_event_rate_lower=variable_mu_event_rate_lower, 
            variable_mu_event_rate_upper=variable_mu_event_rate_upper,
            variable_mu_event_rate_median=variable_mu_event_rate_median,
            variable_mu_event_rate_16pc=variable_mu_event_rate_16pc,
            variable_mu_event_rate_84pc=variable_mu_event_rate_84pc){

            if(all(is.na(peak_stage_kth_gauge))){
               # Skip dry points
               return(list(stage_seq*NA, stage_seq*NA, stage_seq*NA, 
                   stage_seq*NA, stage_seq*NA, stage_seq*NA,
                   stage_seq*NA, stage_seq*NA, stage_seq*NA,
                   stage_seq*NA, stage_seq*NA, stage_seq*NA))
            } 

            # Sort the stages in decreasing order 
            events_sort = sort(peak_stage_kth_gauge, index.return=TRUE, 
                decreasing=TRUE)

            # Value that no stage is larger than for this gauge
            stage_upper_bound = events_sort$x[1] + 1.0e-03
            # Value no stage is larger than for ANY gauge!
            global_stage_upper_bound = 99999999

            # Sort the rates like the stages, and do cumulative sum to get an
            # exceedance rate
            sorted_stages = c(global_stage_upper_bound, stage_upper_bound, 
                events_sort$x)
            sorted_cumulative_rate = cumsum(c(0, 0, 
                event_rate[events_sort$ix]))
            sorted_cumulative_rate_lower = cumsum(c(0, 0, 
                event_rate_lower[events_sort$ix]))
            sorted_cumulative_rate_upper = cumsum(c(0, 0, 
                event_rate_upper[events_sort$ix]))
            sorted_cumulative_rate_median = cumsum(c(0, 0, 
                event_rate_median[events_sort$ix]))
            sorted_cumulative_rate_16pc = cumsum(c(0, 0, 
                event_rate_16pc[events_sort$ix]))
            sorted_cumulative_rate_84pc = cumsum(c(0, 0, 
                event_rate_84pc[events_sort$ix]))
            

            variable_mu_sorted_cumulative_rate = cumsum(c(0, 0, 
                variable_mu_event_rate[events_sort$ix]))
            variable_mu_sorted_cumulative_rate_lower = cumsum(c(0, 0, 
                variable_mu_event_rate_lower[events_sort$ix]))
            variable_mu_sorted_cumulative_rate_upper = cumsum(c(0, 0, 
                variable_mu_event_rate_upper[events_sort$ix]))
            variable_mu_sorted_cumulative_rate_median = cumsum(c(0, 0, 
                variable_mu_event_rate_median[events_sort$ix]))
            variable_mu_sorted_cumulative_rate_16pc = cumsum(c(0, 0, 
                variable_mu_event_rate_16pc[events_sort$ix]))
            variable_mu_sorted_cumulative_rate_84pc = cumsum(c(0, 0, 
                variable_mu_event_rate_84pc[events_sort$ix]))

            output = vector(mode='list', length=12)

            # Rates. Pass rule=2:1 to approx, so that 
            # for stages < min(sorted_stages) {= sorted_stages[length(sorted_stages)]},
            # we return the total event rate. 
            # The case stages > max(sorted_stages) {= sorted_stages[1]}, should
            # never happen. The former case is possible for hazard-points 'on
            # or very near' the source-zone, for small source-zones, which e.g.
            # always have peak_stage > 2cm, for all events
            rate_fun = approx(sorted_stages, sorted_cumulative_rate, 
                xout=stage_seq, rule=2:1)
            output[[1]] = rate_fun$y

            # Upper estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_upper, 
                xout=stage_seq, rule=2:1)
            output[[2]] = rate_fun$y

            # Lower estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_lower, 
                xout=stage_seq, rule=2:1)
            output[[3]] = rate_fun$y

            # Median estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_median, 
                xout=stage_seq, rule=2:1)
            output[[4]] = rate_fun$y

            # 16pc estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_16pc, 
                xout=stage_seq, rule=2:1)
            output[[5]] = rate_fun$y

            # 84pc estimate of rates
            rate_fun = approx(sorted_stages, sorted_cumulative_rate_84pc, 
                xout=stage_seq, rule=2:1)
            output[[6]] = rate_fun$y

            # Rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate, 
                xout=stage_seq, rule=2:1)
            output[[7]] = rate_fun$y

            # Upper estimate of rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate_upper, 
                xout=stage_seq, rule=2:1)
            output[[8]] = rate_fun$y

            # Lower estimate of rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate_lower, 
                xout=stage_seq, rule=2:1)
            output[[9]] = rate_fun$y

            # Median estimate of rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate_median, 
                xout=stage_seq, rule=2:1)
            output[[10]] = rate_fun$y

            # 16pc estimate of rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate_16pc, 
                xout=stage_seq, rule=2:1)
            output[[11]] = rate_fun$y

            # 84pc estimate of rates
            rate_fun = approx(sorted_stages, variable_mu_sorted_cumulative_rate_84pc, 
                xout=stage_seq, rule=2:1)
            output[[12]] = rate_fun$y

            return(output)
        }

        # Main computation
        par_output = parCapply(cl=mycluster, 
            x=peak_stages, 
            FUN=rates_kth_gauge, 
            stage_seq=stage_seq, 
            event_rate=event_rate, 
            event_rate_upper=event_rate_upper, 
            event_rate_lower=event_rate_lower,
            event_rate_median=event_rate_median,
            event_rate_16pc=event_rate_16pc,
            event_rate_84pc=event_rate_84pc,
            variable_mu_event_rate=variable_mu_event_rate, 
            variable_mu_event_rate_upper=variable_mu_event_rate_upper, 
            variable_mu_event_rate_lower=variable_mu_event_rate_lower,
            variable_mu_event_rate_median=variable_mu_event_rate_median,
            variable_mu_event_rate_16pc=variable_mu_event_rate_16pc,
            variable_mu_event_rate_84pc=variable_mu_event_rate_84pc
            )

        # Unpack parallel output to main arrays
        for(i in 1:length(par_output)){
            gi = gauge_indices[i]
            output_rates[,gi] = par_output[[i]][[1]]
            output_rates_upper_ci[,gi] = par_output[[i]][[2]]
            output_rates_lower_ci[,gi] = par_output[[i]][[3]]
            output_rates_median[,gi] = par_output[[i]][[4]]
            output_rates_16pc[,gi] = par_output[[i]][[5]]
            output_rates_84pc[,gi] = par_output[[i]][[6]]

            variable_mu_output_rates[,gi] = par_output[[i]][[7]]
            variable_mu_output_rates_upper_ci[,gi] = par_output[[i]][[8]]
            variable_mu_output_rates_lower_ci[,gi] = par_output[[i]][[9]]
            variable_mu_output_rates_median[,gi] = par_output[[i]][[10]]
            variable_mu_output_rates_16pc[,gi] = par_output[[i]][[11]]
            variable_mu_output_rates_84pc[,gi] = par_output[[i]][[12]]
        }

        rm(par_output)
    }

    nc_close(fid)

    output = list(
        rates = output_rates, 
        rates_upper_ci = output_rates_upper_ci, 
        rates_lower_ci = output_rates_lower_ci,
        rates_median = output_rates_median,
        rates_16pc = output_rates_16pc,
        rates_84pc = output_rates_84pc,
        variable_mu_rates = variable_mu_output_rates, 
        variable_mu_rates_upper_ci = variable_mu_output_rates_upper_ci, 
        variable_mu_rates_lower_ci = variable_mu_output_rates_lower_ci,
        variable_mu_rates_median = variable_mu_output_rates_median,
        variable_mu_rates_16pc = variable_mu_output_rates_16pc,
        variable_mu_rates_84pc = variable_mu_output_rates_84pc)

    return(output)
}

#' This is a 'quick-exit' version of the function above
#'
#'
source_zone_stage_exceedance_rates_for_zero_rate_sources<-function(
    tsunami_file, 
    gauge_points, 
    point_chunk_size, 
    stage_seq){

    lgp = length(gauge_points[,1])
    zero_rate_matrix = matrix( 0, nrow=stage_seq_len, ncol=lgp )

    output = list(
        rates = zero_rate_matrix, 
        rates_upper_ci = zero_rate_matrix, 
        rates_lower_ci = zero_rate_matrix,
        rates_median = zero_rate_matrix,
        rates_16pc = zero_rate_matrix,
        rates_84pc = zero_rate_matrix,
        variable_mu_rates = zero_rate_matrix, 
        variable_mu_rates_upper_ci = zero_rate_matrix, 
        variable_mu_rates_lower_ci = zero_rate_matrix,
        variable_mu_rates_median = zero_rate_matrix,
        variable_mu_rates_16pc = zero_rate_matrix,
        variable_mu_rates_84pc = zero_rate_matrix)

    return(output)
}

#'
#' Take care of saving outputs to netcdf file
#'
#'
create_rate_netcdf_file<-function(
    source_name, 
    gauge_points, 
    stage_seq, 
    uniform_slip_rates, 
    stochastic_slip_rates,
    variable_uniform_slip_rates,
    uniform_slip_tsunami_file, 
    stochastic_slip_tsunami_file,
    variable_uniform_slip_tsunami_file){

    # Dimension for rate curve
    dim_stage_seq = ncdim_def('stage', 'm', vals=stage_seq, unlim=FALSE,
        longname='stages corresponding to tsunami wave height exceedance rates')

    # Dimension for gauges
    dim_station = ncdim_def('station', '', vals=1:length(gauge_points[,1]), 
        unlim=TRUE,
        longname='integer index corresponding to the gauge location')

    # Variables for gauge locations
    gauge_lon_v = ncvar_def(name='lon', units='degrees_east', 
        dim=list(dim_station), missval=NA, longname='station_longitude', 
        prec='float')
    gauge_lat_v = ncvar_def(name='lat', units='degrees_north', 
        dim=list(dim_station), missval=NA, longname='station_latitude', 
        prec='float')
    gauge_elev_v = ncvar_def(name='elev', units='m', dim=list(dim_station), 
        missval=NA, longname='station_ground_elevation_above_mean_sea_level', 
        prec='float')
    gauge_id_v = ncvar_def(name='gaugeID', units='', dim=list(dim_station), 
        missval=NA, longname='real_ID_for_each_station', prec='float')

    all_nc_var = list(gauge_lon_v, gauge_lat_v, gauge_elev_v, gauge_id_v)

    for(slip_type in c('uniform', 'stochastic', 'variable_uniform')){
        for(vary_mu in c('', 'variable_mu_')){

            #
            # Regular rate
            #
            var_name1 = paste0(vary_mu, slip_type, '_rate_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate')

            if(vary_mu==''){
                extra_longname = ''
            }else{
                extra_longname = ' with variable shear modulus'
            }

            assign(var_name1, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate of peak stage for ', slip_type, ' slip events ', extra_longname),
                missval=NA,
                prec='float'))

            ## Variables for rates, uniform slip
            #uniform_rate_v = ncvar_def(
            #    name='uniform_slip_rate', units='events per year',
            #    dim=list(dim_stage_seq, dim_station), 
            #    longname = 'exceedance rate of peak stage for uniform slip events',
            #    missval=NA,
            #    prec='float')

            #
            # Upper credible interval
            # 
            var_name2 = paste0(vary_mu, slip_type, '_rate_upper_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate_upper_ci')

            assign(var_name2, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate (upper credible interval) of peak stage for ', slip_type, 
                    ' slip events ', extra_longname),
                missval=NA,
                prec='float'))

            #uniform_rate_upper_v = ncvar_def(
            #    name='uniform_slip_rate_upper_ci', units='events per year',
            #    dim=list(dim_stage_seq, dim_station), 
            #    longname = 'exceedance rate (upper credible interval) of peak stage for uniform slip events',
            #    missval=NA,
            #    prec='float')

            #
            # Lower credible interval
            # 
            var_name3 = paste0(vary_mu, slip_type, '_rate_lower_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate_lower_ci')

            assign(var_name3, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate (lower credible interval) of peak stage for ', slip_type, 
                    ' slip events ', extra_longname),
                missval=NA,
                prec='float'))

            #uniform_rate_lower_v = ncvar_def(
            #    name='uniform_slip_rate_lower_ci', units='events per year',
            #    dim=list(dim_stage_seq, dim_station), 
            #    longname = 'exceedance rate (lower credible interval) of peak stage for uniform slip events',
            #    missval=NA, prec='float')

            #
            # Median credible interval
            # 
            var_name4 = paste0(vary_mu, slip_type, '_rate_median_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate_median')

            assign(var_name4, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate (median) of peak stage for ', slip_type, 
                    ' slip events ', extra_longname),
                missval=NA,
                prec='float'))
            #
            # 16th percentile
            # 
            var_name5 = paste0(vary_mu, slip_type, '_rate_16pc_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate_16pc')

            assign(var_name5, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate (16th percentile) of peak stage for ', slip_type, 
                    ' slip events ', extra_longname),
                missval=NA,
                prec='float'))

            #
            # 84th percentile
            #

            var_name6 = paste0(vary_mu, slip_type, '_rate_84pc_v')
            var_title = paste0(vary_mu, slip_type, '_slip_rate_84pc')

            assign(var_name6, ncvar_def(
                name=var_title, units='events per year',
                dim=list(dim_stage_seq, dim_station), 
                longname = paste0('exceedance rate (84th percentile) of peak stage for ', slip_type, 
                    ' slip events ', extra_longname),
                missval=NA,
                prec='float'))

            all_nc_var = c(all_nc_var, list(
                eval(as.name(var_name1)), eval(as.name(var_name2)), eval(as.name(var_name3)), 
                eval(as.name(var_name4)), eval(as.name(var_name5)), eval(as.name(var_name6))
                ))
        }
    }
    

    # Make name for output file
    sourcename_dot_nc = paste0(source_name, '.nc')
    output_file_name = paste0(
        dirname(uniform_slip_tsunami_file), '/', 
        'tsunami_stage_exceedance_rates_', sourcename_dot_nc)

    # Create output file
    output_fid = nc_create(output_file_name, vars=all_nc_var)

    # Put attributes on file
    ncatt_put(output_fid, varid=0, attname = 'uniform_slip_tsunami_event_file',
        attval=normalizePath(uniform_slip_tsunami_file), prec='text')
    ncatt_put(output_fid, varid=0, 
        attname = 'stochastic_slip_tsunami_event_file',
        attval=normalizePath(stochastic_slip_tsunami_file), prec='text')
    ncatt_put(output_fid, varid=0, 
        attname = 'variable_uniform_slip_tsunami_event_file',
        attval=normalizePath(variable_uniform_slip_tsunami_file), prec='text')
    ncatt_put(output_fid, varid=0, attname='source_zone_name',
        attval=source_name, prec='text')
    ncatt_put(output_fid, varid=0, attname='parent_script_name',
        attval=parent_script_name(), prec='text')

    # Put gauge info on file
    ncvar_put(output_fid, gauge_lon_v, gauge_points$lon)
    ncvar_put(output_fid, gauge_lat_v, gauge_points$lat)
    ncvar_put(output_fid, gauge_elev_v, gauge_points$elev)
    ncvar_put(output_fid, gauge_id_v, gauge_points$gaugeID)

    # Put uniform slip stage exceedance rates on file
    ncvar_put(output_fid, uniform_rate_v,        uniform_slip_rates$rates)
    ncvar_put(output_fid, uniform_rate_upper_v,  uniform_slip_rates$rates_upper_ci)
    ncvar_put(output_fid, uniform_rate_lower_v,  uniform_slip_rates$rates_lower_ci)
    ncvar_put(output_fid, uniform_rate_median_v, uniform_slip_rates$rates_median)
    ncvar_put(output_fid, uniform_rate_16pc_v,   uniform_slip_rates$rates_16pc)
    ncvar_put(output_fid, uniform_rate_84pc_v,   uniform_slip_rates$rates_84pc)

    # As above with variable shear modulus
    ncvar_put(output_fid, variable_mu_uniform_rate_v,        uniform_slip_rates$variable_mu_rates)
    ncvar_put(output_fid, variable_mu_uniform_rate_upper_v,  uniform_slip_rates$variable_mu_rates_upper_ci)
    ncvar_put(output_fid, variable_mu_uniform_rate_lower_v,  uniform_slip_rates$variable_mu_rates_lower_ci)
    ncvar_put(output_fid, variable_mu_uniform_rate_median_v, uniform_slip_rates$variable_mu_rates_median)
    ncvar_put(output_fid, variable_mu_uniform_rate_16pc_v,   uniform_slip_rates$variable_mu_rates_16pc)
    ncvar_put(output_fid, variable_mu_uniform_rate_84pc_v,   uniform_slip_rates$variable_mu_rates_84pc)

    # Put stochastic slip stage exceedance rates on file
    ncvar_put(output_fid, stochastic_rate_v,        stochastic_slip_rates$rates)
    ncvar_put(output_fid, stochastic_rate_upper_v,  stochastic_slip_rates$rates_upper_ci)
    ncvar_put(output_fid, stochastic_rate_lower_v,  stochastic_slip_rates$rates_lower_ci)
    ncvar_put(output_fid, stochastic_rate_median_v, stochastic_slip_rates$rates_median)
    ncvar_put(output_fid, stochastic_rate_16pc_v,   stochastic_slip_rates$rates_16pc)
    ncvar_put(output_fid, stochastic_rate_84pc_v,   stochastic_slip_rates$rates_84pc)

    # As above with variable shear modulus
    ncvar_put(output_fid, variable_mu_stochastic_rate_v,        stochastic_slip_rates$variable_mu_rates)
    ncvar_put(output_fid, variable_mu_stochastic_rate_upper_v,  stochastic_slip_rates$variable_mu_rates_upper_ci)
    ncvar_put(output_fid, variable_mu_stochastic_rate_lower_v,  stochastic_slip_rates$variable_mu_rates_lower_ci)
    ncvar_put(output_fid, variable_mu_stochastic_rate_median_v, stochastic_slip_rates$variable_mu_rates_median)
    ncvar_put(output_fid, variable_mu_stochastic_rate_16pc_v,   stochastic_slip_rates$variable_mu_rates_16pc)
    ncvar_put(output_fid, variable_mu_stochastic_rate_84pc_v,   stochastic_slip_rates$variable_mu_rates_84pc)

    # Put variable_uniform slip stage exceedance rates on file
    ncvar_put(output_fid, variable_uniform_rate_v,        variable_uniform_slip_rates$rates)
    ncvar_put(output_fid, variable_uniform_rate_upper_v,  variable_uniform_slip_rates$rates_upper_ci)
    ncvar_put(output_fid, variable_uniform_rate_lower_v,  variable_uniform_slip_rates$rates_lower_ci)
    ncvar_put(output_fid, variable_uniform_rate_median_v, variable_uniform_slip_rates$rates_median)
    ncvar_put(output_fid, variable_uniform_rate_16pc_v,   variable_uniform_slip_rates$rates_16pc)
    ncvar_put(output_fid, variable_uniform_rate_84pc_v,   variable_uniform_slip_rates$rates_84pc)

    # As above with variable shear modulus
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_v,        variable_uniform_slip_rates$variable_mu_rates)
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_upper_v,  variable_uniform_slip_rates$variable_mu_rates_upper_ci)
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_lower_v,  variable_uniform_slip_rates$variable_mu_rates_lower_ci)
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_median_v, variable_uniform_slip_rates$variable_mu_rates_median)
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_16pc_v,   variable_uniform_slip_rates$variable_mu_rates_16pc)
    ncvar_put(output_fid, variable_mu_variable_uniform_rate_84pc_v,   variable_uniform_slip_rates$variable_mu_rates_84pc)

    nc_close(output_fid)

    return(invisible(output_file_name))
}

# Get point info
gauge_points = read_lon_lat_elev(all_source_uniform_slip_tsunami[1])

# Names of sources
source_names = basename(dirname(dirname(all_source_uniform_slip_tsunami)))

# Double-check that uniform/stochastic/variable_uniform files are both ordered
# by source_name 
stopifnot(all(source_names == basename(dirname(dirname(all_source_stochastic_slip_tsunami)))))
stopifnot(all(source_names == basename(dirname(dirname(all_source_variable_uniform_slip_tsunami)))))

#
# Useful to know the number of events in the netcdf files
#
get_number_of_events<-function(nc_file){
    fid = nc_open(nc_file, readunlim=FALSE)
    output = fid$dim$event$len
    nc_close(fid)
    return(output)
}
uniform_slip_nevents = sapply(all_source_uniform_slip_tsunami, get_number_of_events)
stochastic_slip_nevents = sapply(all_source_stochastic_slip_tsunami, get_number_of_events)
variable_uniform_slip_nevents = sapply(all_source_variable_uniform_slip_tsunami, get_number_of_events)

# For parallel chunking, we can use larger chunks if there are fewer events.
uniform_chunk_size = round(point_chunk_size_uniform * 
    max(uniform_slip_nevents)/uniform_slip_nevents)
stochastic_chunk_size = round(point_chunk_size_stochastic * 
    max(stochastic_slip_nevents)/stochastic_slip_nevents)
# Note we use 'point_chunk_size_stochastic' here since variable uniform has the
# same number of events as stochastic
variable_uniform_chunk_size = round(point_chunk_size_stochastic * 
    max(variable_uniform_slip_nevents)/variable_uniform_slip_nevents)

for(i in 1:length(source_names)){

    # source information
    source_name = source_names[i]
    uniform_slip_tsunami_file = all_source_uniform_slip_tsunami[i]
    stochastic_slip_tsunami_file = all_source_stochastic_slip_tsunami[i]
    variable_uniform_slip_tsunami_file = all_source_variable_uniform_slip_tsunami[i]

    if(!run_all_source_zones){
        # Call a function to check if the file is ok. If it is ok, then go to
        # the next source.
        is_done = check_if_file_is_ok(uniform_slip_tsunami_file)
        if(is_done) next
    }

    # Get uniform slip outputs
    uniform_slip_rates = source_zone_stage_exceedance_rates(
        uniform_slip_tsunami_file, gauge_points,
        point_chunk_size=uniform_chunk_size[i], stage_seq=stage_seq)

    # Remove memory from cluster
    parLapply(mycluster, as.list(1:MC_CORES), gc)

    # Get stochastic slip outputs
    stochastic_slip_rates = source_zone_stage_exceedance_rates(
        stochastic_slip_tsunami_file, gauge_points,
        point_chunk_size=stochastic_chunk_size[i], stage_seq=stage_seq)

    # Remove memory from cluster
    parLapply(mycluster, as.list(1:MC_CORES), gc)

    # Get variable_uniform slip outputs
    variable_uniform_slip_rates = source_zone_stage_exceedance_rates(
        variable_uniform_slip_tsunami_file, gauge_points,
        point_chunk_size=variable_uniform_chunk_size[i], stage_seq=stage_seq)

    # Remove memory from cluster
    parLapply(mycluster, as.list(1:MC_CORES), gc)

    # Write out to a file 'tsunami_stage_exceedance_rates_SOURCENAME.nc'
    # in the same folder as uniform_slip_tsunami_file
    create_rate_netcdf_file(
        source_name, gauge_points, stage_seq,
        uniform_slip_rates, 
        stochastic_slip_rates, 
        variable_uniform_slip_rates,
        uniform_slip_tsunami_file, 
        stochastic_slip_tsunami_file,
        variable_uniform_slip_tsunami_file)

    # Save some memory
    rm(uniform_slip_rates, stochastic_slip_rates, variable_uniform_slip_rates)
    gc()
}
stopCluster(mycluster)
