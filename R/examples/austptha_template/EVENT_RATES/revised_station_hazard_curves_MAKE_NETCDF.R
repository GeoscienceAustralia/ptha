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

# Files that contain the calculations we need
all_stage_calc_files = Sys.glob(
    './preprocessed_source_rate_revised_stage_exrates_FULL_MERGED/stage_exceedance_rate_percentiles_*.RDS')
# Source names for which we have the above files -- all others will have dummy 'zero rate' files.
SOURCE_NAMES_TO_KEEP = gsub('_p_1_20185.RDS', '',
    gsub('stage_exceedance_rate_percentiles_', '', basename(all_stage_calc_files)))

# NetCDF files with uniform slip max_stage for every point, and also event
# rates
all_source_uniform_slip_tsunami = config$all_source_uniform_slip_tsunami

# NetCDF files with stochastic slip max_stage for every point, and also event
# rates
all_source_stochastic_slip_tsunami = config$all_source_stochastic_slip_tsunami 
# ... and variable_uniform
all_source_variable_uniform_slip_tsunami =
    config$all_source_variable_uniform_slip_tsunami 

# Sequence of stages at which we compute the rate, for every point, for every
# source-zone
stage_seq = config$stage_seq

# Number of cores to use in shared memory parallel
#MC_CORES = config$MC_CORES

#
# END INPUTS
#

stage_seq_len = length(config$stage_seq) # Number of points defining rate curve

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

# Read the updated files, mean exceedance-rates
get_mean_exrates<-function(one_stage_calc, slip_type = 'stochastic', mu_type = 'fixed_mu'){
    # Make one list with stage_mean_exrate
    all_mean_exrates = lapply(one_stage_calc, 
        f<-function(x) x[[slip_type]][[mu_type]]$stage_mean_exrate)
    # Merge to a matrix
    mean_exrates = do.call(cbind, all_mean_exrates)
    return(mean_exrates)
}

# Read the updated files, percentile exceedance-rates
get_percentile_exrates<-function(one_stage_calc, pc, slip_type = 'stochastic', mu_type = 'fixed_mu'){
    row = switch(pc, '0.025' = 1, '0.16' = 2, '0.5' = 3, '0.84' = 4, '0.975' = 5)
    # Make one list with the desired percentile
    all_pc_exrates = lapply(one_stage_calc, 
        f<-function(x) x[[slip_type]][[mu_type]]$stage_percentile_exrates[row,,])
    # Merge to a matrix
    pc_exrates = do.call(cbind, all_pc_exrates)
    return(pc_exrates)
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
    source_name,
    slip_type,
    gauge_points){

    if(!(source_name %in% SOURCE_NAMES_TO_KEEP)){
        # Quick exit without using too much memory
        output = source_zone_stage_exceedance_rates_for_zero_rate_sources(
            source_name, slip_type, gauge_points)
    }else{
        # Standard case
        output = source_zone_stage_exceedance_rates_standard(
            source_name, slip_type, gauge_points)
    }

    return(output)
}

#'
#' Get stage exceedance rates 
#'
#' This is the main workhorse function
#' 
#'
source_zone_stage_exceedance_rates_standard<-function(
    source_name,
    slip_type,
    gauge_points){

    #
    # Get the results
    #
    my_file_index = grep(paste0('_', source_name, '_p_'), all_stage_calc_files)
    one_stage_calc = readRDS(all_stage_calc_files[my_file_index])

    # Mean rates
    output_rates             = get_mean_exrates(one_stage_calc, slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates = get_mean_exrates(one_stage_calc, slip_type = slip_type, mu_type = 'variable_mu')

    # Percentiles
    output_rates_upper_ci             = get_percentile_exrates(one_stage_calc, pc='0.975', slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates_upper_ci = get_percentile_exrates(one_stage_calc, pc='0.975', slip_type = slip_type, mu_type = 'variable_mu')
    output_rates_84pc                 = get_percentile_exrates(one_stage_calc, pc= '0.84', slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates_84pc     = get_percentile_exrates(one_stage_calc, pc= '0.84', slip_type = slip_type, mu_type = 'variable_mu')
    output_rates_median               = get_percentile_exrates(one_stage_calc, pc=  '0.5', slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates_median   = get_percentile_exrates(one_stage_calc, pc=  '0.5', slip_type = slip_type, mu_type = 'variable_mu')
    output_rates_16pc                 = get_percentile_exrates(one_stage_calc, pc= '0.16', slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates_16pc     = get_percentile_exrates(one_stage_calc, pc= '0.16', slip_type = slip_type, mu_type = 'variable_mu')
    output_rates_lower_ci             = get_percentile_exrates(one_stage_calc, pc='0.025', slip_type = slip_type, mu_type = 'fixed_mu')
    variable_mu_output_rates_lower_ci = get_percentile_exrates(one_stage_calc, pc='0.025', slip_type = slip_type, mu_type = 'variable_mu')

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
    source_name, slip_type, gauge_points){

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
        'revised1_tsunami_stage_exceedance_rates_', sourcename_dot_nc)

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

for(i in 1:length(source_names)){

    # source information
    source_name = source_names[i]
    uniform_slip_tsunami_file = all_source_uniform_slip_tsunami[i]
    stochastic_slip_tsunami_file = all_source_stochastic_slip_tsunami[i]
    variable_uniform_slip_tsunami_file = all_source_variable_uniform_slip_tsunami[i]

    # Get uniform slip outputs
    uniform_slip_rates = source_zone_stage_exceedance_rates(
        source_name, slip_type='uniform', gauge_points=gauge_points)

    # Get stochastic slip outputs
    stochastic_slip_rates = source_zone_stage_exceedance_rates(
        source_name, slip_type='stochastic', gauge_points=gauge_points)

    # Get variable_uniform slip outputs
    variable_uniform_slip_rates = source_zone_stage_exceedance_rates(
        source_name, slip_type='variable_uniform', gauge_points=gauge_points)

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
