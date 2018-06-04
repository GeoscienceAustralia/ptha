library(rptha)

#' Append variable-shear-modulus variables to earthquake events file
#'
#' This routine takes the earthquake events file and adds variables
#' related to 'variable shear modulus' event properties. (The file is 
#' initially created with information on Mw & rates, but not variable_mu
#' versions of those).
#'
#' @param netcdf_file name of netcdf file containing earthquake events (but not
#' tsunami data)
#' @return nothing, but append variables to the file 
#' 
append_variable_mu_variables_to_earthquake_event_files<-function(netcdf_file){

    # 
    fid = nc_open(netcdf_file, readunlim=FALSE, write=TRUE)
 
    # Check we've got the right kind of file
    if(!all(sort(names(fid$dim)) == c('max_nchar', 'table_rows'))){
        nc_close(fid)
        stop(paste0(netcdf_file, ' is not the correct type of file for this function (wrong dimensions)'))
    }

    vars_to_expect = c('Mw', 'rate_annual', 'rate_annual_lower_ci', 
        'rate_annual_upper_ci', 'event_index_string', 'sourcename')

    if(!all( vars_to_expect %in% names(fid$var) )){
        nc_close(fid)
        stop(paste0(netcdf_file, ' is not the correct type of file for this function (missing expected variable names)'))
    }

    vars_to_add = paste0('variable_mu_', c('Mw', 'rate_annual', 
        'rate_annual_lower_ci', 'rate_annual_upper_ci'))

    # Check whether we have already appended the variables
    if(any(vars_to_add %in% names(fid$var))){
        if(!all(vars_to_add %in% names(fid$var))){
            nc_close(fid)
            stop(paste0(netcdf_file, ' has some variable_mu_ variables appended, but not all'))
        }
        # Otherwise, the file has already been created with the right stuff,
        # and we don't need to do anything
    }else{
        # Add the variables
        table_rows_dim = fid$dim$table_rows
        for(varname in vars_to_add){
            new_var = ncvar_def(varname, units="", dim=table_rows_dim, prec='double') 
            fid = ncvar_add(fid, new_var)
            ncvar_put(fid, new_var, vals=rep(-1.0, table_rows_dim$len))
        }
    }

    # Flush to file   
    nc_close(fid) 
}

compute_variable_mu_magnitude_and_update_file<-function(netcdf_file, unit_source_statistics_file, 
    update_file=TRUE, assume_fixed_mu=FALSE){

    fid = nc_open(netcdf_file, readunlim=FALSE, write=TRUE)
 
    # Check we've got the right kind of file
    if(!all(sort(names(fid$dim)) == c('max_nchar', 'table_rows'))){
        nc_close(fid)
        stop(paste0(netcdf_file, ' is not the correct type of file for this function (wrong dimensions)'))
    }

    vars_to_expect = c('Mw', 'rate_annual', 'rate_annual_lower_ci', 
        'rate_annual_upper_ci', 'event_index_string', 'sourcename',
        'variable_mu_Mw', 'variable_mu_rate_annual', 
        'variable_mu_rate_annual_lower_ci',
        'variable_mu_rate_annual_upper_ci')

    if(!all( vars_to_expect %in% names(fid$var) )){
        nc_close(fid)
        stop(paste0(netcdf_file, ' is not the correct type of file for this function (missing expected variable names)'))
    }

    variable_mu_Mw = ncvar_get(fid, 'variable_mu_Mw')
    if(all(variable_mu_Mw > 0)){
        true_mu_Mw = ncvar_get(fid, 'Mw')
        if(any(abs(true_mu_Mw - variable_mu_Mw) > 0.4)){
            stop('Problem with variable mu')
        }else{
            # We don't need to do anything -- the file was already updated
        }
    }else{

        if(assume_fixed_mu){
            # Just set variable-mu-mw equal to fixed-mu-mw
            true_mu_Mw = ncvar_get(fid, 'Mw')
            variable_mu_Mw = true_mu_Mw
        }else{
            # Compute magnitude
            eis = ncvar_get(fid, 'event_index_string')
            eis = sapply(eis, f<-function(x) as.numeric(strsplit(x, '-')[[1]]), simplify=FALSE)
            uss = read_table_from_netcdf(unit_source_statistics_file)
            uss_area = uss$length * uss$width
            uss_depth = uss$depth
            uss_mu_variable = shear_modulus_depth(uss_depth)
           
            variable_mu_Mw = rep(NA, length(eis)) 

            if('event_slip_string' %in% names(fid$var)){
                # Stochastic or variable uniform slip
                ess = ncvar_get(fid, 'event_slip_string')
                ess = sapply(ess, f<-function(x) as.numeric(strsplit(x, '_')[[1]]), simplify=FALSE)
            }else{
                # Uniform slip without size variability
                ess = as.list(ncvar_get(fid, 'slip'))
            }

            for(i in 1:length(variable_mu_Mw)){
                inds = eis[[i]]
                slip = ess[[i]]
                moment = sum(uss_area[inds] * 1e+06 * slip * uss_mu_variable[inds])
                variable_mu_Mw[i] = M0_2_Mw(moment)
            }
        }

        if(update_file){
            ncvar_put(fid, 'variable_mu_Mw', variable_mu_Mw)

        }
    }
    nc_close(fid)

    return(invisible(variable_mu_Mw))
}


config_env = new.env()
source('config.R', local=config_env)

for(i in 1:length(config_env$unit_source_statistics_netcdf_files)){

    # Get the unit-source statistics + event files (without tsunami) 
    # for each slip generation method
    uss = config_env$unit_source_statistics_netcdf_files[i]
    stoc_file = gsub('_tsunami', '', config_env$all_source_stochastic_slip_tsunami[i], fixed=TRUE)
    unif_file = gsub('_tsunami', '', config_env$all_source_uniform_slip_tsunami[i], fixed=TRUE)
    vary_unif_file = gsub('_tsunami', '', config_env$all_source_variable_uniform_slip_tsunami[i], fixed=TRUE)

    # For the normal faults, we currently do not implement a varible shear modulus treatment
    # In our schematization they often only have 1 unit source!
    uss_table = read_table_from_netcdf(uss)
    if(isTRUE(all.equal(uss_table$rake[1], -90))){
        fixed_mu_only = TRUE
    }else if(isTRUE(all.equal(uss_table$rake[1], 90))){
        fixed_mu_only = FALSE
    }else{
        print(c(uss, range(uss_table$rake)))
        stop('The code allows rake = -90 or 90 ONLY, but this file violates that constraint')
    }

    # Update stochastic slip
    print(paste0('Updating ', stoc_file))
    append_variable_mu_variables_to_earthquake_event_files(stoc_file)
    compute_variable_mu_magnitude_and_update_file(stoc_file, uss, assume_fixed_mu=fixed_mu_only)

    # Update uniform slip
    print(paste0('Updating ', unif_file))
    append_variable_mu_variables_to_earthquake_event_files(unif_file)
    compute_variable_mu_magnitude_and_update_file(unif_file, uss, assume_fixed_mu=fixed_mu_only)
    
    # Update variable-uniform slip
    print(paste0('Updating ', vary_unif_file))
    append_variable_mu_variables_to_earthquake_event_files(vary_unif_file)
    compute_variable_mu_magnitude_and_update_file(vary_unif_file, uss, assume_fixed_mu=fixed_mu_only)

}


