#library(ncdf4)

#' Write a data.frame to a netcdf file
#'
#' The unlimited dimension corresponds to the rows of the data.frame, which 
#' should enable efficient row-wise access to the data.
#'
#' @param dataframe the input data.frame to be written to file
#' @param filename the output filename
#' @param global_attributes_list a named list of global attributes to add to the file, or NULL
#' @param units character vector with units for each column of the dataframe, or NULL
#' @param long_names character vector with long_name for each column of the dataframe, or NULL
#' @param var_prec character vector with data type for each column (possible
#' values are 'short', 'integer', 'float', 'double', 'char', 'byte'). If NULL it is made automatically.
#' @param add_session_info_attribute if true, add a global attribute containing information from sessionInfo()
#' @param force_v4 logical Do we create a netcdf4 file (TRUE) or not (FALSE, default)
#' @return nothing, but save the file
#' @import ncdf4
#' @export
#'
#' @examples
#'
#'  test_table = data.frame(x = c(1,2,3), y=c('a', 'b', 'csdf'), z=c(1.1, 1.2, 1.3))
#'
#'  write_table_to_netcdf(test_table, 
#'      file='test.nc', 
#'      units=c('m', '', 'kg'), 
#'      long_names=c('head count', 'mychar', 'asdfasdfa'), 
#'      var_prec=c('float', 'char', 'float'))
#'
#'  # Clean up
#'  unlink('test.nc')
#'
write_table_to_netcdf<-function(dataframe, filename, global_attributes_list=NULL, 
    units=NULL, long_names=NULL, var_prec = NULL, add_session_info_attribute=FALSE,
    force_v4=FALSE){

    # Make the rows an unlimited dimension
    rowdim = ncdim_def('table_rows', units='', vals = 1:nrow(dataframe), unlim=TRUE,
        longname='dimension for rows of table')

    # Units -- empty by default
    if(is.null(units)){
        units = rep("", ncol(dataframe))
    }

    # long_names -- empty by default
    if(is.null(long_names)){
        long_names = rep("", ncol(dataframe))
    }

    # Var-prec -- by default, try to classify as 'integer', 'double', or 'character'
    if(is.null(var_prec)){
        var_prec = rep("", ncol(dataframe))
        for(i in 1:ncol(dataframe)){
            if(class(dataframe[1,i]) %in% c('numeric', 'integer')){
                # Check if they can be represented as integers
                if(all(floor(dataframe[,i]) == dataframe[,i])){
                    var_prec[i] = 'integer' 
                }else{
                    var_prec[i] = 'double'
                }
            }else if(class(dataframe[1,i]) %in% c('character', 'factor')){
                var_prec[i] = 'char'
            }else{
                stop(paste0('Cannot automatically detect var_prec for column ', 
                    i,' based on first entry ', dataframe[i,1]))
            }
        }
    }

    # Check that var_prec is reasonable
    allowed_var_precs = c('short', 'integer', 'float', 'double', 'char', 'byte')
    if(!(all(var_prec %in% allowed_var_precs))){
        stop(paste0('illegal var_prec value ', setdiff(var_prec, allowed_var_precs)))
    }

    # Check lengths of inputs
    stopifnot(length(long_names) == ncol(dataframe))
    stopifnot(length(var_prec) == ncol(dataframe))
    stopifnot(length(units) == ncol(dataframe))


    # Make dimension for any character variables
    if(any(var_prec == 'char')){
        chardim = 0
        for(i in 1:ncol(dataframe)){
            if(var_prec[i] == 'char'){
                chardim = max(chardim, max(nchar(as.character(dataframe[,i]))))
            }
        }
        chardim = ncdim_def('max_nchar', units='', vals=1:chardim, 
            longname='dimension giving maximum number of string characters')
    }

    # Make variables
    var_list = list()
    for(i in 1:ncol(dataframe)){
        name_i = names(dataframe)[i]
        if(var_prec[i] == 'char'){
            var_list[[name_i]] = ncvar_def(name_i, units=units[i], 
                dim=list(chardim, rowdim), 
                longname=long_names[i], prec=var_prec[i])
            # Ensure factors are also characters
            dataframe[,i] = as.character(dataframe[,i])
        }else{
            var_list[[name_i]] = ncvar_def(name_i, units=units[i], 
                dim=rowdim, 
                longname=long_names[i], prec=var_prec[i])
        }
    }

    # Make file
    output_nc_file = nc_create(filename, vars=var_list, force_v4 = force_v4)

    # Add global attributes
    if(!is.null(global_attributes_list)){
        for(i in 1:length(global_attributes_list)){
            ncatt_put(output_nc_file, varid=0, attname=names(global_attributes_list)[i], 
                attval=global_attributes_list[[i]])
        }
    }

    if(add_session_info_attribute){
        ncatt_put(output_nc_file, varid=0, attname='R_session_info', 
            attval=paste(capture.output(sessionInfo()), collapse=" ; "))
    }

    # Add variables
    for(i in 1:ncol(dataframe)){
        ncvar_put(output_nc_file, var_list[[i]], dataframe[,i]) 
    }

    nc_close(output_nc_file)

}

#' Read a netcdf file into a data.frame
#'
#' It is assumed that the file was produced by \code{write_table_to_netcdf},
#' or is compatible with it.
#'
#' @param filename netcdf file
#' @param desired_rows integer vector giving the rows to extract. If null, read everything.
#'   Note that rows are read in contiguous chunks (up to 100 at once) for efficiency. For
#'   small datasets it may well be faster to read everything. However, for large datasets
#'   over remote connections this is not possible.
#' @return data.frame with the data
#' @import ncdf4
#' @export
#'
read_table_from_netcdf<-function(filename, desired_rows = NULL){

    fid = nc_open(filename)
    nc_var_names = unlist(lapply(fid$var, f<-function(x) x$name))

    var_list = list()
    for(i in 1:length(nc_var_names)){
        #print(i)
        if(is.null(desired_rows)){
            # Get the entire variable
            var_list[[nc_var_names[i]]] = c(ncvar_get(fid, varid=nc_var_names[i]))
        }else{
            #
            # Get only a few rows
            # Read them in chunks of size chunk_size (defined below),
            # because that tends to be faster. 
            #

            temp_var = rep(NA, length(desired_rows))
            n = fid$var[[nc_var_names[i]]]$ndim

            max_rows = fid$var[[nc_var_names[i]]]$varsize[n]
            chunk_size = 100

            for(j in seq(min(desired_rows), max(desired_rows), by=chunk_size)){

                search_inds = j:(min(max_rows, j+chunk_size-1))

                k = which(desired_rows %in% search_inds)
                if(length(k) == 0) next
                index_match = match(desired_rows[k], search_inds)
                search_inds = search_inds[1:max(index_match)]

                # Characters will have > 1 dimension, so we need this trick
                # to specify start/end
                start = rep(1, length=n)
                count = rep(-1, length=n)

                start[n] = j
                count[n] = diff(range(search_inds)) + 1

                myvar = c(ncvar_get(fid, varid=nc_var_names[i], start=start, count=count))

                temp_var[k] = myvar[index_match]
            }
            var_list[[nc_var_names[i]]] = temp_var
        }
    }

    output_df = as.data.frame(var_list, stringsAsFactors=FALSE)

    return(output_df)
}

