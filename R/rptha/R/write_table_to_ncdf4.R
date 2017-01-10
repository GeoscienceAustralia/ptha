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
#' @return nothing, but save the file
#'
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
    units=NULL, long_names=NULL, var_prec = NULL){

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
                    i,' based on first entry ', datafame[i,1]))
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
    output_nc_file = nc_create(filename, vars=var_list)

    # Add global attributes
    if(!is.null(global_attributes_list)){
        for(i in 1:length(global_attributes_list)){
            ncatt_put(output_nc_file, varid=0, attname=names(global_attributes_list)[i], 
                attval=global_attributes_list[[i]])
        }
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
#' @return data.frame with the data
#' @export
#'
read_table_from_netcdf<-function(filename){

    fid = nc_open(filename)
    nc_var_names = unlist(lapply(fid$var, f<-function(x) x$name))

    var_list = list()
    for(i in 1:length(nc_var_names)){
        var_list[[nc_var_names[i]]] = c(ncvar_get(fid, varid=nc_var_names[i]))
    }

    output_df = as.data.frame(var_list, stringsAsFactors=FALSE)

    return(output_df)
}

