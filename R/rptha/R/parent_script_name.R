#' Find the file used to execute this comment
#'
#' If the call to \code{parent_script_name} is placed inside a script which is
#' source'd from R or executed with Rscript, then it should return the full path
#' to the corresponding file. See discussion on http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
#'
#' @return full path to the script from which parent_script_name was called
#' @export
parent_script_name <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        file_prefix <- "--file="
        file_prefix_index <- grep(file_prefix, cmdArgs)
        if (length(file_prefix_index) > 0) {
                # Rscript
                return(normalizePath(sub(file_prefix, "", cmdArgs[file_prefix_index])))
        } else {
                # 'source'd via R console
                output = try(normalizePath(sys.frames()[[1]]$ofile), silent=TRUE)
                if(class(output) != 'try-error'){
                    return(normalizePath(sys.frames()[[1]]$ofile))
                }else{
                    return('Error: Could not determine path. Possibly parent_script_name() was called from the console?')
                }
        }
}
