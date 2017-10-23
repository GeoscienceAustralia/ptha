#
# For the PTHA, we want netcdf files in */TSUNAMI_EVENTS/*.nc to
# be on gdata (so they can be remotely accessed)
#
# However, the code is structured assuming they are in local files
#
# A solution is to move the files to gdata, and replace the local files
# with symbolic links pointing to the right file on gdata
#

new_basedir = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/'

nc_files = Sys.glob('*/TSUNAMI_EVENTS/*.nc')

for(i in 1:length(nc_files)){

    nc_file = nc_files[i]
    print(nc_file)
    new_nc_name = paste0(new_basedir, nc_file) 
    dir.create(dirname(new_nc_name), recursive=TRUE, showWarnings=FALSE)

    # Copy to gdata
    file.copy(nc_file, new_nc_name)
    # Rename local version (can delete later)
    file.rename(nc_file, paste0(nc_file, 'BACKUP'))
    # Make symlink
    file.symlink(new_nc_name, nc_file)
    
}
