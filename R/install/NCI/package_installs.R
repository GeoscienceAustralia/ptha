#
# How to setup R installation on raijin, NCI, to use the rptha package
#
# Be sure to read the comments in this script, they ask you to do essential
# things outside of R
#
#

# Before starting R (everytime), you need to: 
#     source R_modules.sh 
# on the command line, to load the required modules.
# The 'R_modules.sh' file is in the current directory.  
#
# I suggest putting a copy in your home directory, so you can do
#      source ~/R_modules.sh 
# wherever you are.
#
# You should also copy the hidden '.R' folder in the current directory to your
# home directory. In contains default configuration information for R which
# will help with package compilation. You only need to do this once.
#


# Now we start installing packages

# Start R and run this command. The packages should install cleanly if
# R_modules.sh has been sourced beforehand
install.packages(c('sp', 'geosphere', 'FNN', 'minpack.lm', 'geometry', 'raster', 'testthat', 'ncdf4'))

# Now we will install rgdal, which is R's gdal interface. It is problematic because
# it depends on various other software which is installed in slightly non-standard places
# on NCI.
#
# Here we need to explicitly link to proj_api.h and libproj.so
# I was not able to get it to work with the versions of proj available
# on NCI, so built proj.4.9.1 separately here
proj4_build_dir = '/short/w85/gxd547/PTHA_Aust/SOURCE/proj4/proj.4.9.1/proj-4.9.1/build'
#
# Make sure you do not have another version of proj loaded. 
# Use 'module list'
# to see the modules you have, and 'module unload XXXX' to unload XXXX
install.packages('rgdal',
    configure.args=c(paste0('--with-proj-include=', proj4_build_dir, '/include'),
                     paste0('--with-proj-lib=', proj4_build_dir, '/lib')))


# rgeos should build nicely
install.packages('rgeos')


# devtools needs a few tweaks to build
zlib_build_dir = '/apps/zlib/1.2.8'
install.packages('devtools',
    configure.args=c(paste0('--with-zlib-include=', zlib_build_dir, '/include'),
                     paste0('--with-zlib-lib=', zlib_build_dir, '/lib'),
                     '--disable-cxx11'))

# roxygen2 also needs some tweaks, related to the stringi package
install.packages('stringi', configure.args='--disable-cxx11')
install.packages('roxygen2')


# Then you should be able to build rptha and install it just as described in the package README.md
