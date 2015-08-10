rptha is the main workhorse R package. To build it, go inside 'rptha', start R, and do

    source('build_package.R')

This will make an R package file in the directory above the package, which can be installed on the command line with:

    sudo R CMD INSTALL rptha_XXXXX.tar.gz

where the XXXX are adapted to match the file name.
