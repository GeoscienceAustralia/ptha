#
# Shortcut to getting elevation filesnames from QGIS, which I have already ordered well to make a nice DEM.
#

# I saved a QGIS file that included the layers in desired order (and other things) to the old format, which is scrapable, then renamed the file.
#x = readLines('FILE_TO_SCRAPE_ELEVATION_RASTERS_FROM.txt')
x = readLines('Revised_file_to_scrape_elevation_files_from.qgs')

# Find text corresponding to files in the 'BATHY' group
# From experience, the information I want is between these two indices
b1 = grep('BATHY', x)[1]
b2 = grep('reference_bathy', x)[1]

# Subset the relevant parts
y = x[b1:b2]
z = y[grep('source=', y)] # Something with a filename

# Extract filenames
filenames = sapply(z, function(x){
    x_split = strsplit(x, split="=")[[1]]
    k = which(endsWith(x_split, "source"))
    if(length(k) != 1) stop('error')
    k = k+1
    output = gsub(' legend_exp', '', x_split[k], fixed=TRUE)
    return(output)
    }, USE.NAMES=FALSE)
# Clean up quotes
filenames = gsub('\"', '', filenames, fixed=TRUE)

# The first entry is a shapefile that I don't want
if(grepl('OEH_singlebeam', filenames[1])) filenames = filenames[-1]

filenames_exist = file.exists(filenames)
if(!all(filenames_exist)){
    print(cbind(filenames, filenames_exist))
    stop('Some files do not exist! Halting')
}

# Easier to work with full paths
normal_filenames = normalizePath(filenames)

# Get the md5sums
library(tools)
files_md5 = md5sum(normal_filenames)
names(files_md5) = ""

output = data.frame(files = normal_filenames, files_md5 = files_md5)

write.csv(output, file='Elevation_rasters_scraped_from_QGIS_file.csv', row.names=FALSE)
