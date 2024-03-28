#
# Shortcut to getting elevation filesnames from QGIS, which I have already ordered well to make a nice DEM.
#

# I saved a QGIS file that included the layers in desired order (and other things) to the old format, which is scrapable, then renamed the file.
#x = readLines('FILE_TO_SCRAPE_ELEVATION_RASTERS_FROM.txt')
x = readLines('NSW_Tsunami_Inundation_2023_2025.qgz.qgs')

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
# Remove a word at the end
for(i in 1:length(filenames)){
    # Remove whatever comes after the final space
    filenames_split = strsplit(filenames[i], split="*")[[1]]
    space_ind = which(filenames_split == " ")
    if(length(space_ind) > 0){
        nc = max(space_ind)
        filenames[i] = substring(filenames[i], 1, nc-1)
    }
    ## remove " id"
    #if(endsWith(filenames[i], " id")) filenames[i] = substring(filenames[i], 1, nchar(filenames[i]) - 3)
    ## remove " checked
    #if(endsWith(filenames[i], " checked")) filenames[i] = substring(filenames[i], 1, nchar(filenames[i]) - 8)
    ##
    #if(endsWith(filenames[i], " patch_size")) filenames[i] = substring(filenames[i], 1, nchar(filenames[i]) - 11)
}

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
