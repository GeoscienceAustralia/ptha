## Batch download all the historical DART data.

output_dir = 'dart_historical'
data_site = 'http://www.ndbc.noaa.gov/data/historical/dart/'

# Place for saved files
dir.create(output_dir, showWarnings=FALSE)
setwd(output_dir)

# First we get a file 'index.html' from the page, which has a list of the filenames
download_command = paste0('wget ', data_site)
system(download_command)

# Extract the filenames
dart_list = readLines('index.html')
dart_list = dart_list[grep('.txt.gz', dart_list)]
dart_names = lapply(as.list(dart_list), 
    f<-function(x) strsplit(strsplit(x, split='>', fixed=TRUE)[[1]][3], '<', fixed=TRUE)[[1]][1])
dart_names = unlist(dart_names)

# Download
for(nm in dart_names){
    download.file(paste0(data_site, nm), destfile=nm)
}

# I also separately downloaded 'dartmeta_public.xls' from their website, which gives station
# metadata including locations. (Beware that this is not really an xls file)!
