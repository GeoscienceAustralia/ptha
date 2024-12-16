#
# Extract an elevation contour
#

library(stars)
library(terra)

desired_elevation_contour = 0.


# High resolution elevation rasters on the mainland and Lord Howe Island
all_highres_elevation_rasts = paste0(
    './OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/elevation0_domain_', 
    c(seq(2, 488)), '.tif')
stopifnot(all(file.exists(all_highres_elevation_rasts)))

# Get points on the desired elevation contour
get_contour<-function(raster_file){
    rast = read_stars(raster_file)
    contour = st_contour(rast, breaks=desired_elevation_contour, contour_lines=TRUE)
    names(contour) = c('value', 'geometry')
    return(contour)
}

all_contours = lapply(all_highres_elevation_rasts, get_contour)
merged_contours = do.call(rbind, all_contours)

output_dir = paste0('elevation_contour_level_', desired_elevation_contour)

st_write(merged_contours, dsn=output_dir, layer=paste0(output_dir, '.shp'), driver='ESRI Shapefile', append=FALSE)
