#
# Extract max-stage at a desired elevation contour, as an indication of nearshore wave heights
#

library(stars)
library(terra)

desired_elevation_contour = 1.2
initial_stage = 1.1
desired_elevation_contour_name = '+1.2m'
desired_elevation_contour_filename_stub = 'pos1.2m'


# High resolution elevation rasters on the mainland and Lord Howe Island
all_highres_elevation_rasts = paste0(
    './OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/elevation0_domain_', 
    c(seq(261, 457), seq(476, 488)), '.tif')
stopifnot(all(file.exists(all_highres_elevation_rasts)))

# Get points on the desired elevation contour
get_contour<-function(raster_file){
    rast = read_stars(raster_file)
    contour = st_contour(rast, breaks=desired_elevation_contour, contour_lines=TRUE)
    coords = st_coordinates(contour)
    return(coords)
}

all_contour_points = lapply(all_highres_elevation_rasts, get_contour)
# Combine, ignoring 'empty' rasters (which do not contain the desired elevation)
keep = unlist(lapply(all_contour_points, function(x) nrow(x) > 0))
merged_pts = do.call(rbind, all_contour_points[keep])

# Get the max stage rasters
# Use terra, because currently 'stars' has some problems in correctly interpreting NA values
max_stage_rast_convergence = rast(Sys.glob(
    './OUTPUTS/run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1/RUN_20241112_173726773/all_max_stage_CONVERGENCE.vrt'))
max_stage_rast_regular = rast(Sys.glob(
    './OUTPUTS/run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1/RUN_20241112_155633335/all_max_stage.vrt'))

max_stage_values_convergence = extract(max_stage_rast_convergence, merged_pts[,1:2], method='bilinear')
max_stage_values_regular     = extract(max_stage_rast_regular    , merged_pts[,1:2], method='bilinear')

outputs = cbind(merged_pts, max_stage_values_convergence, max_stage_values_regular)

# NSW mainland, points where the max-stage exceeded the initial stage
k = which((outputs[,2] > -37.6) & (outputs[,2] < -28.16) & outputs[,4] > initial_stage)

# Error statistics
relerr = (outputs[k,4] - outputs[k,5])/(0.5*(outputs[k,4]+outputs[k,5]))
err = (outputs[k,4] - outputs[k,5])
print('> mean(abs(err) < 0.05, na.rm=TRUE)')
print(mean(abs(err) < 0.05, na.rm=TRUE))
print('> mean(abs(err) < 0.3, na.rm=TRUE)')
print(mean(abs(err) < 0.3, na.rm=TRUE))
print('> mean(abs(relerr) < 0.05, na.rm=TRUE)')
print(mean(abs(relerr) < 0.05, na.rm=TRUE))
print('> mean(abs(relerr) < 0.15, na.rm=TRUE)')
print(mean(abs(relerr) < 0.15, na.rm=TRUE))

## PRINT OUTS FROM 2024/11/22
#[1] "> mean(abs(err) < 0.05, na.rm=TRUE)"
#[1] 0.6954728
#[1] "> mean(abs(err) < 0.3, na.rm=TRUE)"
#[1] 0.9541697
#[1] "> mean(abs(relerr) < 0.05, na.rm=TRUE)"
#[1] 0.8369443
#[1] "> mean(abs(relerr) < 0.15, na.rm=TRUE)"
#[1] 0.9695503

png(paste0('tsunami_maxima_error_statistics_', desired_elevation_contour_filename_stub, '_elevation.png'), 
    width=12, height=5, units='in', res=200)
par(mfrow=c(1,2))
#hist(relerr, breaks=201); abline(v=seq(-1,1, by=0.1), col='orange', lty='dotted')
hist(err, breaks=201, xlim=c(-2,2), freq=FALSE,
    main=paste0('Maximum water level difference between high-resolution \n and regular resolution along ', 
                desired_elevation_contour_name, ' elevation contour'), 
    xlab='Maximum water level difference (m)',
    cex.main=1.2, cex.axis=1.5, cex.lab=1.5); abline(v=seq(-5,5, by=0.25), col='orange', lty='dotted') #grid(col='orange')
#plot(outputs[k,4], outputs[k,5], asp=1, pch='.', col=rgb(0, 0, 0, alpha=0.05, maxColorValue=1))
plot(outputs[k,4], outputs[k,5], asp=1, pch='.', col=rgb(0, 0, 0, alpha=0.2, maxColorValue=1), cex=2,
    xlim=c(0, max(outputs[k,4:5], na.rm=TRUE)), 
    xlab='High-resolution max water level (m)', ylab='Regular resolution max water level (m)',
    main='Maximum water level along +1m elevation contour', cex.main=1.3, cex.axis=1.3, cex.lab=1.3)
grid(col='orange')
abline(0, 1, col='red')
dev.off()

png(paste0('tsunami_maxima_at_', desired_elevation_contour_filename_stub, '_elevation.png'), 
    width=12, height=7, units='in', res=200)
par(mfrow=c(1,2))
par(mar=c(4,4,4,2))
plot(outputs[k,4], outputs[k,2], pch='.', cex=2, col='blue', 
     xlim=c(max(outputs[k,4], na.rm=TRUE), 0), main="", xlab="", ylab="", 
     cex.lab=1.5, cex.axis=1.5)
mtext(side=3, paste0('High-resolution max water level (m) \n along ', 
    desired_elevation_contour_name, ' elevation contour'), cex=1.5, line=1)
mtext(side=2, 'Latitude (degrees)', cex=1.5, line=2.5)
mtext(side=1, 'Maximum water level (m)', cex=1.5, line=2.5)
grid(col='orange')
plot(outputs[k,5], outputs[k,2], pch='.', cex=2, col='brown', 
     xlim=c(0, max(outputs[k,4], na.rm=TRUE)), main="", xlab="", ylab="", 
     cex.lab=1.5, cex.axis=1.5)
mtext(side=3, paste0('Standard resolution max water level (m) \n along ', 
    desired_elevation_contour_name, ' elevation contour'), cex=1.5, line=1)
mtext(side=2, 'Latitude (degrees)', cex=1.5, line=2.5)
mtext(side=1, 'Maximum water level (m)', cex=1.5, line=2.5)
grid(col='orange')
dev.off()

