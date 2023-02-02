#
# Code to reorder the merged hazard points, such that nearby points are close
# to each other.
#
# It is not essential to use this if the points are already 'well ordered'.
# But it can be used to make nearby points closer to each other in the file,
# which should improve read performance [the effect could be dramatic if the
# points were not initially 'well ordered'].
#
# Note we do not include the DART points in the re-ordering. They are kept
# separate at the end of the hazard_points. The reason is that typically,
# users can separate their needs into extraction of either A) a cluster of
# nearshore points, or B) a set of dart points. We want to be able to efficiently
# extract both these cases [in particular, to easily get all dart points], so
# it's best to keep the dart points separate, not merge them into the remaining
# hazard points.
#
# If the original points are generated along contour lines, then the ordering
# might be quite ok -- e.g. a given bounding box will typically only contain a few
# clusters of points. Other methods (e.g. leading to latitude-sorting of points) could
# be far worse
#
# One approach is to reorder with a travelling salesman type
# ordering. However, this does not necessarily ensure that nearby points are close
# to each other in the file [or at least in a few distinct chunks].
# 
# We can also order into N-degree x N-degree boxes. This might be convenient and
# quite predictable
#
#
# Also, in any case, we need to implement an algorithm in the data getting routines,
# which can split indices of points into 'nearly contiguous chunks' and read that way.
#
library(rptha) # provides readOGR, writeOGR

hp = read.csv('merged_hazard_points.csv')

#hp = hp[1:2000,]

# Get indices of dart and non-dart points
# Dart points have a label ending in '.4'
dart_points = which( abs(abs(hp[,3] - floor(hp[,3])) - 0.4) < 1.0e-04)
non_dart = setdiff(1:length(hp[,1]), dart_points)

#
# Re-ordering code here 
#

use_travelling_salesman = FALSE

if(use_travelling_salesman){

    library(TSP)
    library(geosphere)

    # Compute spherical coordinates distance matrix
    big_distm = distm(hp[non_dart,1:2])

    # Convert to TSP object
    hp_tsp = TSP(big_distm, labels=hp[non_dart,3])

    # Get an approximate solution. Run time estimate about 3.5 hours for all
    # points, on my machine [single core only]
    hp_tour = solve_TSP(hp_tsp, method = 'cheapest_insertion')

    # Reorder the points
    hp_reordered_nd = hp[non_dart[hp_tour],]

    if(length(dart_points) > 0){
        # Add the dart points on the end
        hp_reordered = rbind(hp_reordered_nd, hp[dart_points,])
    }
 
    write.csv(hp_reordered, 'merged_hazard_points_ts.csv', row.names=FALSE)

}else{
    # A natural alternative would be to group in e.g. 2deg x 2deg blocks

    hp_reordered = rbind(hp[non_dart,], hp[dart_points,])
    n = length(non_dart)

    # Put an order on the 'non-dart' points
    ndeg = 2
    hpA = floor(hp_reordered[1:n,1]/ndeg)
    hpB = floor(hp_reordered[1:n,2]/ndeg)

    ord = order(hpA, hpB)

    # Reorder the non-dart
    hp_reordered[1:n,] = hp_reordered[ord,]
    
    write.csv(hp_reordered, 'merged_hazard_points_blocked.csv', row.names=FALSE)
}



#
#
#


# For plotting, do not depict lines joining the right/left extremes of the
# map as connected, since they are 'close' and would wrap from left to right
png('hp_before_after_ordering.png', width=12, height=12, units='in', res=300)
par(mfrow=c(2,1))

wrap_cut_threshold = 200

plot(range(hp[,1]), range(hp[,2]), asp=1, col=0, main='Original order', cex.main=2, xlab='Lon', 
    ylab = 'Lat')

# ID dart bupys
ii = which(abs( (hp[,3] - trunc(hp[,3])) - 0.4) < 1.0e-05)
n = length(hp[-ii,1])
line_splits_hp = c(1, which(abs(diff(hp[-ii,1])) > wrap_cut_threshold) + 1, n+1)

LinesList_orig = list()
for(i in 1:(length(line_splits_hp)-1)){
    start = line_splits_hp[i]
    end = line_splits_hp[i+1] - 1
    points(hp[start:end,1:2], t='l')

    LinesList_orig[[i]] = Lines(list(Line(hp[start:end,1:2])), ID=as.character(i))
}
# Add dart buoys [which are at the end]
points(hp[ii,1:2], col='blue',pch=19, cex=0.5)

plot(range(hp_reordered[,1]), range(hp_reordered[,2]), asp=1, col=0, cex.main=2,
    main='Updated order [separate DART buoys]', xlab='Lon', ylab='Lat')

# ID dart bupys
ii = which(abs( (hp_reordered[,3] - trunc(hp_reordered[,3])) - 0.4) < 1.0e-05)

n = length(hp_reordered[-ii,1])
line_splits_hp = c(1, which(abs(diff(hp_reordered[-ii,1])) > wrap_cut_threshold) + 1, n+1)

LinesList_new = list()
for(i in 1:(length(line_splits_hp)-1)){
    start = line_splits_hp[i]
    end = line_splits_hp[i+1] - 1
    points(hp_reordered[start:end,1:2], t='l', col='green')
    LinesList_new[[i]] = Lines(list(Line(hp_reordered[start:end,1:2])), ID=as.character(i))
}
# Add dart buoys [which are at the end]
points(hp_reordered[ii,1:2], col='blue',pch=19, cex=0.5)

dev.off()


# Save as shapefiles
hp_lines1 = SpatialLines(LinesList_orig, proj4string=CRS("+init=epsg:4326"))
hp_lines1b = SpatialLinesDataFrame(hp_lines1, data=data.frame(ID=1:length(hp_lines1)))
writeOGR(hp_lines1b, dsn='hp_before_ordering', layer='hp_before_ordering', 
    driver='ESRI Shapefile', overwrite=TRUE)

hp_lines2 = SpatialLines(LinesList_new, proj4string=CRS("+init=epsg:4326"))
hp_lines2b = SpatialLinesDataFrame(hp_lines2, data=data.frame(ID=1:length(hp_lines1)))
writeOGR(hp_lines2b, dsn='hp_after_ordering', layer='hp_after_ordering', 
    driver='ESRI Shapefile', overwrite=TRUE)

#
# Graphical example of localisation
#
if(FALSE){

    # Find which points inside bounding box -- example near Sydney
    bbox_lower_left = c(149, -37)
    bbox_upper_right = c(155, -30)
    hp_inside = which(
        hp[,1] > bbox_lower_left[1] & 
        hp[,2] > bbox_lower_left[2] & 
        hp[,1] < bbox_upper_right[1] & 
        hp[,2] < bbox_upper_right[2])
    range(diff(hp_inside))
    
    hpnew_inside = which(
        hp_reordered[,1] > bbox_lower_left[1] & 
        hp_reordered[,2] > bbox_lower_left[2] & 
        hp_reordered[,1] < bbox_upper_right[1] & 
        hp_reordered[,2] < bbox_upper_right[2])
    range(diff(hpnew_inside))
    
    par(mfrow=c(1,2))
    plot(hp_inside)
    plot(hpnew_inside)
    
    # Find which points inside bounding box -- example near Mentawai
    bbox_lower_left = c(93, -7)
    bbox_upper_right = c(104, -3)
    hp_inside = which(
        hp[,1] > bbox_lower_left[1] & 
        hp[,2] > bbox_lower_left[2] & 
        hp[,1] < bbox_upper_right[1] & 
        hp[,2] < bbox_upper_right[2])
    range(diff(hp_inside))
    
    hpnew_inside = which(
        hp_reordered[,1] > bbox_lower_left[1] & 
        hp_reordered[,2] > bbox_lower_left[2] & 
        hp_reordered[,1] < bbox_upper_right[1] & 
        hp_reordered[,2] < bbox_upper_right[2])
    range(diff(hpnew_inside))
    
    par(mfrow=c(1,2))
    plot(hp_inside)
    plot(hpnew_inside)
    
    
    # Find which points inside bounding box -- example near Perth
    bbox_lower_left = c(110, -35)
    bbox_upper_right = c(118, -28)
    hp_inside = which(
        hp[,1] > bbox_lower_left[1] & 
        hp[,2] > bbox_lower_left[2] & 
        hp[,1] < bbox_upper_right[1] & 
        hp[,2] < bbox_upper_right[2])
    range(diff(hp_inside))
    
    hpnew_inside = which(
        hp_reordered[,1] > bbox_lower_left[1] & 
        hp_reordered[,2] > bbox_lower_left[2] & 
        hp_reordered[,1] < bbox_upper_right[1] & 
        hp_reordered[,2] < bbox_upper_right[2])
    range(diff(hpnew_inside))
    
    par(mfrow=c(1,2))
    plot(hp_inside)
    plot(hpnew_inside)
}
