#
# Combine subuction zone trace information from multiple sources to a single table
#

sources  = c('bird_traces_table.csv', 'Arakan_traces_table.csv', 'Seramsouth_traces_table.csv', 'traces_table.csv')

#
# Data to help tweak outer-rise parameters from traces_table.csv
#

#
# Read each of the sources csv files into a big array
#
for(i in 1:length(sources)){
    newtab = read.csv(sources[i])

    if(i == 1){
        bigtab = newtab
    }else{
        bigtab = rbind(bigtab, newtab)
    }
}

# Flip Azimuth in some regions. This is done because Bird (2003) does not always put the arrows
# in the up-dip direction. We want consistency for plotting, so make this change.
library(rptha) # Provides readOGR, writeOGR
library(geosphere)
flip_regions = readOGR('Bird_reverse_vector_regions', layer='Bird_reverse_vector_regions')
pointloc = cbind(0.5*(bigtab$Long1 + bigtab$Long2), 0.5*(bigtab$Lat1 + bigtab$Lat2))
pointloc = SpatialPoints(coords=pointloc, proj4string=CRS(proj4string(flip_regions)))

to_flip = which(!is.na(over(pointloc, flip_regions)))

bigtab_test = bigtab
bigtab_test$Azi_Vel[to_flip] = (bigtab_test$Azi_Vel[to_flip] - 180)%%360

write.csv(bigtab_test, 'sourcezone_traces_table_merged_uncorrected.csv', row.names=FALSE)
