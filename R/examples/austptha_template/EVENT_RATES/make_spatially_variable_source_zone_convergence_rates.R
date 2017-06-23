#
# Find the 'Bird' convergence rates along our unit-source top edges
#
suppressPackageStartupMessages(library(rptha))

#
# Parse bird's data 
# 

bird_data = '../DATA/BIRD_PLATE_BOUNDARIES/PB2002_steps.dat.txt'
bd = read.table(bird_data)

# See table 2 of the bird paper for definitions of the data columns
names(bd) = c('id', 'plateboundary', 'lon1', 'lat1', 'lon2', 'lat2', 
    'length', 'azi', 'vel_L2R', 'vel_azi', 'vel_div', 'vel_rl', 
    'elev', 'age', 'class')
bird_centroid = midPoint(as.matrix(bd[,3:4]), as.matrix(bd[,5:6]), f = 0)


#
# Parse unit-source top edges
#
unit_source_files = Sys.glob(
    '../SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source_statistics*.nc')
unit_source_tables = lapply(as.list(unit_source_files), read_table_from_netcdf)
names(unit_source_tables) = unit_source_files

#
# Make 'top-edge-only' tables
#
top_edge_tables = unit_source_tables
for(i in 1:length(unit_source_tables)){
    dd = top_edge_tables[[i]]$downdip_number
    top_edge_tables[[i]] = top_edge_tables[[i]][which(dd==1),]
}

#
# Find Bird centroid nearest to top_edges
#
nearest_bird_point<-function(p){
    p_mat = bird_centroid*0
    p_mat[,1] = p[1]
    p_mat[,2] = p[2]

    distances = distHaversine(p_mat, bird_centroid)
    k = which.min(distances)
    output = c(k, distances[k])
    return(output)
}

for(i in 1:length(top_edge_tables)){
    ti = top_edge_tables[[i]]
    di = ti[,1]*0 # Store distances to nearest bird point
    ki = ti[,1]*0 # Store index of nearest bird point
    for(j in 1:nrow(ti)){
        output = nearest_bird_point(ti[j,1:2]) 
        di[j] = output[2]
        ki[j] = output[1]
    }
    top_edge_tables[[i]] = cbind(ti, 
        data.frame('distance_bird' = di, 'bird_index' = ki) )
}

