#'
#' Make the gauge locations for the Gladstone model
#'
#' Include gauges at:
#' - DARTS
#' - NZ DARTS
#' - Tide gauges in eastern Australia [HTHH paper + other TG data]
#' - PTHA18 points near eastern Australia (but not all global pts, to manage file sizes)

library(terra)
library(geosphere)

find_close_locations<-function(ALL_DARTS, threshold){
    # Handy routine for combining gauges.
    # ALL_DARTS must have columns lon/lat, giving the location in degrees
    # Find a subset of sites in ALL_DARTS such that the full set of sites
    # is always within a distance "threshold" of the subset. 
    # Result is to_keep, an array of row indices from ALL_DARTS to keep
    to_keep = rep(NA, nrow(ALL_DARTS))
    to_keep[1] = TRUE
    for(i in 2:nrow(ALL_DARTS)){
        w_t_k = which(to_keep)
        p0 = cbind(ALL_DARTS$lon[w_t_k], ALL_DARTS$lat[w_t_k]) # Points we know we will keep
        p1 = p0*0 # point i
        p1[,1] = ALL_DARTS[i,1]
        p1[,2] = ALL_DARTS[i,2]
        dists = distHaversine(p0, p1)
        to_keep[i] = all(dists > threshold)
    }
    return(to_keep)
}


#
# DARTS
#
merge_darts<-function(){
    # NZ darts
    nz_darts = read.csv('nz_dart_locations/dart_locations.csv')
    nz_darts = data.frame(lon=nz_darts$lon, lat=nz_darts$lat, id=70000 + 1:nrow(nz_darts)+0.1)

    # Regular DARTS
    darts_2 = vect('dart_locations/dart_locations2.shp')
    darts_2 = data.frame(lon=darts_2$lon + 360.0*(darts_2$lon < -40), lat=darts_2$lat, id=as.numeric(darts_2$WMO_ID)+0.2)

    # Higher-res DART records. These mostly have similar or nearby locations as
    # darts2, but we see the gauges moving over time.
    darts_hr = read.csv('dart_locations/all_DART_nc_files_and_coordinates.csv')
    darts_hr = data.frame(lon=darts_hr$lon, lat=darts_hr$lat, 
        # Use a second vector to break ties in ID names
        id=as.numeric(darts_hr$dartID)+rep(c(0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39), length.out=nrow(darts_hr)))

    all_darts = rbind(nz_darts, darts_2, darts_hr)

    to_keep = find_close_locations(all_darts, threshold=5000)
    all_darts = all_darts[which(to_keep),]

    return(all_darts)
}

#
# Coastal tide gauges
#
merge_coastal_tidegauges<-function(bounding_box){
    # HTHH paper locations
    hthh_gauges = read.csv('hthh_paper_tide_gauges/01_tide_gauge_locations.csv')
    hthh_gauges = data.frame(lon=hthh_gauges$lon, lat=hthh_gauges$lat, id=100000 + 1:nrow(hthh_gauges) + 0.1)

    # Tide gauges from our previous studies with gauge data
    gauge_coords = read.csv('gauge_data_coords/gauge_coords.csv')
    gauge_coords = data.frame(lon=gauge_coords$lon, lat=gauge_coords$lat, id=200000 + 1:nrow(gauge_coords) + 0.2)    

    # QLD tide gauges
    qld_gauges = read.csv('qld_gov_locations/gauge_coords.csv')
    qld_gauges = data.frame(lon=qld_gauges$lon, lat=qld_gauges$lat, id=300000 + 1:nrow(qld_gauges) + 0.3)

    # Points of interest around Gladstone itself
    gladstone_coords = read.csv('gladstone_bay/gladstone_coords.csv')
    gladstone_coords = data.frame(lon=gladstone_coords$lon, lat=gladstone_coords$lat, id=400000 + 1:nrow(gladstone_coords) + 0.4)

    # Historic QLD tide gauges
    qld_historic_gauges = read.csv('/g/data/w85/tsunami/DATA/TIDES/QLD_TIDAL_GAUGES/gauge_locations/ArlieBeach_to_FraserIsand.csv')
    qld_historic_gauges = data.frame(lon=qld_historic_gauges$Instrument.Longitude.Decimal, lat=qld_historic_gauges$Instrument.Latitude.Decimal, id=500000 + 1:nrow(qld_historic_gauges) + 0.5)
    
    # Combine
    ALL_GAUGES = rbind(hthh_gauges, gauge_coords, qld_gauges, gladstone_coords, qld_historic_gauges)
    
    # Limited to QLD coast region bounding box:
    k = which(ALL_GAUGES$lon > bounding_box$lon_range[1] & ALL_GAUGES$lon < bounding_box$lon_range[2] &
              ALL_GAUGES$lat > bounding_box$lat_range[1] & ALL_GAUGES$lat < bounding_box$lat_range[2])

    ALL_GAUGES = ALL_GAUGES[k,]

    tokeep = find_close_locations(ALL_GAUGES, threshold=10)
    ALL_GAUGES = ALL_GAUGES[which(tokeep),]

    # Consider adding locations shifted by in neighbouring cell directions. This is to
    # protect against the situation where a model discretization puts a site on
    # land, but there is a nearby site that is water
    elev_rast_path = '/g/data/w85/tsunami/MODELS/inundation/east_australian_coast_2021_02/swals/rasters/Fine_KT43731_EdenToFraser_120321-full-ambient_sea_level_1.0/RUN_20210312_182512362/elevation_all.vrt' # File on NCI
    distoffset = 65 # m
    elev_rast <- rast(elev_rast_path)
    for(i in 1:9){ 
        # 3x3 grid, centre = gauge
        # Position defined by lx, ly
        lx = (i-1)%%3 -1 # -1, 0, 1
        ly = floor((i-1)/3) - 1 # -1, 0, 1 for every value of lx
        if((lx == 0) & (ly == 0)) next

        deg = atan2(ly, lx)/pi*180 # Bearing to shift point
        newpt = destPoint(as.matrix(ALL_GAUGES[,1:2]), b=deg, d=distoffset, f=0) # New point, spherical earth has f=0
        elev = extract(elev_rast, newpt) # Elevation at the new point

        if(i == 1){
            keep_pt = newpt
            keep_elev = elev[,1]
        }else{
            # Update the value of keep_pt and keep_elev if the elevation is
            # below the minimum so far.
            # This means we are keeping the deepest point (judged more likely to be in
            # a deep region that won't be spuriously dry on a coarse grid.)
            k = which(elev[,1] < keep_elev)
            if(length(k) > 0){
                keep_pt[k,] = newpt[k,]
                keep_elev[k] = elev[k,1]
            }
        }
    }

    ALL_GAUGES = rbind(ALL_GAUGES, data.frame(lon=keep_pt[,1], lat=keep_pt[,2], id=1000000 + ALL_GAUGES$id))
    return(ALL_GAUGES)
}

get_relevant_PTHA18_points<-function(bounding_box){
    # Get the PTHA access routines, where on NCI or at home
    ptha_access_script  = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R'
    ptha_env = new.env()
    source(ptha_access_script, local=ptha_env, chdir=TRUE)

    ALL_GAUGES = ptha_env$get_all_gauges()
    # Limited to QLD Gladstone coast region (roughly Mackay to Coffs harbour)
    # as a box, less generously than for nearshore tide gauges
    k = which(ALL_GAUGES$lon > bounding_box$lon_range[1] & ALL_GAUGES$lon < bounding_box$lon_range[2] &
              ALL_GAUGES$lat > bounding_box$lat_range[1] & ALL_GAUGES$lat < bounding_box$lat_range[2])

    ALL_GAUGES = ALL_GAUGES[k,]

    output = data.frame(lon=ALL_GAUGES$lon, lat=ALL_GAUGES$lat, id=ALL_GAUGES$gaugeID)
    return(output)
}

# Create a bounding box named list for the region of interest
bounding_box <- list(lon_range=c(147, 156), lat_range=c(-29, -17))

ALL_DARTS = merge_darts()
ALL_GAUGES = merge_coastal_tidegauges(bounding_box)
PTHA_gauges = get_relevant_PTHA18_points(bounding_box)

# When combining the points, make the IDs for non-PTHA18 points negative to avoid clashes
ALL_GAUGES$id = -1*ALL_GAUGES$id
ALL_DARTS$id = -1*ALL_DARTS$id

# combine all points
ALL_PTS = rbind(PTHA_gauges, ALL_DARTS, ALL_GAUGES)

# remove duplicates and log which points were removed
duplicates = ALL_PTS[duplicated(ALL_PTS[,1:2]),]
if(nrow(duplicates) > 0){
    warning_msg <- paste('Warning: Duplicates found and removed:\n', paste(duplicates$lon, duplicates$lat, duplicates$id, collapse=',\n'))
    print("Duplicates found and removed. Check the log file.")
    # log the warning
    log_filename <- "make_gauge_locations.log"
    fid <- file(log_filename, open="a")
    writeLines(warning_msg, fid)
    close(fid)
}
ALL_PTS = ALL_PTS[!duplicated(ALL_PTS[,1:2]),]


# write all points to file
write.csv(ALL_PTS, file='locations.csv', row.names=FALSE)

# log the time
log_filename <- "make_gauge_locations.log"
fid <- file(log_filename, open="a")
today <- format(Sys.time(), "%Y-%m-%d")
writeLines(today, fid)
close(fid)

# Check for unique IDS
if(length(ALL_PTS$id) != length(unique(ALL_PTS$id))){
    warning_msg <- 'Warning: IDS are not unique'
    print(warning_msg)
    # log the warning
    fid <- file(log_filename, open="a")
    writeLines(warning_msg, fid)
    close(fid)
}
