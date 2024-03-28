##
## Make gauge locations for the NSW model
##
## Include gauges at
# - DARTS
# - NZ DARTS
# - Tide gauges in eastern Australia [HTHH paper + other TG data]
# - PTHA18 points near eastern Australia (but not all global pts, to manage file sizes)
# - Other?
library(terra)
library(geosphere)

find_close_locations<-function(ALL_DARTS, threshold){
    # Handy routine for combining gauges.
    # ALL_DARTS must have columns lon/lat, giving the location in degrees
    # Find a subset of sites in ALL_DARTS such that the full set of sites
    # is always within a distance "threshold" of the subset. 
    # Result is tokeep, an array of row indices from ALL_DARTS to keep
    tokeep = rep(NA, nrow(ALL_DARTS))
    tokeep[1] = TRUE
    for(i in 2:nrow(ALL_DARTS)){
        w_t_k = which(tokeep)
        p0 = cbind(ALL_DARTS$lon[w_t_k], ALL_DARTS$lat[w_t_k]) # Points we know we will keep
        p1 = p0*0 # point i
        p1[,1] = ALL_DARTS[i,1]
        p1[,2] = ALL_DARTS[i,2]
        dists = distHaversine(p0, p1)
        tokeep[i] = all(dists > threshold)
    }
    return(tokeep)
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

    ALL_DARTS = rbind(nz_darts, darts_2, darts_hr)

    tokeep = find_close_locations(ALL_DARTS, threshold=5000)
    ALL_DARTS = ALL_DARTS[which(tokeep),]

    return(ALL_DARTS)
}
ALL_DARTS = merge_darts()


#
# Coastal tide gauges
#
merge_coastal_tidegauges<-function(){
    # HTHH paper locations
    hthh_gauges = read.csv('hthh_paper_tide_gauges/01_tide_gauge_locations.csv')
    hthh_gauges = data.frame(lon=hthh_gauges$lon, lat=hthh_gauges$lat, id=100000 + 1:nrow(hthh_gauges) + 0.1)

    # Tide gauges from my gauge data
    gauge_coords = read.csv('gauge_data_coords/gauge_coords.csv')
    gauge_coords = data.frame(lon=gauge_coords$lon, lat=gauge_coords$lat, id=200000 + 1:nrow(gauge_coords) + 0.2)    

    # MHL tidal plane sites
    mhl_tidal_planes = read.csv('mhl_tidal_plane_coords/MHL_Tidal_Planes_combined_summary.csv')
    mhl_tidal_planes = data.frame(lon=mhl_tidal_planes$lon, lat=mhl_tidal_planes$lat, id=300000 + 1:nrow(mhl_tidal_planes) + 0.3)

    # New gauges covering Port Authority + additional Lord Howe Island locations
    pa_lhi_gauges = vect('extra_gauges_2023_11_14/extra_gauges.shp')
    pa_lhi_gauge_coords = crds(pa_lhi_gauges)
    pa_lhi_coords = data.frame(lon=pa_lhi_gauge_coords[,1], lat=pa_lhi_gauge_coords[,2], id=400000 + 1:nrow(pa_lhi_gauge_coords) + 0.4)
    

    # Combine
    ALL_GAUGES = rbind(hthh_gauges, gauge_coords, mhl_tidal_planes, pa_lhi_coords)
    
    # Limited to NSW coast region, as a box:
    # (149, -39) to (170, -27)
    k = which(ALL_GAUGES$lon > 149 & ALL_GAUGES$lon < 170 &
              ALL_GAUGES$lat > -39 & ALL_GAUGES$lat < -27)

    ALL_GAUGES = ALL_GAUGES[k,]

    tokeep = find_close_locations(ALL_GAUGES, threshold=10)
    ALL_GAUGES = ALL_GAUGES[which(tokeep),]

    # Consider adding locations shifted by in neighbouring cell directions. This is to
    # protect against the situation where a model discretization puts a site on
    # land, but there is a nearby site that is water
    nci_path = '/g/data/w85/tsunami/MODELS/inundation/east_australian_coast_2021_02/swals/rasters/Fine_KT43731_EdenToFraser_120321-full-ambient_sea_level_1.0/RUN_20210312_182512362/elevation_all.vrt' # File on NCI
    home_path = '/media/gareth/Data3/DATA/SWALS/2021_east_coast_model/swals/rasters/Fine_KT43731_EdenToFraser_120321-full-ambient_sea_level_1.0/RUN_20210312_182512362/elevation_all.vrt' # File at home
    elev_rast_path = ifelse(file.exists(nci_path), nci_path, home_path) # Works both on NCI, and at home
    elev_rast = rast(elev_rast_path)
    distoffset = 65 # m
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
ALL_GAUGES = merge_coastal_tidegauges()

get_relevant_PTHA18_points<-function(){
    # Get the PTHA access routines, where on NCI or at home
    ptha_access_script_NCI  = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
    ptha_access_script_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
    ptha_access_script = ifelse(file.exists(ptha_access_script_NCI), ptha_access_script_NCI, ptha_access_script_home)
    ptha_env = new.env()
    source(ptha_access_script, local=ptha_env, chdir=TRUE)

    ALL_GAUGES = ptha_env$get_all_gauges()
    # Limited to NSW coast region, as a box, less generously than for nearshore tide gauges:
    # (149, -39) to (157, -27)
    k = which(ALL_GAUGES$lon > 149 & ALL_GAUGES$lon < 157 &
              ALL_GAUGES$lat > -39 & ALL_GAUGES$lat < -27)

    ALL_GAUGES = ALL_GAUGES[k,]

    output = data.frame(lon=ALL_GAUGES$lon, lat=ALL_GAUGES$lat, id=ALL_GAUGES$gaugeID)
    return(output)
}
PTHA_gauges = get_relevant_PTHA18_points()


# When combining the points, make the IDs for non-PTHA18 points negative to avoid clashes
ALL_GAUGES$id = -1*ALL_GAUGES$id
ALL_DARTS$id = -1*ALL_DARTS$id

ALL_PTS = rbind(PTHA_gauges, ALL_DARTS, ALL_GAUGES)
#write.csv(ALL_PTS, file='point_gauges_2023_08_08.csv', row.names=FALSE)
write.csv(ALL_PTS, file='point_gauges_2023_11_14.csv', row.names=FALSE)

# Check for unique IDS
if(length(ALL_PTS$id) != length(unique(ALL_PTS$id))){
    print('Warning: IDS are not unique')
}
