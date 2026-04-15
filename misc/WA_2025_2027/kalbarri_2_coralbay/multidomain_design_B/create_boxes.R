#
# Try creating domain bounding boxes to cover a polygonal region
#
library(rptha)

# Utilities to cover a polygon with boxes (to easily place domains within a
# desired irregular region).
source('_aggregate_boxes.R') 

##
## INPUTS
##

# Get a DEM covering all regions. We use this to estimate the min elevation in each
# domain, which is useful to guide local time-stepping.
dem_file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
dem_file_nci = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
GLOBAL_DEM = raster(ifelse(file.exists(dem_file_nci), dem_file_nci, dem_file_home))

# To align the nesting correctly it is useful to know this
GLOBAL_DOMAIN_LOWER_LEFT = c(10, -79)

# Split the global domain into this many pieces. The current code doesn't do
# the split, but writes a load balance file that will cause it.
GLOBAL_DOMAIN_NUM_PARTITIONS = 48*4 #16

# Buffer domains by this fraction of their longest dimension before searching
# for the depth range (which influences timestepping). This mimics halos. 
# Also used to increase the size of the polygon used to create geometries by
# merging with a buffered version of the next-finest nesting polygon inside it.
domain_buffer_factor = 0.25 

# Do we include a version of outputs with aggregated domains (i.e. combining
# neighbouring boxes if they seem to have similar timestepping)? This needs work,
# to date when I've tried it was more efficient without.
# Probably the box-merging scheme needs to be improved (e.g. by scanning logfiles from an 
# initial model run to get the typical timestep that the domain can take).
include_box_aggregation = FALSE

# If include_box_aggregation = TRUE, then when aggregating domains we enforce
# the following minimum depth (irrespective of what was found). This accounts
# for the fact that here we are only using global data to estimate elevation
# ranges -- but the model will use higher-res data.
aggregation_min_depth = 30

# Outputs will be stored in a folder having named that begins with the following
outdir_startname = 'kalbarri2coralbay_B_domains_20260324_'

#
# Define the nesting regions with multiple shapefiles. 
# Each shapefile should have one or more simple polygons (no interior holes,
# although in-principle that case could be treated) 
#
NUMBER_NESTING_LEVELS = 5 # Not counting global domain
nesting = vector(mode='list', length=NUMBER_NESTING_LEVELS)
for(i in 1:length(nesting)){
    # For each level of nesting, store
    # - SpatialPolygonsDataFrame defining its region, 
    # - the desired size of each domain (for when we convert the region to boxes) 
    # - the index (in nesting[[ ]]) of another geometry that is inside this
    #   level and has finer resolution (so we can remove boxes in this geometry 
    #   geometry that are completely covered by finer ones).
    # Also make workspace (initial_region_sp)
    nesting[[i]] = vector(mode='list', length=4)
    names(nesting[[i]]) = c('region', 'domain_size', 'buffer_with_index', 'initial_region_sp')

    # Union nesting[[i]]$region with a buffered version of geometries in
    # nesting[[buffer_with_index]] before computing boxes, unless
    # buffer_with_index <= 0. 
    # If buffer_with_index is positive then it must be > i (since we loop
    # over nesting regions in reverse order, and need to operate on finer
    # levels of nesting first)
    nesting[[i]]$buffer_with_index = ifelse(i < length(nesting), i+1, -1)
    # Store some modified geometries
    nesting[[i]]$initial_region_sp = NA
}

#
# Kalbarri to Onslow with Xmas and Cocos Islands
#

# Kalbarri to Onslow nesting
#nesting[[1]]$region = readOGR('kalbarri2onslow_C_level1/kalbarri2onslow_C_level1.shp', layer='kalbarri2onslow_C_level1')
#names(nesting)[1] = 'kalbarri2onslow_C_level1'
#nesting[[1]]$domain_size = c(3,3)

nesting[[1]]$region = readOGR('kalbarri2onslow_C_level2/kalbarri2onslow_C_level2.shp', layer='kalbarri2onslow_C_level2')
names(nesting)[1] = 'kalbarri2onslow_C_level2'
nesting[[1]]$domain_size = c(1,1) #nesting[[1]]$domain_size/3

nesting[[2]]$region = readOGR('kalbarri2onslow_C_level3/kalbarri2onslow_C_level3.shp', layer='kalbarri2onslow_C_level3')
names(nesting)[2] = 'kalbarri2onslow_C_level3'
nesting[[2]]$domain_size = nesting[[1]]$domain_size/3

nesting[[3]]$region = readOGR('kalbarri2onslow_C_level4/kalbarri2onslow_C_level4.shp', layer='kalbarri2onslow_C_level4')
names(nesting)[3] = 'kalbarri2onslow_C_level4'
nesting[[3]]$domain_size = nesting[[2]]$domain_size/7
nesting[[3]]$buffer_with_index = -1 # This is the finest domain, don't buffer with the next one before making the geometry

## Extras for Xmas/Cocos Island
#nesting[[5]]$region = readOGR('Xmas_cocos_level1/Xmas_cocos_level1.shp', layer='Xmas_cocos_level1')
#names(nesting)[5] = 'Xmas_cocos_level1'
#nesting[[5]]$domain_size = c(1,1) # Domains are smaller since our focus tends to be narrow

nesting[[4]]$region = readOGR('Xmas_cocos_level2/Xmas_cocos_level2.shp', layer='Xmas_cocos_level2')
names(nesting)[4] = 'Xmas_cocos_level2'
nesting[[4]]$domain_size = c(1,1) #nesting[[4]]$domain_size/3

nesting[[5]]$region = readOGR('Xmas_cocos_level3/Xmas_cocos_level3.shp', layer='Xmas_cocos_level3')
names(nesting)[5] = 'Xmas_cocos_level3'
nesting[[5]]$domain_size = nesting[[4]]$domain_size/5 # Keep this small since large domains spill into deep water (short timestep)
nesting[[5]]$buffer_with_index = -1

#
#
###...more nesting levels here.

##
## END INPUTS
##

nesting_ID = paste(unlist(lapply(nesting, function(x) (1.0/x$domain_size[1]))), collapse='_')
outdir = paste0(outdir_startname, nesting_ID)
dir.create(outdir, showWarnings=FALSE)

#
# Make the domains, starting with finer nesting regions
#
for(i in seq(length(nesting), 1, by=-1)){

    # Index in nestting[[ ]] for a geometry that is inside this one,
    # representing the next-finest nesting layer (or -1 if there is none)
    buffer_with_index = nesting[[i]]$buffer_with_index

    #
    # Get polygon coordinates; assumes each geometry has a single loop without holes.
    #
    if(buffer_with_index <= 0){
        region_polys = lapply(nesting[[i]]$region@polygons, function(x) x@Polygons[[1]]@coords)
    }else{
        stopifnot(buffer_with_index > i)

        # Merge region polygon with a buffered version of the finer nesting region.
        # This ensures that the level-i nesting box boundaries are not too
        # close to the level-(i+1) nesting box boundaries.
        tmp = gUnion(nesting[[i]]$region, 
                     gBuffer(nesting[[buffer_with_index]]$initial_region_sp, 
                            width=max(nesting[[buffer_with_index]]$domain_size)*domain_buffer_factor), 
                     byid=TRUE)
        # The first geometry in each Polygon is the one we want (the others
        # might relate to parts of bufered region_sp that don't actually touch
        # this part of nesting[[i]]$region, e.g. if the polygons are not perfectly nested)
        region_polys = lapply(tmp@polygons, function(x) x@Polygons[[1]]@coords)
    }
    region_size = nesting[[i]]$domain_size

    # Convert the polygon to a set of boxes    
    region_boundaries = lapply(region_polys, function(x){ 
        cover_polygon_with_rectangular_domains(x, domain_size=region_size, 
            grid_alignment_point=GLOBAL_DOMAIN_LOWER_LEFT, make_plot=FALSE)
        })
    region_boundaries = do.call(rbind, region_boundaries)
    # Ensure each box only occurs once (could be more than once with multi-part polygons).
    region_boundaries = unique(region_boundaries)

    # Extract min raster elevations in each domain (with some buffering to
    # account for nesting halos) by converting to spatialpolygons
    region_sp = domain_boundaries_to_spatialpolygons(region_boundaries, ID_flag=paste0('N', i, '_'))
    proj4string(region_sp) = CRS(proj4string(nesting[[i]]$region))
    region_elev_minmax = extract(GLOBAL_DEM, 
        gBuffer(region_sp, width=max(region_size)*domain_buffer_factor, byid=TRUE), 
        fun=range)

    # Useful to store the boxes before further editing, so we can use them to remove
    # boxes on other nesting levels.
    nesting[[i]]$initial_region_sp = region_sp

    if(buffer_with_index > 0){
        #
        # Remove boxes that are entirely covered by the geometry of nesting[[buffer_with_index]]
        #
        is_covered = rep(FALSE, length(region_sp))

        # Add a tiny (1/1000 x width) buffer to the finer domain before
        # checking if it covers, as gCovers might not detect exact overlap.
        buffered_initial_region_sp = gBuffer(nesting[[buffer_with_index]]$initial_region_sp,
            width = nesting[[buffer_with_index]]$domain_size*1e-05)

        # Check each box in region_sp to see if it is covered
        for(j in 1:length(region_sp)){
            is_covered[j] = gCovers(buffered_initial_region_sp, region_sp[j,])
        }

        # Remove regions that are covered
        k = which(!is_covered)
        if(length(k) == 0) stop(paste0('ERROR: nesting[[', i, ']] is eliminated by finer domains'))
        region_boundaries = region_boundaries[k,]
        region_sp = region_sp[k,]
        region_elev_minmax = region_elev_minmax[k,]
    }

    # Write-out for SWALS
    output_df = cbind(region_boundaries, 
        data.frame('min_elev_coarse_dem' = region_elev_minmax[,1], 'max_elev_coarse_dem' = region_elev_minmax[,2]))
    write.csv(output_df, file=paste0(outdir, '/', names(nesting)[i], '.csv'), row.names=FALSE, quote=FALSE)
    nesting[[i]]$output_df = output_df

    # Write to shapefile
    tmp = domain_boundaries_to_shapefile(output_df,
        file_dsn=paste0(outdir, '/', names(nesting)[i], '_domains'),
        file_layer=paste0(names(nesting)[i], '_domains'), 
        ID_flag=paste0('Level_', i, '_'))

    if(include_box_aggregation){
        # Make alternative version with aggregated boxes
        alternate_fine_boxes = aggregate_boxes(region_boundaries, region_elev_minmax[,1], depth_lower_bound=aggregation_min_depth)
        output_df = cbind(alternate_fine_boxes$boxes, data.frame('min_elev_coarse_dem' = alternate_fine_boxes$min_elev))
        write.csv(output_df, file=paste0(outdir, '/', names(nesting)[i], '_with_box_merging.csv'), row.names=FALSE, quote=FALSE)
        nesting[[i]]$output_df_box_merging = output_df

        # Write to shapefile
        tmp = domain_boundaries_to_shapefile(output_df,
            file_dsn=paste0(outdir, '/', names(nesting)[i], '_domains_with_box_merging'),
            file_layer=paste0(names(nesting)[i], '_domains_with_box_merging'), 
            ID_flag=paste0('Level_', i, '_'))
    }

}

#
# Make a default "load_balance_partition.txt" file. It will assume a single
# global domain [subdivided into one-piece-per-mpi-rank]. The other domains
# will not be subdivided. 
#
number_nested_domains = sum(unlist(lapply(nesting, function(x) nrow(x$output_df))))
file_lines = c(as.character(paste0(seq(1,GLOBAL_DOMAIN_NUM_PARTITIONS), collapse=" ")),
               # Do not split the other domains further [i.e. one number per row]
               seq(1, number_nested_domains))
cat(file_lines, file=paste0(outdir, '/', 'load_balance_default.txt'), sep="\n")

##
## Make a default load balance partition file for the 'merged_box' alternative
##
#number_nested_domains_merged = sum(unlist(lapply(nesting, function(x) nrow(x$output_df_box_merging))))
#file_lines = c(as.character(paste0(seq(1,GLOBAL_DOMAIN_NUM_PARTITIONS), collapse=" ")),
#               # Do not split the other domains further [i.e. one number per row]
#               seq(1, number_nested_domains_merged))
#cat(file_lines, file=paste0(outdir, '/', 'load_balance_default_box_merging.txt'), sep="\n")
