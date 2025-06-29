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
dem_file_nci = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
GLOBAL_DEM = raster(dem_file_nci)

# To align the nesting correctly it is useful to know this
GLOBAL_DOMAIN_LOWER_LEFT = c(100, -79.0)

# Split the global domain into this many pieces
GLOBAL_DOMAIN_NUM_PARTITIONS = 48

# Buffer domains by this fraction of their longest dimension before searching
# for the depth range (which influences timestepping). This mimics halos.
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
outdir_startname = 'domains_'

#
# Define the nesting regions with multiple shapefiles. 
# Each shapefile should have one or more simple polygons (no interior holes,
# although in-principle that case could be treated) 
#
NUMBER_NESTING_LEVELS = 4 # Not counting global domain
nesting = vector(mode='list', length=NUMBER_NESTING_LEVELS)
for(i in 1:length(nesting)){
    # For each level of nesting, store the SpatialPolygonsDataFrame defining
    # its region, and the desired size of each domain.
    nesting[[i]] = vector(mode='list', length=2)
    names(nesting[[i]]) = c('region', 'domain_size')
}

# One level finer than global domain
nesting[[1]]$region = readOGR('input/nesting_level_1.shp', layer='nesting_level_1')
# Tile size (extent of one domain prior to halo buffering)
nesting[[1]]$domain_size = c(1.0, 1.0)
names(nesting)[1] = 'nesting_level_1'

# Two levels finer than global domain (typically inside first level)
nesting[[2]]$region = readOGR('input/nesting_level_2.shp', layer='nesting_level_2')
# Must be an integer divisor of the coarser domain size.
# Also, to ensure nesting behaves, each 'piece' should contain an integer number of coarser domain cells.
nesting[[2]]$domain_size = nesting[[1]]$domain_size/3
names(nesting)[2] = 'nesting_level_2'

# Three levels finer than global domain (typically inside second level)
nesting[[3]]$region = readOGR('input/nesting_level_3.shp', layer='nesting_level_3')
# Must be an integer divisor of the coarser domain size.
# Also, to ensure nesting behaves, each 'piece' should contain an integer number of coarser domain cells.
nesting[[3]]$domain_size = nesting[[2]]$domain_size/5
names(nesting)[3] = 'nesting_level_3'

## Fourth level is inside the second level, not the third.
nesting[[4]]$region = readOGR('input/nesting_level_4.shp', layer='nesting_level_4')
nesting[[4]]$domain_size = nesting[[2]]$domain_size/5 # Must be an integer divisor of the coarser domain size. Also,
                                                     # to ensure nesting behaves, each 'piece' should contain an 
                                                     # integer number of coarser domain cells.
names(nesting)[4] = 'nesting_level_4'
#
###...more nesting levels here.

nesting_ID = paste(unlist(lapply(nesting, function(x) (1.0/x$domain_size[1]))), collapse='_')
outdir = paste0(outdir_startname, nesting_ID)
dir.create(outdir, showWarnings=FALSE)

##
## END INPUTS
##


#
# Make the domains, starting with finer nesting regions
#
for(i in seq(length(nesting), 1, by=-1)){

    if(i == length(nesting)){
        # Separate individual polygons in region_polys
        region_polys = lapply(nesting[[i]]$region@polygons, function(x) x@Polygons[[1]]@coords)
    }else{
        # For all but the finest level of nesting, merge with a buffered
        # version of the finer nesting region ("region_sp" exists from the last loop iteration)
        # This ensures that the level-i nesting box boundaries are not too close to the
        # level-(i+1) nesting box boundaries.
        tmp = gUnion(nesting[[i]]$region, 
                     gBuffer(region_sp, width=max(nesting[[i+1]]$domain_size)*domain_buffer_factor), 
                     byid=TRUE)
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

    # NOTE: Consider removing boxes that are entirely covered by finer domains

    # Extract min raster elevations in each domain (with some buffering to
    # account for nesting halos) by converting to spatialpolygons
    region_sp = domain_boundaries_to_spatialpolygons(region_boundaries, ID_flag=paste0('N', i, '_'))
    proj4string(region_sp) = CRS(proj4string(nesting[[i]]$region))
    region_elev_minmax = extract(GLOBAL_DEM, gBuffer(region_sp, width=max(region_size)*domain_buffer_factor, byid=TRUE), fun=range)

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
