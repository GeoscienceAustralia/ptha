#
# Create polygons depicting our inundation models for various JATWC zones, limited by a high PTHA percentile:
#     - No threat zone
#     - Marine and immediate foreshore threat zone
#     - Land threat zone
#     - [Not yet supported by JATWC] Minor land threat zone (0.55 < H < 1.5)
#     - [Not yet supported by JATWC] Major land threat zone (H > 1.5)
#
# We limit the zones with PTHA, because that gives a more stable upper limit than does choosing the largest (random) scenario.
#     - Mean or 84th percentile? Is there much difference at 1/10,000?
#         + From visual inspection of the Greater Perth model -- no, at 1/10,000 they tend
#         not to be very different. While minor differences are common, I
#         didn't see notice differences that I thought woudl affect the action.
#    - But at 1/2500 the differences are greater.
#
source('make_vrt.R')
library(terra)
library(stars)

#
# INPUTS
#
asfm = new.env()
source('application_specific_metadata.R', local=asfm)

# The ATWS Zone name should be set on the command line, and match a zone in the
# ATWS_ZONES shapefile for which we already did raster computations.
ATWS_ZONE_NAME = commandArgs(trailingOnly=TRUE)[1]

# The zone type should be set on the command line, and match one of the allowed types, like
# 'no_threat', 'marine_warning', 'land_warning', 'minor_land_warning', 'major_land_warning'
zone_type = commandArgs(trailingOnly=TRUE)[2]

# Max-stage values below 'AMBIENT_SEA_LEVEL + WET_TOL' are treated as dry. (WET_TOL > 0)
# can remove ponded areas that are flooded by the initial condtion, but never
# reached by the tsunami. However it cannot exclude such ponded areas if they
# are flooded by even 1 scenario in the JATWC warning category -- so some ponds
# are included, others excluded. On balance that can be confusing, so maybe
# best set to zero.
WET_TOL = 0.0 # Not the same concept as asfm$WETTOL

#
# END INPUTS
#

# Check the "ATWS Zone" name
ATWS_Zone_name_nospace = gsub(' ', '-', ATWS_ZONE_NAME)
working_dir = paste0('Inundation_zones/', ATWS_Zone_name_nospace)
stopifnot(file.exists(working_dir))

# Check the zone type
allowed_zone_types = names(asfm$JATWC_H_ranges) 
stopifnot(zone_type %in% allowed_zone_types)

# Initial model sea level
AMBIENT_SEA_LEVEL = asfm$SCENARIO_AMBIENT_SEA_LEVEL


# We once used a special treatment of the initial stage in the Vasse estuary
# FIXME: Generalise this in future applications if it is needed.
VASSE_ESTUARY_SPECIAL_TREATMENT = FALSE # FALSE will skip these calculations
if(VASSE_ESTUARY_SPECIAL_TREATMENT){
    VASSE_ESTUARY_POLYGON = vect('../../elevation/initial_stage_40cmAHD/initial_stage_40cmAHD.shp')
    VASSE_ESTUARY_INITIAL_STAGE = 0.4
}

# The upper limit of the zones is limited based on the PTHA inundation
# exceedance-rate results.
# This represents a probabilistic attempt to say "how big is the largest
# reasonable scenario" -- which must be set somehow to define the largest
# inundation zone. 
# Consider using either the logic-tree-mean, or the 84th percentile. 
ptha_exceedance_rate_rasts = asfm$ptha_exceedance_rate_rasts 
stopifnot(length(ptha_exceedance_rate_rasts) > 0)
# Make a name for the PTHA raster. Used in output folder name to give
# information on PTHA raster used. Avoid spaces
ptha_exceedance_rate_rast_tag = asfm$ptha_exceedance_rate_rast_tag # 'Logic_tree_mean'
# The exceedance-rate defining the upper zone.
PTHA_EXRATE_TOL = asfm$PTHA_EXRATE_TOL #1/10000

#
# END INPUTS
#

setwd(working_dir)

# Make a VRT with the JATWC zones defined from a union of all inundation scenarios.
# At this stage, there is no limitation from PTHA.
all_files = Sys.glob(paste0(zone_type, '_max_stage_domain_*.tif'))
stopifnot(length(all_files) > 1)

#jatwc_zones_all_scenarios = vrt(all_files, options=c('-resolution', 'highest'))
jatwc_zones_all_scenarios = rast(make_vrt(all_files))

# Get files in the PTHA raster that match these files.
# Not all are included, because the JATWC zone will not use all the model tifs
all_files_endofname = sapply(all_files, function(x) {y = strsplit(x, '_')[[1]]; return(y[length(y)])})
ptha_rasts_endofname = sapply(ptha_exceedance_rate_rasts, function(x) {y = strsplit(x, '_')[[1]]; return(y[length(y)])})
matching_files = match(ptha_rasts_endofname, all_files_endofname)
# Check we got a 1:1 match
stopifnot(all(sort(na.omit(matching_files)) == (1:length(all_files))  ))
k = which(!is.na(matching_files))

#ptha_vrt = vrt(ptha_exceedance_rate_rasts[k], options=c('resolution', 'highest'))
ptha_vrt = rast(make_vrt(ptha_exceedance_rate_rasts[k]))


# Logical raster that is true in the regions that should be included in the zone.
# This includes a PTHA based limit which ensures the 'largest reasonable' scenario
# is not a random quantity.
included_sites = ( jatwc_zones_all_scenarios >= (AMBIENT_SEA_LEVEL + WET_TOL) ) & (ptha_vrt >= PTHA_EXRATE_TOL)

if(VASSE_ESTUARY_SPECIAL_TREATMENT){
    # Special treatment for Vasse estuary, where we set the initial sea-level to a different value
    vasse_sites = rasterize(VASSE_ESTUARY_POLYGON, jatwc_zones_all_scenarios, background=0)
    included_sites = included_sites | (
        (vasse_sites == 1) & 
        (jatwc_zones_all_scenarios >= (VASSE_ESTUARY_INITIAL_STAGE + WET_TOL)) & 
        (ptha_vrt >= PTHA_EXRATE_TOL)) 
}

# FALSE --> NA
included_sites[!included_sites] = NA

zone_polygon = as.polygons(included_sites)

# Make the output directory have an informative name
output_dir = paste0(ATWS_Zone_name_nospace, '_', zone_type, '_with-PTHA-exrate-limit_', 
                    ptha_exceedance_rate_rast_tag, '_', PTHA_EXRATE_TOL) 
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

writeVector(zone_polygon, file=paste0(output_dir, '/', basename(output_dir), '.shp'), overwrite=TRUE)
