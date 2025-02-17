#
# Get JATWC parameters for a single scenario [typically a source inversion]
# - This is ONLY used to compare JATWC warning statistics to an inverted scenario.
# - The code used to compute JATWC-warning zones from hazard scenarios is in ../analysis_..../jatwc_to_inundation/
#

library(sf)
library(raster)

#
# INPUTS
#
asfm = new.env()
source('application_specific_metadata.R', local=asfm)
source('make_vrt.R')

# Points in ATWS zones used to compute their offshore wave height statistic
offshore_points_RDS = Sys.glob('Inundation_zones/*/comparison*.RDS')
atws_zones_to_compute = basename(dirname(offshore_points_RDS))

source_inversion_max_stage_tifs = list(
    'Puysegur2009' = "../../swals/OUTPUTS/run_Puysegur2009_PTHA18HS1788_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20240813_130956242/max_stage*.tif",
    'Chile2010'    = "../../swals/OUTPUTS/run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_183023159/max_stage*.tif",
    'Tohoku2011'   = "../../swals/OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/max_stage*.tif")

# Overwrite the initial stage [since it was 0 for the source inversions, different to hazard scenarios]
asfm$SCENARIO_AMBIENT_SEA_LEVEL = 0.0

# Sanity check
file_counts = unlist(lapply(source_inversion_max_stage_tifs, function(x) length(Sys.glob(x))))
if(any(file_counts == 0)) stop('Could not find max-stage files')

# Create vrt files holding max-stage for each scenario
source_inversion_max_stage = lapply(source_inversion_max_stage_tifs, function(x) raster(make_vrt(x)))
names(source_inversion_max_stage) = names(source_inversion_max_stage_tifs)
#
# Get the 95th percentile max-stage at sample points in the warning zone,
# excluding dry sites and locations that the tsunami never reaches.
#
get_max_stage_percentile_in_zone<-function(scenario_max_stage, sample_points_file, prob=0.95){

    sample_points = readRDS(sample_points_file)

    # Extract stage/elevation at the sampled points
    sample_stage = extract(scenario_max_stage, sample_points[,1:2])
    sample_elev = sample_points[,3]

    # Exclude points with NA stage, or dry points, or points that never exceed their initial stage
    keep = which( (sample_stage > (sample_elev + asfm$WETTOL)) &
                  (sample_stage > (asfm$SCENARIO_AMBIENT_SEA_LEVEL + asfm$WETTOL)) &
                  !is.na(sample_stage))

    # JATWC H parameter
    jatwc_H_parameter = quantile(sample_stage[keep], probs=prob, type=6) - asfm$SCENARIO_AMBIENT_SEA_LEVEL

    rm(sample_stage, sample_elev, keep); gc()

    return(c(jatwc_H_parameter))
}

# Compute results for each source inversion
JATWC_H = vector(mode='list', length=length(source_inversion_max_stage))
names(JATWC_H) = names(source_inversion_max_stage)
for(i in 1:length(source_inversion_max_stage)){
    source_inversion = names(source_inversion_max_stage)[i]
    scenario_max_stage = source_inversion_max_stage[[i]]
    JATWC_H[[source_inversion]] = vector(mode='list', length=length(atws_zones_to_compute))
    names(JATWC_H[[source_inversion]]) = atws_zones_to_compute

    # Compute results for each ATWS czone
    for(j in 1:length(atws_zones_to_compute)){
        JATWC_H[[source_inversion]][[atws_zones_to_compute[j]]] = 
            get_max_stage_percentile_in_zone(scenario_max_stage, offshore_points_RDS[j])
    }
}

saveRDS(JATWC_H, 'test_scenario_JATWC_statistics.RDS')

y = lapply(JATWC_H, unlist)
z = do.call(cbind, y)
write.csv(z, 'test_scenario_JATWC_statistics.csv')
