#
# For each site we want
# - Pairwise plots of 'depth' at the site using both 'elevation variation' and 'real static sea level' for every scenario
# - Scatterplot [with all scenarios] of the depth at a chosen site
#
source('sites_at_sea_levels.R')

# Get the SWALS post-processing scripts
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
swals = new.env()
source(file_nci, local=swals)

# Scenarios with elevation varying to represent the sea level (so the model's
# sea level of zero corresponds to HAT)
scenarios_sea_level_vary = Sys.glob('../../swals/OUTPUTS/ptha18-kalbarri2onslow-tidal_testing/random_*/*ambient_sea_level_0')

# What other static sea levels did we run the scenarios at?
static_sea_levels = c(0.65, 0.95, 1.17, 1.42)

OUTPUT_BASEDIR = './rasters_at_sites'

# Some scenarios went unstable at one or more sea levels (because I wasn't
# careful enough to give them a short timestep, as was done in the production
# runs to work around these issues). Remove them -- we'll still have enough
# scenarios for this test.
remove_scenarios_that_failed_in_some_cases<-function(scenarios_sea_level_vary){

    scenarios_to_drop = c(
        'sunda2_row_0107985_Mw_94_HS',
        'sunda2_row_0108100_Mw_94_HS',
        'sunda2_row_0110660_Mw_96_HS',
        'sunda2_row_0110734_Mw_96_HS',
        'sunda2_row_0110907_Mw_96_HS',
        'sunda2_row_0109369_Mw_95_HS')

    scenario_keep_flag = rep(TRUE, length(scenarios_sea_level_vary))
    for(i in 1:length(scenarios_to_drop)){
        k = which(grepl(scenarios_to_drop[i], scenarios_sea_level_vary, fixed=TRUE))
        if(length(k) > 0) scenario_keep_flag[k] = FALSE
    }

    scenarios_sea_level_vary[scenario_keep_flag]
}
scenarios_sea_level_vary = remove_scenarios_that_failed_in_some_cases(scenarios_sea_level_vary)

# Given a simulation at a site, extract the relevant info
process_simulation<-function(scenario_i, domain_indices, site_name, sea_level){

        stopifnot(length(scenario_i) == 1)
        print(paste0('Processing ', scenario_i))

        # Find relevant grid.nc files for varying elevation runs
        # Make a raster and store it locally
        grid_nc_files = Sys.glob(paste0(scenario_i, '/RUN*/RUN_ID*00000', domain_indices, '_*/Grid*.nc'))
        if(length(domain_indices) != length(grid_nc_files)){
            print(grid_nc_files)
            stop('Did not find the right grid nc files for domain_indices ', paste0(domain_indices, collapse=" "))
        }
        elevation_vary = swals$merge_domains_nc_grids(grid_nc_files, desired_var='elevation0', return_raster=TRUE, proj4string='EPSG:4326')
        max_stage_vary = swals$merge_domains_nc_grids(grid_nc_files, desired_var='max_stage', return_raster=TRUE, proj4string='EPSG:4326')
        depth_vary = max_stage_vary - elevation_vary

        # Save the files
        if(sea_level == 0){
            dirtag = '/varying_elevation/'
        }else{
            dirtag = paste0('/static_sea_level_', sea_level, '/')
        }
        output_folder = paste0(OUTPUT_BASEDIR, '/', site_name, dirtag)
        dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)
        output_scenario_info = paste0(basename(scenario_i), '_', site_name)
        writeRaster(elevation_vary, paste0(output_folder, output_scenario_info, '_elevation0.tif'))
        writeRaster(max_stage_vary, paste0(output_folder, output_scenario_info, '_max_stage.tif'))
        writeRaster(depth_vary, paste0(output_folder, output_scenario_info, '_depth.tif'))

}

# Extract the relevant rasters for a site
make_rasters_for_plotting<-function(site_name, site_data, scenarios_sea_level_vary){

    sea_level = site_data$sealevel
    plot_bbox = matrix(site_data$plot_bbox, ncol=2, byrow=TRUE)
    domain_indices = site_data$domain_indices

    # Filenames for varying elevation runs
    scenarios_sea_level_vary = scenarios_sea_level_vary 
    # Filenames for static sea level 
    scenarios_sea_level_static = gsub("ambient_sea_level_0", paste0('ambient_sea_level_', sea_level), scenarios_sea_level_vary)

    for(i in 1:length(scenarios_sea_level_vary)){
        process_simulation(scenarios_sea_level_vary[i], domain_indices, site_name, sea_level=0)
        process_simulation(scenarios_sea_level_static[i], domain_indices, site_name, sea_level)
    }

    return(0)
}

# Make the rasters for all sites. This takes a few minutes, could easily
# be run in parallel but not yet.
for(site_name in names(sites_at_sea_levels)){
    nothing = make_rasters_for_plotting(site_name, sites_at_sea_levels[[site_name]], scenarios_sea_level_vary)
}
