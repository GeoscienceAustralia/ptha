# Setup threshold epistemic uncertainty jobs that failed the first attempt (because I didn't request enough walltime).

# Which variables to do calculation for
varnames = c('depth', 'max_stage', 'max_speed', 'max_flux')

# Where did the previous script save tifs for each variable
var_output_dirs = list(
    'depth' = 'kalbarri2coralbay_highres_domains_depth_percentile_0.84_exrate_0.0004_hazard',
    'max_flux' = 'kalbarri2coralbay_highres_domains_max_flux_percentile_0.84_exrate_0.0004_hazard',
    'max_speed' = 'kalbarri2coralbay_highres_domains_max_speed_percentile_0.84_exrate_0.0004_hazard',
    'max_stage' = 'kalbarri2coralbay_highres_domains_max_stage_percentile_0.84_exrate_0.0004_hazard')

# Do calculations for domains between these indices (inclusive)
MIN_INDEX_ALL = 2 # Avoid the global domain -- it's big and expensive and not very informative.
MAX_INDEX_ALL = 146 #

# Th template_script_name file contains placeholders for the min_domain_index
# (__MIN_DOMAIN_INDEX__) and the max_domain_index (__MAX_DOMAIN_INDEX__) and
# the variable name (__VARNAME__).
template_script_name = 'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile___VARNAME_____MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh'

for(varname in varnames){

    template_script = readLines(template_script_name)

    variable_outdir = var_output_dirs[[varname]] # The existing tifs for this variable should be inside here.

    for(i in MIN_INDEX_ALL:MAX_INDEX_ALL){

        # Search for the tif with this index. If it exists we don't have to do
        # anything
        this_case_worked_already = (length(grep(paste0('_domain_index_', i, '.tif'), Sys.glob(paste0(variable_outdir, '/*.tif')), fixed=TRUE)) > 0)
        if(this_case_worked_already) next

        min_domain_index = i #job_minmax[[i]][1]
        max_domain_index = i #job_minmax[[i]][2]

        # Name of the script, modified to contain the min/max domain indices in the filename
        output_filename = template_script_name
        output_filename = gsub('__VARNAME__', varname, output_filename)
        output_filename = gsub('__MIN_DOMAIN_INDEX__', min_domain_index, output_filename)
        output_filename = gsub('__MAX_DOMAIN_INDEX__', max_domain_index, output_filename)

        # Make the script contents, modified to refer to the min/max domain indices for the i'th job.
        output_script = template_script
        output_script = gsub('__VARNAME__', varname, output_script)
        output_script = gsub('__MIN_DOMAIN_INDEX__', min_domain_index, output_script)
        output_script = gsub('__MAX_DOMAIN_INDEX__', max_domain_index, output_script)

        # Write the modified output script to a file
        cat(output_script, file=output_filename, sep="\n")
    }
}
