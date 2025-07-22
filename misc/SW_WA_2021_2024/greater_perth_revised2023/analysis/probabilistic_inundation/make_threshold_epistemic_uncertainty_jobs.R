# Split up jobs to compute thresholds corresponding to epistemic uncertainty exceedance-rates

# Do calculations for domains between these indices (inclusive)
MIN_INDEX_ALL = 293
MAX_INDEX_ALL = 571

# How many jobs to create?
N_JOBS = 20 # 8

library(parallel)
inds = MAX_INDEX_ALL - MIN_INDEX_ALL + 1
jobinds = splitIndices(inds, N_JOBS)
job_minmax = lapply(jobinds, function(x) range(x) + MIN_INDEX_ALL - 1)

#
# Th template_script_name file contains placeholders for the min_domain_index
# (__MIN_DOMAIN_INDEX__) and the max_domain_index (__MAX_DOMAIN_INDEX__).
#
# This makes it easy to edit in an automated way
#
# template_script_name = 'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_DEPTH___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh' # Depth
#template_script_name = 'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_STAGE___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh' # Max stage
# template_script_name = 'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_FLUX___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh' # Flux
# template_script_name = 'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_SPEED___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh' # Speed
#
# Loop over multiple template_script_name values
all_template_scripts = c(
    # Depth
    'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_DEPTH___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh',
    # max stage 
    'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_STAGE___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh',
    # max flux
    'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_FLUX___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh',
    # max speed
    'run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile_MAX_SPEED___MIN_DOMAIN_INDEX_____MAX_DOMAIN_INDEX__.sh')

for(template_script_name in all_template_scripts){

    template_script = readLines(template_script_name)

    for(i in 1:length(job_minmax)){

        min_domain_index = job_minmax[[i]][1]
        max_domain_index = job_minmax[[i]][2]

        # Name of the script, modified to contain the min/max domain indices in the filename
        output_filename = template_script_name
        output_filename = gsub('__MIN_DOMAIN_INDEX__', min_domain_index, output_filename)
        output_filename = gsub('__MAX_DOMAIN_INDEX__', max_domain_index, output_filename)

        # Make the script contents, modified to refer to the min/max domain indices for the i'th job.
        output_script = template_script
        output_script = gsub('__MIN_DOMAIN_INDEX__', min_domain_index, output_script)
        output_script = gsub('__MAX_DOMAIN_INDEX__', max_domain_index, output_script)

        # Write the modified output script to a file
        cat(output_script, file=output_filename, sep="\n")
    }
}
