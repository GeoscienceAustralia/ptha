# Split up jobs to compute thresholds corresponding to epistemic uncertainty exceedance-rates

# Do calculations for domains between these indices (inclusive)
MIN_INDEX_ALL = 2
MAX_INDEX_ALL = 133

# How many jobs to create?
N_JOBS = 40 # 8

library(parallel)
inds = MAX_INDEX_ALL - MIN_INDEX_ALL + 1
jobinds = splitIndices(inds, N_JOBS)
job_minmax = lapply(jobinds, function(x) c(min(x), max(x)))

# inds = c(56, 60, 72, 78)
# nx = length(inds)
# jobinds = splitIndices(nx, N_JOBS)
# jobinds = inds[as.numeric(jobinds)]
# job_minmax = lapply(jobinds, function(x) c(min(x), max(x)))


# Loop over multiple template_script_name values
flow_variables = c('min_stage', 'max_stage', 'max_speed', 'max_flux', 'depth')

template_script_name = 'TEMPLATE_run_compute_thresholds_at_exceedance_rate_of_epistemic_uncertainty_percentile.txt'

for(variable in flow_variables){
    for(i in 1:length(job_minmax)){

        min_domain_index = job_minmax[[i]][1]
        max_domain_index = job_minmax[[i]][2]

        # Name of the script, modified to contain the min/max domain indices in the filename
        output_filename = template_script_name
        output_filename = gsub('TEMPLATE_', '', output_filename)
        output_filename = gsub('.txt', '', output_filename)
        output_filename = paste(output_filename, variable, sep='_')
        output_filename = paste(output_filename, min_domain_index, max_domain_index, sep='_')
        output_filename = paste(output_filename, '.pbs', sep='')

        # Make the script contents, modified to refer to the min/max domain indices for the i'th job.
        template_script = readLines(template_script_name)
        output_script = template_script
        output_script = gsub('__MIN_DOMAIN_INDEX__', min_domain_index, output_script)
        output_script = gsub('__MAX_DOMAIN_INDEX__', max_domain_index, output_script)
        output_script = gsub('__FLOW_VARIABLE__', variable, output_script)

        # Write the modified output script to a file
        cat(output_script, file=output_filename, sep="\n")
    }
}
