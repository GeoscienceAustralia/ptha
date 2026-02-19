#
# Spread the epistemic uncertainty calculations over multiple qsub scripts.
# This code creates the qsub scripts.
#

asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)

##
## INPUTS
##

# Number of qsub jobs for each percentile and source
NJOBS = 12 
# The template qsub file (that will be modified)
template_qsub_file = "compute_exceedance_rates_at_epistemic_uncertainty_VARIABLE_SOURCEZONE_PERCENTILE_LOWER_UPPER_EXCEEDANCETHRESHOLD.pbs"
# Source zones we model
source_zones = names(asfm$source_zone_modelled_tsunami_scenario_basedirs)
# A single raster tarfile (just so we can count the number of domains)
reference_raster_tarfile = asfm$reference_raster_tar_file
# A single variable
variable = 'depth'     # 'depth' or 'max_stage'
# The exceedance_threshold
exceedance_threshold = 0.001 
# The epistemic uncertainties percentiles to compute
percentiles = c("0.16", "0.84")
# Which domain to start from? Note domain 1 is often too big to run
# effectively
min_domains = 88 # Only include high-res domains

##
## END INPUTS
##

# Determine the number of domains in the multidomain by counting depth raster files in one case.
max_domains = as.numeric(system(paste0('tar --list -f ', reference_raster_tarfile, ' | grep "depth_" | wc -w'), intern=TRUE))

# Ensure there is enough work to avoid repeating runs
stopifnot(NJOBS <= max_domains/2)

# Define sets of domains by upper/lower domain indices. Each set will be run on a separate job.
# Notice we skip domain 1 -- it is too large, we get memory failures. Solution
# would be to not merge subdomains during raster creation.
dbounds = round(seq(min_domains, max_domains, length=NJOBS))
n = length(dbounds)
uppers = dbounds[-1]
lowers = dbounds[-n] + 1
lowers[1] = min_domains

stopifnot(length(unique(lowers)) == length(lowers))
stopifnot(length(unique(uppers)) == length(uppers))

template_script = readLines(template_qsub_file)

#
# Make a set of template scripts
#

for(sz in source_zones){
    random_scenarios_dir = asfm$source_zone_modelled_tsunami_scenario_basedirs[[sz]]

    for(pp in percentiles){
        # output_folder should be something like:
        # ptha18-BunburyBusseltonRevised-sealevel60cm-depth_exrate_0.001_0.84_outerrisesunda
        output_folder = paste0(
            basename(dirname(random_scenarios_dir)), '-', 
            variable, '_exrate_', exceedance_threshold, 
            '_', pp, '_', sz)

        for(i in 1:length(lowers)){

            script = template_script
            script = gsub('_VARIABLE_', variable, script)
            script = gsub('_EXCEEDANCETHRESHOLD_', exceedance_threshold, script)
            script = gsub('_PERCENTILE_', pp, script)
            script = gsub('_LOWER_', lowers[i], script)
            script = gsub('_UPPER_', uppers[i], script)
            script = gsub('_RANDOMSCENARIOS_', random_scenarios_dir, script, fixed=TRUE)
            script = gsub('_OUTPUTDIR_', output_folder, script)
            

            outfile = template_qsub_file
            outfile = gsub('VARIABLE', variable, outfile)
            outfile = gsub('SOURCEZONE', sz, outfile)
            outfile = gsub('PERCENTILE', pp, outfile)
            outfile = gsub('LOWER', lowers[i], outfile)
            outfile = gsub('UPPER', uppers[i], outfile)
            outfile = gsub('EXCEEDANCETHRESHOLD', exceedance_threshold, outfile)

            cat(script, file=outfile, sep="\n")

        }
    }
}
