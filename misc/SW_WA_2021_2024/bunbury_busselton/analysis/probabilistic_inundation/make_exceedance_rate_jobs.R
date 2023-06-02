template_file = 'run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh'
template_script = readLines(template_file)

rundir = 'ptha18-BunburyBusseltonRevised-sealevel60cm'

# Determine the number of domains in the multidomain by counting depth raster files in one case.
testfile = Sys.glob(paste0('../../swals/OUTPUTS/', rundir, '/random_outerrisesunda/*/raster_output_files.tar'))[1]
max_domains = as.numeric(system(paste0('tar --list -f ', testfile, ' | grep "depth_" | wc -w'), intern=TRUE))

# Define sets of domains by upper/lower indices. Each set will be run on a separate job.
# Notice we skip domain 1 -- it is too large, we get memory failures. Solution
# would be to not merge subdomains during raster creation.
NJOBS = 12 # Number of jobs (for each percentile)
dbounds = round(seq(1, max_domains, length=NJOBS))
n = length(dbounds)
uppers = dbounds[-1]
lowers = dbounds[-n] + 1  # Skip domain 1

percentiles = c("0.16", "0.84")

for(pp in percentiles){
    for(i in 1:length(lowers)){

        script = template_script
        script = gsub('_PERCENTILE_', pp, script)
        script = gsub('_LOWER_', lowers[i], script)
        script = gsub('_UPPER_', uppers[i], script)
        script = gsub('_RUNDIR_', rundir, script)
        

        outfile = template_file
        outfile = gsub('PERCENTILE', pp, outfile)
        outfile = gsub('LOWER', lowers[i], outfile)
        outfile = gsub('UPPER', uppers[i], outfile)
        outfile = gsub('RUNDIR', rundir, outfile)

        cat(script, file=outfile, sep="\n")

    }
}
