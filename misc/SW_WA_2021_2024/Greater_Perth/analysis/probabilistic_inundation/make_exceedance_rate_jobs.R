template_file = 'run_compute_exceedance_rates_at_epistemic_uncertainty_RUNDIR_PERCENTILE_LOWER_UPPER.sh'
template_script = readLines(template_file)

print('CHECK THAT THE LARGEST VALUE OF uppers IS AS DESIRED')
print('In future edit the code to derive the required number')

lowers = c(2 , 51, 101, 151, 201, 251, 301, 351, 401, 451, 501, 544) # NOTE: Added small jobs to the end (missed by accident)
uppers = c(50,100, 150, 200, 250, 300, 350, 400, 450, 500, 543, 559) # In future, use code to get the largest value of 'uppers'
percentiles = rep('0.16', length(lowers))
rundir = 'ptha18-GreaterPerth-sealevel60cm-reviseddomain-highres'

for(i in 1:length(lowers)){

    script = template_script
    script = gsub('_PERCENTILE_', percentiles[i], script)
    script = gsub('_LOWER_', lowers[i], script)
    script = gsub('_UPPER_', uppers[i], script)
    script = gsub('_RUNDIR_', rundir, script)
    

    outfile = template_file
    outfile = gsub('PERCENTILE', percentiles[i], outfile)
    outfile = gsub('LOWER', lowers[i], outfile)
    outfile = gsub('UPPER', uppers[i], outfile)
    outfile = gsub('RUNDIR', rundir, outfile)

    cat(script, file=outfile, sep="\n")

}
