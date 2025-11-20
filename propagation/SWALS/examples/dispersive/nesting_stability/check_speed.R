library(ncdf4)
source('../../../plot.R')
run_tag = commandArgs(trailingOnly=TRUE)[1]

# Get the simulations
md_dir = rev(sort(Sys.glob('OUTPUTS/RUN*')))[1]

# Get a log
md_log = Sys.glob(paste0(md_dir, '/*.log'))[1]

# Get the maximum speed over time
x = readLines(md_log)
k = grep("Global speed range", x, fixed=TRUE)
y = as.numeric(x[k+1])

# Choose a threshold that will notice early models with artefacts here
# It distinguished cases with/without long-time instabilities that we tested,
# without having to run for much longer.
threshold = 7e-12 # See above before changing this.
if(max(y) < threshold){
    print(c('PASS', max(y)))
}else{
    print(c('FAIL', max(y), run_tag))
}

