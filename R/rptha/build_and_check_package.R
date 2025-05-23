#
# Run the test suite and build
#
# This currently doesn't work on Windows (as it requires unix shell commands). But the tests can
# be run in Windows once everything is installed by starting R in the 'tests' directory and running:
#
#     source('testthat.R')
#


library(devtools)

# Modify the 'version'
try({
    description_lines = readLines('DESCRIPTION')
    dvi = grep('Version:', description_lines)
    new_version_number = as.numeric(strsplit(description_lines[dvi], '\\.')[[1]][3]) + 1
    description_lines[dvi] = paste0('Version: 0.1.', new_version_number)
    cat(description_lines, file='DESCRIPTION', sep="\n")
    system('git commit -a -m "......increment version number"')
})

document('.')
document('.')
system('mkdir -p ../old_builds')
system('mv ../rptha_0.0*.tar.gz ../old_builds')
build('.')
setwd('..')
system('R CMD check rptha_0.1*.tar.gz')
system('mv rptha_0.1*.tar.gz old_builds')
