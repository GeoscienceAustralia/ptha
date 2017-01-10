library(devtools)

# Modify the 'version'
try({
    description_lines = readLines('DESCRIPTION')
    dvi = grep('Version:', description_lines)
    new_version_number = as.numeric(strsplit(description_lines[dvi], '\\.')[[1]][3]) + 1
    description_lines[dvi] = paste0('Version: 0.0.', new_version_number)
    cat(description_lines, file='DESCRIPTION', sep="\n")
    system('git commit -a -m "    increment version number"')
})

document('.')
document('.')
build('.')
setwd('..')
system('R CMD check rptha_0.0.tar.gz')
