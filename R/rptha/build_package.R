library(devtools)
document('.')
document('.')

# Modify the 'version' to have the git revision number
version_extra = system('git describe --abbrev=50 --always --tags --dirty', intern=TRUE)
description_lines = readLines('DESCRIPTION')
dvi = grep('Version:', description_lines)
description_lines[dvi] = paste0('Version: 0.0.', version_extra)
cat(description_lines, file='DESCRIPTION', sep="\n")
system('git commit -a -m "autobuild"')

build('.')
