#
# Put together gauges
#
library(rptha)
ptha18 = new.env()
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
file_nci  = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

older_gauges = ptha18$get_all_gauges()
# Drop the 'elevation' from older_gauges
older_gauges = older_gauges[,c(1,2,4)]

# Only keep gauges in Indian Ocean / WA region
global_ll = c(16.0  , -75.0)
global_ur = c(130.0 ,  33.0)
keep = which(older_gauges$lon >= global_ll[1] & older_gauges$lon <= global_ur[1] &
             older_gauges$lat >= global_ll[2] & older_gauges$lat <= global_ur[2])
older_gauges = older_gauges[keep,]

# Add in some extras
extra_gauges = coordinates(readOGR('gauges_extra', layer='gauges_extra'))

# From inspecting the gaugeID values for teh older gauges, I see the following numbers will not overlap
new_gauge_ID = 80000 + 1:nrow(extra_gauges) + 0.6
extra_gauges = cbind(extra_gauges, new_gauge_ID)
extra_gauges = as.data.frame(extra_gauges)
names(extra_gauges) = names(older_gauges)

final_gauges = rbind(older_gauges, extra_gauges)

write.csv(final_gauges, file='point_gauges_2023_11_30.csv', row.names=FALSE)

