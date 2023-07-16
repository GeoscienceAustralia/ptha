#
# Find gauges that are 'near' to one another
# By plotting nearby gauges we can compare the results and sometimes see differences due to instrumentation. 
#
source('global_variables.R')
library(sp)
DISTANCE_THRESHOLD = 0.35 # km

tg_data = read.csv(TIDEGAUGE_METADATA_TABLE_FILE)
coordinates(tg_data) = c('lon', 'lat')
proj4string(tg_data) = CRS("+init=epsg:4326")

# Find gauges close to each other
nearby_inds = zerodist(tg_data, zero=DISTANCE_THRESHOLD)

print(paste0('Sites within ', DISTANCE_THRESHOLD, 'km of each other'))
options(width=120)
cbind(tg_data$name[nearby_inds[,1]], tg_data$name[nearby_inds[,2]])

# Plot nearby sites
pdf(paste0(OUTPUT_GRAPHICS_DIR, '/nearby_tide_gauges_comparison.pdf'), width=12, height=7)
par(mfrow=c(2,1))
for(i in 1:nrow(nearby_inds)){

    i1 = nearby_inds[i,1]
    i2 = nearby_inds[i,2]

    file_1 = tg_data$postprocessed_file[i1]
    file_2 = tg_data$postprocessed_file[i2]

    site_1 = read.csv(file_1, comment.char='#')
    site_2 = read.csv(file_2, comment.char='#')

    site_names = paste0(tg_data$name[i1], ' / ', tg_data$name[i2])

    if(!all(is.na(site_1$stage))){
    plot(site_1$juliant, site_1$stage - median(site_1$stage, na.rm=TRUE), t='l', main=site_names, 
        xlab='Julian Day', ylab='Stage (m above median)')
    points(site_2$juliant, site_2$stage - median(site_2$stage, na.rm=TRUE), t='l', col='red')
    abline(v=julian(explosion_start_time), lty='dashed', col='blue', lwd=2)
    grid()
    }else{
        plot(c(0, 1), c(0, 1), main=paste0(site_names, '\n MISSING STAGE DATA (licence limitations)'))
    }

    YLIM = c(-1,1)*max(c(abs(site_1$resid3h), abs(site_2$resid3h)), na.rm=TRUE) + c(-0.1, 0.1)
    plot(site_1$juliant, site_1$resid3h + 0.1, xlim=c(19007, 19009), t='l', main='Residual near HT explosion (blue dashed line)', 
        xlab='Julian Day', ylab='Residual (3h or shorter)', ylim=YLIM)
    points(site_2$juliant, site_2$resid3h - 0.1, t='l', col='red')
    legend('topright', c(tg_data$name[i1], tg_data$name[i2]), lty=c(1, 1), col=c('black', 'red'), pch=c(NA, NA), horiz=TRUE, bty='n')
    abline(v=julian(explosion_start_time), lty='dashed', col='blue', lwd=2)
    grid()
}
dev.off()
