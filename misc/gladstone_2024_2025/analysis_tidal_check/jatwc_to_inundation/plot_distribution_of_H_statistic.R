# Look at the distributions of H95 statistics in each ATWS Coastal Zone
files = Sys.glob('Inundation_zones/*/all-JATWC*.RDS')
all_H = lapply(files, readRDS)
names(all_H) = basename(dirname(files))

# Convert to a single data.frame for each zone
all_H_df = lapply(all_H, function(x) do.call(rbind, x))
names(all_H_df) = names(all_H)

png('all_H_distributions_in_coastal_zones.png', width=15, height=5, units='in', res=200)
par(mfrow=c(2,5))
par(mar=c(2,2.5,2,2))
par(oma=c(4, 0, 0, 0)) # Space for legend at bottom
for(i in 1:length(all_H_df)){
    H_df = all_H_df[[i]]
    no_threat = mean(H_df[,2] < 0.2)
    marine_warning = mean(H_df[,2] >= 0.2 & H_df[,2] < 0.55)
    land_warning_minor = mean(H_df[,2] >= 0.55 & H_df[,2] < 1.5)
    land_warning_major = mean(H_df[,2] >= 1.5)

    barplot(c(no_threat, marine_warning, land_warning_minor, land_warning_major),
        #names.arg = c('no_threat', 'marine_warning', 'land_warning_minor', 'land_warning_major'),
        #las = 2,
        col = c('green', 'blue', 'red', 'black'),
        ylim=c(0, 0.6),
        cex.names=1.5,
        cex.axis=1.4)
    title(main=names(all_H_df)[i], cex.main=2)
}
# Put a legend at the bottom
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(2, 2, 2, 2), new=TRUE, xpd=NA)
legend(-0.60, 0.017, 
    c('No threat (< 0.2)', 'Marine W. (0.2 - 0.55)', 'Minor Land W. (0.55 - 1.5)', 'Major Land W. (>=1.5)'),
    fill=c('green', 'blue', 'red', 'black'), horiz=TRUE, bty='n', cex=2.5)
dev.off()
