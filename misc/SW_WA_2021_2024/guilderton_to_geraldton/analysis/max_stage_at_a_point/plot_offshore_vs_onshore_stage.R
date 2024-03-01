# Offshore
x = readRDS('rdata/run_stage_at_target_point_114.3333_-31.3333.RDS')

# Hillarys
y = readRDS('rdata/run_stage_at_target_point_115.7393_-31.8225.RDS') 

png('results/Offshore_max_stage_vs_nearshore_max_stage.png', width=6, height=4, units='in', res=300)
plot(x$max_stage-0.6, y$max_stage-0.6, 
     xlab='Offshore tsunami maxima (m)', ylab='Nearshore tsunami maxima (m)',
     cex.axis=1.5, cex.lab=1.5, pch=19, xlim=c(0,1), ylim=c(0, 3))
grid()
dev.off()
