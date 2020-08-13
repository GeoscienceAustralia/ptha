if(FALSE){
    # Extract data from the models using the older post-processing that had more gauge locations
    # Explanation:
    #    Originally I made this plot with post-processed files that included DARTs and other points.
    #    Later the outputs were been stripped down to only that used in the nearshore testing paper,
    #      and an inconvenience is that the DART model-record required for this plot is missing.
    #    For simplicity, I used this code to write the required model data using the files that
    #      also contain the DART (the SWALS model is unchanged -- it is the post-processed outputs that have been reduced)

    x_lf = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_linear_friction-0-highres_NSW/RUN*/gauge*.RDS'))
    x_nf = readRDS(Sys.glob('../gauge_RDS_files/OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_no_friction-0-highres_NSW/RUN*/gauge*.RDS'))

    # Find index of DART_52406, which we will plot
    ind = grep('DART_52406', names(x_nf))
    ind2 = grep('DART_52406', names(x_lf))
    if(ind != ind2) stop('error -- models do not have same gauges')

    file1 = 'Model_output_52406_linear_with_no_friction.csv'
    write.csv(data.frame(time = x_nf[[ind]]$model_time, model_stage = x_nf[[ind]]$model_stage), file=file1, row.names=FALSE)
    file2 = 'Model_output_52406_linear_with_linear_friction.csv'
    write.csv(data.frame(time = x_lf[[ind]]$model_time, model_stage = x_lf[[ind]]$model_stage), file=file2, row.names=FALSE)

}

# Read the stored files (simpler than hacking the post-processing script to include the DART in this case only)
x_nf = read.csv('Model_output_52406_linear_with_no_friction.csv', colClasses='numeric')
x_lf = read.csv('Model_output_52406_linear_with_linear_friction.csv', colClasses='numeric')


png('Compare_linear_friction_with_exponential_transform.png', width=10, height=9, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(5.1, 5.1, 4.1, 2.1))

options(scipen=0)
plot(x_nf$time/3600, x_nf$model_stage, t='l', xlab='Time (hours post-earthquake)', ylab='Stage (m)',
     main='DART 52406: LSWE, and LSWE + linear-friction', 
     cex.main=2, cex.axis=1.6, cex.lab=2, ylim=c(-0.15, 0.2))
points(x_lf$time/3600, x_lf$model_stage, t='l', col='red')
grid(col='orange'); abline(h=0, col='orange')
legend('topleft', c('LSWE', 'LSWE + linear-friction'), col=c('black', 'red'), 
       lwd=c(1, 1), pch=c(NA, NA), cex=2, horiz=TRUE, bty='n')

approximate_linear = x_nf$model_stage * exp(-1e-05/2 * x_nf$time)

plot(x_nf$time/3600, approximate_linear, t='l', xlab='Time (hours post-earthquake)', ylab='Stage (m)',
     main='As above, applying the approximate linear-friction transformation to the LSWE', 
     cex.main=2, cex.axis=1.6, cex.lab=2, ylim=c(-0.15, 0.2))
points(x_lf$time/3600, x_lf$model_stage, t='l', col='red')
grid(col='orange'); abline(h=0, col='orange')
legend('topleft', c('Approximate-linear-friction-solution', 'LSWE + linear-friction'), col=c('black', 'red'), 
       lwd=c(1, 1), pch=c(NA, NA), cex=2, horiz=TRUE, bty='n', xjust=1)

options(scipen=5)
plot(x_nf$time/3600, approximate_linear - x_lf$model_stage, t='l',
     xlab='Time (hours post-earthquake)', ylab='Error (m)', main='Error in the approximate linear-friction solution (y-axis range of 0.001 m)',
     cex.main=2, cex.axis=1.6, cex.lab=2, ylim=c(-1,1)*1e-03)
grid(col='orange'); abline(h=0, col='orange')
dev.off()
