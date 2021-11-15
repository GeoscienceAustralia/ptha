source('../../plot.R')

# We ran the NS-aligned channel case first. The tests below focus on that, but later
# we will check the EW-aligned results.
NS_model = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[2])
EW_model = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])

# Solution from model
ind_y = round(4/16*length(NS_model[[1]]$ys)) # Index that is most of the way downstream
ind_t = 67
model_x = NS_model[[1]]$xs
model_vh = NS_model[[1]]$vd[,ind_y,ind_t]
model_depth = NS_model[[1]]$stage[,ind_y, ind_t] - NS_model[[1]]$elev[,ind_y,ind_t]

# Get the equivalent values from the transposed model 
model_transpose_x = EW_model[[1]]$ys
model_transpose_vh = EW_model[[1]]$ud[ind_y, , ind_t]
model_transpose_depth = EW_model[[1]]$stage[ind_y, , ind_t] - EW_model[[1]]$elev[ind_y, , ind_t]

# Check for consistent results in transposed/regular model
if(all(abs(model_transpose_vh - model_vh) <= 1.0e-05*abs(model_vh))){
    print('PASS')
}else{
    print('FAIL')
}



# Solution from theory
source('shiono_knight_model.R')
theory_x = solution_coarse_matching_discharge$y
theory_vh = - solution_coarse_matching_discharge$U * solution_coarse_matching_discharge$depth
theory_depth = solution_coarse_matching_discharge$depth

# Solution from theory without momentum exchange
theoryB_x = solution_coarse_nolambda$y
theoryB_vh = - solution_coarse_nolambda$U * solution_coarse_nolambda$depth
theoryB_depth = solution_coarse_nolambda$depth

# Make a plot
png('Model_vs_theory_shionoknight.png', width=8, height=4, units='in', res=300)
plot(model_x, model_vh/model_depth, t='o', pch=19, cex=0.3,
     main=' Downstream velocity with eddy-viscosity \n Model vs analytical solution',
     xlab='Cross-channel coordinate (m)', 
     ylab = 'Velocity (m/s)', cex.main=1.4, cex.lab=1.4)

points(theory_x, theory_vh/theory_depth, t='l', col='red')

points(theoryB_x, theoryB_vh/theoryB_depth, t='l', col='darkgreen', 
       lty='dashed')

legend('top', 
       c('Model (with eddy viscosity)', 'Theory (with eddy viscosity)', 'Theory (no eddy viscosity)'), 
       lty=c('solid', 'solid', 'dashed'), pch=c(19,NA, NA), pt.cex=c(0.3, NA, NA),
       col=c('black', 'red', 'darkgreen'), bty='n')
dev.off()

#
# Report a pass/fail statistic
#
model_vel = model_vh/(model_depth + 1.0e-20)
N = length(model_vel)
# Note the model has a 'dry wall' boundary condition on the sides -- so we only use the interior points
theory_at_model_points = approx(theory_x, theory_vh/theory_depth, xout=model_x[2:(N-1)])$y

error_stat = mean(abs(model_vel[2:(N-1)] - theory_at_model_points)/abs(theory_at_model_points))

if(error_stat < 0.02){
    print('PASS')
}else{
    print('FAIL')
}

