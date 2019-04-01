source('../../plot.R')

#
# Read the model log and see if the test passed
# (just a mass conservation check)
#
model_log = readLines('outfile.log')
if(model_log[length(model_log)] == ' PASS'){
    print('PASS')
}else{
    print('FAIL')
}

lx = 320
ly = 415

x = get_all_recent_results()
nx = dim(x$stage)[1]
ny = dim(x$stage)[2]
ts = length(x$time)


depth = x$stage[,,ts] - x$elev0
elev = x$elev0
xs = matrix(1:nx, nrow=nx, ncol=ny)*(lx/nx)
ys = matrix(1:ny, nrow=nx, ncol=ny, byrow=TRUE)*(ly/ny)

pdf('merewether_plots.pdf', width=10, height=10)
# Close up arrows plot
image(xs[,1], ys[1,], elev, asp=1, col=rainbow(255), xlim=c(100, 250), ylim=c(100, 250),
    main=paste0('Velocity up close: Time = ', round(x$time[ts],3)))
arrows(xs, ys, xs+x$ud[,,ts]/depth, ys+x$vd[,,ts]/depth, length=0.001)

# Froude plot
vel = sqrt(x$ud[,,ts]**2 + x$vd[,,ts]**2)/(depth + 1.0e-10)
froude = vel/sqrt(9.81 * depth + 1.0e-10)

image(xs[,1], ys[1,], (froude > 1) + (froude > 0) , asp=1, main=paste0('Froude number: Time = ', round(x$time[ts],3)))


hazard = vel * depth
hazard_cat = (hazard > 0) + (hazard > 0.2) + (hazard > 0.5) + (hazard > 1) + (hazard > 2)
image(xs[,1], ys[1,],  hazard_cat, col=c('white', 'skyblue', 'lightgreen', 'yellow', 'orange', 'red'), asp=1,
    main=paste0('Hazard: Time =', round(x$time[ts],3)))
dev.off()
