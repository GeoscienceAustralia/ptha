source('../../plot.R')
md_dir = sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1]
md = get_multidomain(md_dir)
x = md[[1]]
nx = dim(x$stage)[1]
ny = dim(x$stage)[2]
ts = length(x$time)

#
# Check mass conservation -- this model should have an 'unexplained change'
# due to the discharge source (which SWALS doesn't track). In the model code
# we wrote the volume to the file
#
model_log = readLines(Sys.glob(paste0(md_dir, '/multi*.log'))[1])
# Get 'unexplained change'
unexplained_change = as.numeric(strsplit(model_log[max(grep('unexplained change', model_log))], ':')[[1]][2])
# Get the volume we added
expected_change = as.numeric(model_log[grep('Expected mass change due to inflows', model_log)+1])

if(abs(expected_change - unexplained_change) < 1.0e-06*abs(unexplained_change)){
    print('PASS')
}else{
    print('FAIL')
}

# Length of x/y domain extents
lx = 320
ly = 415

depth = x$stage[,,ts] - x$elev0
elev = x$elev0

# x,y coordinate, translated with lower-left = 0,0
xs = (matrix(1:nx, nrow=nx, ncol=ny)-0.5)*(lx/nx)
ys = (matrix(1:ny, nrow=nx, ncol=ny, byrow=TRUE)-0.5)*(ly/ny)

#pdf('merewether_plots.pdf', width=10, height=10)
png('Velocity_vector_plot.png', width=10, height=10, units='in', res=300)
# Close up arrows plot
image(xs[,1], ys[1,], elev, asp=1, col=rainbow(255), xlim=c(100, 250), ylim=c(100, 250),
    main=paste0('Velocity up close: Time = ', round(x$time[ts],3)))
arrows(xs, ys, xs+x$ud[,,ts]/depth, ys+x$vd[,,ts]/depth, length=0.001)
dev.off()

# Froude plot
png('Froude_number_plot.png', width=10, height=10, units='in', res=300)
vel = sqrt(x$ud[,,ts]**2 + x$vd[,,ts]**2)/(depth + 1.0e-10)*(depth > 1.0e-05)
froude = vel/sqrt(9.81 * depth + 1.0e-10)*(depth > 1.0e-05)

image(xs[,1], ys[1,], (froude > 1) + (froude > 1.0e-05) , asp=1, 
      main=paste0('Froude number: Time = ', round(x$time[ts],3)))
dev.off()

png('Hazard_plot.png', width=10, height=10, units='in', res=300)
hazard = vel * depth
hazard_cat = (hazard > 1.0e-3) + (hazard > 0.2) + (hazard > 0.5) + (hazard > 1) + (hazard > 2)
image(xs[,1], ys[1,],  hazard_cat, col=c('white', 'skyblue', 'lightgreen', 'yellow', 'orange', 'red'), asp=1,
    main=paste0('Hazard: Time =', round(x$time[ts],3)))
dev.off()

# Compare observed peak flow depths with model, and reports of other models.
#
# I obtained the observed and modelled values from:
#   Smith, G. & Wasko, C. Revision Project 15: Two Dimensional Simulations In
#   Urban Areas - Representation of Buildings in 2D Numerical Flood Models
# Note the above report doesn't contain coordinates -- but the coordinates for
# some points are contained in the ANUGA test-case, and those are the ones
# tested here.
#
observed = read.csv('point_observations.csv')
tmp = make_max_stage_raster(x, proj4string='', return_elevation=TRUE)
max_stage = tmp$max_stage
max_stage[tmp$max_stage < tmp$elevation + 1.0e-03] = NA

# The model is 'just dry' at some observed points -- the observations may well be
# inundation limits, so this is not surprising. To get observations, slightly buffer the
# point locations (by 1.5 m) and take the mean of the non-missing values 
model_stage = extract(max_stage, cbind(observed$x, observed$y), buffer=1.5)
model_stage = unlist(lapply(model_stage, f<-function(x) mean(x, na.rm=TRUE)))

# The results are comparable to TUFLOW at various grid resolutions (also in 'point_observations.csv')

model_err = mean(abs(model_stage - observed$Observed))
tuflow1m_err = mean(abs(observed$TUFLOW_1m_Table9 - observed$Observed))
tuflow2m_err = mean(abs(observed$TUFLOW_2m_Table9 - observed$Observed))
tuflow0p5m_err = mean(abs(observed$TUFLOW_0.5m_Table9 - observed$Observed))

# Ad-hoc error criteria. Note that the TUFLOW 1m was probably the 'calibrated'
# version, and it does significantly better than the 0.5 m, and slightly better than
# 2m. 
if(model_err < max(tuflow1m_err, tuflow2m_err, tuflow0p5m_err)){
    print('PASS')
}else{
    print('FAIL')
}
