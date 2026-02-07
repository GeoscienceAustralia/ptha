#
# Check the funwave results
#
library(fields)

## Refined version of the original test problem, with dispersion switched off.
input_negative_elev_file = '../input/depth.txt_refined_2.txt'
input_max_stage_file = './results_higher_res/hmax_00029'
nx = 702*2
ny = 602*2

# Initial elevation (need to flip sign to make ocean negative)
elev_start = -matrix(scan(input_negative_elev_file), ncol=ny)

# Max stage 
hmax = matrix(scan(input_max_stage_file), ncol=ny)

# Get max island runup
is_dry = (hmax == 0)
is_near_dry = is_dry * 0
# Max-runup zone == wet-cells with at least 1 dry neighbour
is_near_dry[2:(nx-1),2:(ny-1)] = (1 - is_dry[2:(nx-1), 2:(ny-1)])*pmax(
    is_dry[3:nx    , 2:(ny-1)], 
    is_dry[1:(nx-2), 2:(ny-1)], 
    is_dry[2:(nx-1), 3:ny    ], 
    is_dry[2:(nx-1), 1:(ny-2)]) 
wd_peak_inds = which(is_near_dry == 1, arr.ind=TRUE)
wd_peak_stage = apply(wd_peak_inds, 1, f<-function(x) hmax[x[1], x[2]])

# Convert to a radial coordinate
island_centre_ind = apply(which(elev_start == max(elev_start), arr.ind=TRUE), 2, mean)
wd_theta = atan((wd_peak_inds[,2] - island_centre_ind[2])/(wd_peak_inds[,1] - island_centre_ind[1])) -
    pi*(wd_peak_inds[,1] - island_centre_ind[1] < 0)
wd_theta = (wd_theta/pi*180 + 90)%%360
plot(wd_theta, wd_peak_stage)
output=data.frame(theta=wd_theta, max_runup=wd_peak_stage)
write.csv(output, file='max_island_runup.csv', row.names=FALSE)


#
# Check distance between initial wave peak and the island (useful reference point for setting up other models)
#
eta0 = matrix(scan('results_higher_res/eta_00000'), ncol=ny)
hny = ceiling(ny/2)
wave_peak = which.max(eta0[,hny]*(eta0[,hny] > elev_start[,hny]+0.001))
cells_wavepeak2island = (island_centre_ind[1] - wave_peak)
distance_wavepeak2island = cells_wavepeak2island * 0.025
