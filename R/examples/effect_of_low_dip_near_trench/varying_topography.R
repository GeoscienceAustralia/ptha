#
# Tests with approaches like in Zhu et al 2025, methods 'B' and 'C'
#
library(rptha)

width = 50
length = 50
strk = 0

dip = 5
topo_angle = 5 # Slope of topography above horizontal where x > 0

# Methods B and C vary in terms of their 'equivalent' fault downdip width
# changes and centroid depth. But the horizontal width stays the same
horizontal_width = width/cos(dip/180*pi)

# centroid horizontal position in units of meters
centroid_x = horizontal_width / 2 * 1000
centroid_y = 0

xs = seq(-100.5,100.5,len=100)
ys = seq(-100.5,100.5,len=100)
# Output solution at mypts (units of m)
mypts=expand.grid(xs, ys)*1000


#
# Their method B
#
method_b_width = width
method_b_centroid_depth = method_b_width/2*sin(dip/180*pi) # Rupture to trench
method_b_dip = dip
method_b = okada_tsunami(centroid_x, centroid_y, method_b_centroid_depth, 
    strk, method_b_dip, length, method_b_width, 0, 1, mypts[,1],mypts[,2], dstmx_min = 100)
method_b_dz = method_b$zdsp # Vertical component
dim(method_b_dz) = c(length(xs), length(ys))

# Say the slope begins at the trench
method_b_horiz = -1.0 * tan(topo_angle/180*pi) * (mypts[,1] > 0) * method_b$edsp # Due to sloping topography
dim(method_b_horiz) = c(length(xs), length(ys))

method_b_combined = method_b_dz + method_b_horiz # Vertical + sloping topography

#
# Their method C
#
method_c_width = horizontal_width/cos((dip+topo_angle)/180*pi)
method_c_centroid_depth = method_c_width/2*sin((dip+topo_angle)/180*pi) # Rupture to trench
method_c_dip = dip + topo_angle
method_c = okada_tsunami(centroid_x, centroid_y, method_c_centroid_depth, 
    strk, method_c_dip, length, method_c_width, 0, 1, mypts[,1],mypts[,2], dstmx_min = 100)
method_c_dz = method_c[[3]]
dim(method_c_dz) = c(length(xs), length(ys))

j = 50
YLIM = c(-1,1)*max(abs(range(c(range(method_c_dz), range(method_b_combined), range(method_b_dz)))))
plot(xs, method_c_dz[,j], t='l', ylim=YLIM)
points(xs, method_b_combined[,j], t='l', col='red')
points(xs, method_b_dz[,j], t='l', col='darkgreen')
legend('bottomleft', c('C', 'B + horizontal', 'B'), lty='solid', col=c('black', 'red', 'darkgreen'), bty='n', cex=1.5)
title(main='The Okada computation of Method C involves a longer fault width to cover the same horizontal distance. \n So it is a bit like having a larger earthquake, although would not be accounted that way (just a trick to deal with the slope).')

