#
# Compute flood hazard categories from Flood risk management guideline FB03
#

library(terra)

#' Hazard categories as per Figure 1 of Flood risk management guideline FB03
#'
#' All inputs are terra:rast objects.
#'
#' @param vd The flux (velocity x depth, m^2 / s)
#' @param depth The depth (m)
#' @param speed The speed (m/s)
#' @return hazard category as a terra:rast
compute_flood_hazard_categories<-function(vd, depth, speed){

    if(!(compareGeom(vd, depth) & compareGeom(vd, speed))){
        stop('Incompatible raster extents')
    }

    hazard_category = depth * 0

    # Firstly compute hazard categories based purely on VD.
    # Later we update based on thresholds for depth and speed
    hazard_category[vd > 0   & vd <= 0.3] = 1
    hazard_category[vd > 0.3 & vd <= 0.6] = 2
    # Category 3 has the same flux limit as category 2, but different depth/speed limits
    hazard_category[vd > 0.6 & vd <= 1.0] = 4
    hazard_category[vd > 1.0 & vd <= 4.0] = 5
    hazard_category[vd > 4.0] = 6

    # Now apply limits based on velocity
    sl = speed > 4.0
    hazard_category[sl] = 6
    sl = speed > 2.0
    hazard_category[sl] = pmax(hazard_category[sl], 5)

    # Now apply limits based on depth
    dl = depth > 4.0
    hazard_category[dl] = 6
    dl = depth > 2.0
    hazard_category[dl] = pmax(hazard_category[dl], 5)
    dl = depth > 1.2
    hazard_category[dl] = pmax(hazard_category[dl], 4)
    dl = depth > 0.5
    hazard_category[dl] = pmax(hazard_category[dl], 3)
    dl = depth > 0.3
    hazard_category[dl] = pmax(hazard_category[dl], 2)

    hazard_category[hazard_category==0] = NA

    return(hazard_category)
}

test_compute_flood_hazard_categories_with_plot<-function(){
    # This makes a plot depicting how the hazard categories are related to speed and depth
    # It is created by making rast objects with depth, speed, and vd, then doing the calculation.
    #
    # It is easy to visually compare the plot with Figure 1 from Flood risk management guideline FB03,
    # to see that everything is OK.

    dx = 0.01
    depths = seq(0, 6, by=dx)
    speeds = seq(0, 5, by=dx)

    rast_template = rast(nrows=length(depths), ncols=length(speeds), 
        xmin=0-dx/2, xmax=max(speeds)+dx/2, ymin=0-dx/2, ymax=max(depths)+dx/2)
   
    # Make a raster with depths varying vertically. Need to rev(depths) to meet
    # the rast() conventions
    depthmat = matrix(rev(depths), nrow=length(depths), ncol=length(speeds), byrow=FALSE)
    depthrast = rast_template
    values(depthrast) = depthmat

    # Make a raster with speeds varying horizontally
    speedmat = matrix(speeds, nrow=length(depths), ncol=length(speeds), byrow=TRUE)
    speedrast = rast_template
    values(speedrast) = speedmat

    # Flux
    vdrast = depthrast * speedrast

    # This should lead to a raster matching the Figure 1 of "Flood risk
    # management guideline FB03"
    hc = compute_flood_hazard_categories(vdrast, depthrast, speedrast)    

    # Make a plot
    plot(range(speeds), range(depths), asp=0.5, col='white', xlab='Speed (m/s)', 
        ylab='Depth (m)', cex.axis=1.5, cex.lab=1.5, 
        main='Flood hazard categories', cex.main=2, 
        xaxs='i', yaxs='i')
    cols = c('blue', 'lightblue', 'turquoise', 'green', 'yellow', 'red')
    plot(hc, add=TRUE, col=cols, alpha=0.5, legend=FALSE)

    # Add labels
    text_xs = rep(0.3, 6)
    text_ys = c(0.15, 0.4, 0.8, 1.5, 3.0, 4.5)
    text_vals = paste0('H', seq(1, 6))
    text(text_xs, text_ys, text_vals, cex=1.5)

    legend('topright', fill=rev(cols), legend=c(
        'H6: unsafe for vehicles and people. \n All building types considered vulnerable to failure',
        'H5: unsafe for vehicles and people. Buildings require \n special engineering design and construction',
        'H4: unsafe for vehicles and people',
        'H3: unsafe for vehicles, children and the elderly',
        'H2: unsafe for small vehicles',
        'H1: generally safe for people, vehicles and buildings'),
        title = 'NSW Flood risk management guideline FB03, 2023',
        title.cex=1.4, title.font=4,
        y.intersp = c(1.5,2,1,1,1,1),
        bty='n', cex=1.3, text.font=3)
        

}

make_flood_hazard_categories_png_plot<-function(){
    png('Flood_hazard_categories_FB03.png', width=10.62, height=7.49, units='in', res=200)
    test_compute_flood_hazard_categories_with_plot()
    dev.off()
}

