##
## Various utility functions for working with points
##

haz_pts_2_ursga_format<-function(haz_pts, outfile='haz_pts_gebco08_trimmed.txt'){
    # Write points into the ursga text format
    # INPUT: haz_pts= hazard points in SpatialPointsDataFrame format
    #        outfile = output file name
    # OUTPUT: None, but writes the output file we need

    haz_pt_keep=1:length(haz_pts)

    # Write to output file for URSGA
    cat(length(haz_pt_keep), 
        file=outfile,fill=1)

    outdata=cbind(coordinates(haz_pts[haz_pt_keep,])[,2:1], 
                  1:length(haz_pt_keep))
    options(digits=12)
    for(i in 1:length(outdata[,1])){
        cat(paste(c(outdata[i,],'\n'), collapse=" "), 
            file=outfile, append=T)
    }

    return()
}

###################################################################

cut_hazpts_in_poly<-function(haz_orig, 
                             haz_cull_poly_name, 
                             new_haz_pts_name, 
                             lower_left=-65,
                             outdir='OUTPUTS/',
                             indir='./'){
    ### Code to remove hazard points from region defined by haz_cull_poly_name
    ### Also, the points can have the lower-left adjusted to lower_left
    ###
    ### INPUT: haz_orig: SpatialPointsDataFrame of the hazard points
    ###        haz_cull_poly_name: Directory/Layer name for the polygon where we cut points
    ###        new_haz_pts_name : Name for output file with adjusted hazard points
    ###        lower_left: Translate the points so this is the lower left
    ###        outdir: Parent directory for output file
    ###        indir: Parent directory for haz_cull_poly_name
    ###
    ### OUTPUT: new hazard points, + various important side effects (IO)
    ###
    library(sp)
    library(rgdal)
    set_ll_TOL(5000.0) # Don't worry if long < -180
    set_ll_warn(TRUE) # Avoid errors related to long < -180

    # Get hazard removal region
    haz_cull=readOGR(paste0(indir,haz_cull_poly_name),
                     layer=haz_cull_poly_name)

    # Cut from polygon
    cutpts=over(haz_orig,haz_cull)
    haz_keep=which(is.na(cutpts[[1]]))

    # Shift lower left
    xlong=coordinates(haz_orig)[,1]
    ylong=coordinates(haz_orig)[,2]
    xlong=xlong*(xlong>=lower_left) + (360+xlong)*(xlong<lower_left)

    haz_new2=SpatialPointsDataFrame(cbind(xlong,ylong)[haz_keep,], 
                              data=haz_orig@data[haz_keep,],
                              proj4string=CRS(proj4string(haz_orig)))

    writeOGR(haz_new2,dsn=paste0(outdir, new_haz_pts_name),
            layer=new_haz_pts_name,driver='ESRI Shapefile',overwrite=TRUE)

    return(haz_new2)
}
