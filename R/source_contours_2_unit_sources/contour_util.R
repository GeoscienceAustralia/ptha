#
# Various utilities for working with contours
#
# Gareth Davies, Geoscience Australia 2013+
#


suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(rgdal))

## FIXME: At the moment we use a script to fix a bug in geosphere's antipodal
## I have reported the bug, when it it fixed remove this
#library(geosphere)
source('override_antipodal_geosphere.R')


#' Interface to gdal_contour
#'
#' Use gdal to make a contour (more reliable/efficient than R's algorithms)
#'
#' @param dem = raster for which contours are desired
#' @param contour_levels = elevations for the contours
#' @param out_dsn = Folder to write outputs to
#' @param contour_shp_name = name for the contour shapefile (inside out_dsn)
#' @param delete_tmp_files = TRUE/FALSE -- delete temp files??
#' @param verbose = TRUE/FALSE -- print more or less
#' @param clean_dsn = TRUE/FALSE -- delete files matching contour_shp_name
#'                    in out_dsn. If we dont do this, gdal will throw 
#'                    errors which may not be detected
#' @return A contour, which seems to be of higher quality (and execute faster)
#'       than the native R routines
#'
#' @export
gdal_contour<-function(dem, contour_levels=0.0, 
            out_dsn=NA, contour_shp_name=NA,
            delete_tmp_files=is.na(out_dsn),
            verbose=TRUE, clean_dsn=TRUE){


    # Make directory to write file in
    if(is.na(out_dsn)){
        out_dsn='.GCONTOUR_TMP'
    } 
    dir.create(out_dsn,showWarnings=FALSE)
    
    # Make sure the dem is a tif file 
    if(class(dem)=='RasterLayer'){ 
        # Write the dem to a file
        if(verbose) print('Writing dem to file ...')
        tmp_raster_filename=paste0(out_dsn,'/temp_dem.tif')
        writeRaster(dem,tmp_raster_filename, driver='GTiff',overwrite=T)
    }else{
        if(file.exists(dem)){
            tmp_raster_filename=dem
        }else{
            stop('ERROR: dem should either be a raster layer, or a .tif filename')
        }
    }
    
    # Make the contour layer name
    if(!is.na(contour_shp_name)){
        contour_layer=contour_shp_name
    }else{
        contour_layer='contour'
    }

    # Delete all files in out_dsn called contour_shape_name -- or there will be problems
    if(clean_dsn) file.remove(dir(out_dsn, pattern=contour_layer,full.names=T))


    gdal_contour_command=paste('gdal_contour -a elev -fl', paste(contour_levels,collapse=" "), 
                               '-nln ', contour_layer, tmp_raster_filename, out_dsn)
    if(verbose){
        print('Running gdal_contour ...')
        print(gdal_contour_command)
    } 
    system(gdal_contour_command)

    # Read the file
    out=readOGR(out_dsn,layer=contour_layer,verbose=verbose)
  
    # Clean up 
    if(delete_tmp_files){
        unlink(out_dsn, recursive=TRUE)
    } 

    return(out)
}

#' Return a set of points with spacing='spacing' along a spatialLines object SL
#'
#' The spacing is calculated with linear interpolation (no great circles),
#' so it is assumed that the points defining lines in SL are close enough
#' together that 'great circle' and euclidean interpolation are practically identical.
#' However distances are computed assuming longlat (and spherical earth) if longlat = TRUE.
#'
#' The units of spacing depend on the projection of SL, and the value of longlat
#' If longlat=FALSE, spacing is in units of SL's coordinates (e.g. metres
#'   for UTM projections, degrees for longlat projections)
#' If longlat=TRUE, SL is in longlat degrees, and spacing is in kilometres
#'
#' INPUTS:
#' @param SL = spatialLines object
#' @param spacing = approximate spacing between points
#' @param n = Number of points in output (if spacing = NULL)
#' @param longlat = TRUE/FALSE
#' @param verbose
#' @return SpatialPointsDataFrame along SL, with data identifying the
#'   corresponding line index in SL
#'
#' @export
approxSpatialLines<-function(SL, spacing=NULL, n=NULL, longlat=FALSE,verbose=FALSE){

    if(is.null(n)) stopifnot(!is.null(spacing))
    if(is.null(spacing)) stopifnot(!is.null(n))

    all_points = c()
    for(j in 1:length(SL@lines)){

        if(verbose) print(paste0('Line ', j))

        # Get the segment distance between points defining each line in SL
        SL_seglength = lapply(SL@lines[[j]]@Lines, 
            myfun<-function(x) LineLength(x,longlat, sum=F))

        # Store x,y points on each segment
        newxlist = list()
        newylist = list()
        for(i in 1:length(SL_seglength)){
            seg_coords = SL@lines[[j]]@Lines[[i]]@coords
            seg_lnth = c(0,cumsum(SL_seglength[[i]]))

            if (!is.null(spacing)){
                # Space out_lnth_points evenly around the segment
                out_lnth_num = ceiling(max(seg_lnth)/spacing)

            }else if(!is.null(n)){
                out_lnth_num = n-1
            }

            out_lnth = seq(0, max(seg_lnth), len=out_lnth_num+1) #[1:out_lnth_num]

            # Interpolate the new coordinates
            new_x = approx(seg_lnth, seg_coords[,1], out_lnth, rule=2)$y 
            new_y = approx(seg_lnth, seg_coords[,2], out_lnth, rule=2)$y

            newxlist[[i]] = new_x
            newylist[[i]] = new_y
        
            all_points_local = cbind(unlist(newxlist), unlist(newylist))
            all_points_local = cbind(all_points_local, rep(j, length(all_points_local[,1])))
        
        }
        
        all_points = rbind(all_points,all_points_local)
    }

    # Convert to spatial points
    if(length(all_points) > 0){
        out = SpatialPointsDataFrame(all_points[,1:2], data=data.frame(SLID=all_points[,3]),
            proj4string=CRS(proj4string(SL)))
    }else{
        out = NULL
    }
    return(out)
}

##########################################################################################


#####################################################################################

#' Convert SpatialLinesDataFrame to SpatialPolygonsDataFrame -- so we can use 'over' to
#' identify those not containing coastal points
#'
#' @param sldf = SpatialLinesDataFrame object (e.g. from contouring) where the lines can be viewed as polygons
#' @return A spatial polygons data frame
#' 
#' @export
SpatialLinesDF2Polygons<-function(sldf){

    mypolylist = lapply(sldf@lines, 
        myfun<-function(x){
            mycrds = x@Lines[[1]]@coords
            l = length(mycrds[,1])
            if(any(mycrds[1,1:2] != mycrds[l,1:2])){
                mycrds=rbind(mycrds, mycrds[1,])
            }
            return(Polygons(list(Polygon(mycrds, hole=F)), ID=1))
        }
    )

    for(i in 1:length(mypolylist)){
       mypolylist[[i]]@ID = as.character(i)
    }

    sldf_poly = SpatialPolygons(mypolylist, proj4string=CRS(proj4string(sldf)))
    sldf_poly = SpatialPolygonsDataFrame(sldf_poly, 
                                         data=data.frame(ID=1:length(sldf_poly)))

    return(sldf_poly)
}


##################################################################################
#
# DRAFT code to split spatialLines appropriately
#

library(sp)
rewrap_SpatialLines<-function(SL, newmin=0, split_eps=1e-05, progress=TRUE, ll_TOL=180){
    # Split lines at the min/max of the newrange
    # FIXME: Fails for some values of newmin / overly conservative
    # CAN SUPRESS FAILURE HERE
    store_ll_TOL = get_ll_TOL()
    store_ll_warn = get_ll_warn()

    set_ll_TOL(ll_TOL)
    set_ll_warn(TRUE)
    
    nwSL = SL
    #proj4string(nwSL)<-CRS(as.character(NA))
    #proj4string(nwSL)<-CRS('+proj=longlat')
    nwSL<- nowrapSpatialLines(nwSL, offset=newmin, eps=c(1,1)*split_eps)
    nwSL<- nowrapSpatialLines(nwSL, offset=newmin+360, eps=c(1,1)*split_eps)
    nwSL<- nowrapSpatialLines(nwSL, offset=newmin-360, eps=c(1,1)*split_eps)
    #nwSL=SL
    #nwSL@lines <- lapply(nwSL@lines, .nowrapLines, offset = newmin+360, eps = eps)

    l1 = length(nwSL)
    # Set up progress bar
    if(progress) pb = txtProgressBar(min=0,max=l1,style=3)

    for(j in 1:l1){
        if(progress) setTxtProgressBar(pb,j)
        l2 = length(nwSL@lines[[j]]@Lines)
        for(i in 1:l2){
            if(nwSL@lines[[j]]@Lines[[i]]@coords[1,1] <= newmin){
                # Fix lines (< newmin)
                  
                nwSL@lines[[j]]@Lines[[i]]@coords[,1] = 
                   (nwSL@lines[[j]]@Lines[[i]]@coords[,1]-newmin)%%(360)+newmin

            }else if(nwSL@lines[[j]]@Lines[[i]]@coords[1,1] >= (newmin+360)){
                # Fix lines (> newmin+360)
                nwSL@lines[[j]]@Lines[[i]]@coords[,1] = 
                    (nwSL@lines[[j]]@Lines[[i]]@coords[,1]-newmin)%%(360)+newmin
            }
         }
     }
     if(progress) close(pb)

     output = SpatialLines(nwSL@lines, proj4string=CRS(proj4string(SL)))

     set_ll_TOL(store_ll_TOL)
     set_ll_warn(store_ll_warn)
     return(output)
}

#
# library(mapdata)
# library(maptools)
# m1=map(database='world')
# m2=map2SpatialLines(m1)
# proj4string(m2)=CRS("+proj=longlat +datum=WGS84")
# m3=rewrap_SpatialLines(m2,-30,eps=1.0e-03)

#rewrap_SpatialLines2<-function(SL, newmin=0, split_eps=1e-05,progress=TRUE){
#    #browser()
#    # Split lines at the min/max of the newrange
#    # FIXME: Fails for some values of newmin / overly conservative
#    # HOWEVER, THIS REPRESENTS AN ALTERNATIVE APPROACH TO PREVIOUS
#    nwSL=SL
#    #proj4string(nwSL)<-CRS(as.character(NA))
#    #proj4string(nwSL)<-CRS('+proj=longlat')
#    #nwSL=SL
#    #nwSL@lines <- lapply(nwSL@lines, .nowrapLines, offset = newmin+360, eps = eps)
#
#    l1=length(nwSL)
#    # Set up progress bar
#    if(progress) pb=txtProgressBar(min=0,max=l1,style=3)
#
#    for(j in 1:l1){
#        if(progress) setTxtProgressBar(pb,j)
#        l2=length(nwSL@lines[[j]]@Lines)
#        for(i in 1:l2){
#                nwSL@lines[[j]]@Lines[[i]]@coords[,1] = 
#                   (nwSL@lines[[j]]@Lines[[i]]@coords[,1]-newmin)%%(360)+newmin
#
#         }
#     }
#   
#     nwSL<-SpatialLines(nwSL@lines,proj4string=CRS(proj4string(SL))) 
#     nwSL<- nowrapSpatialLines(nwSL, offset=newmin, eps=c(1,1)*split_eps)
#     nwSL<- nowrapSpatialLines(nwSL, offset=newmin+360, eps=c(1,1)*split_eps)
#
#     if(progress) close(pb)
#
#     output=SpatialLines(nwSL@lines, proj4string=CRS(proj4string(SL)))
#     return(output)
#}
#

#}


