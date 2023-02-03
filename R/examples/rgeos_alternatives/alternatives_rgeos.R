#
# Limited replacements for rgeos functionality
#
# The rptha package historically used rgeos for many geometry operations,
# but in 2023 this is being retired (with sf providing comparable functionality).
# Here we try to make "rgeos-like" interfaces that cover the functionality needed
# in rptha.
#

#' limited replacement for rgeos::gBuffer using sf functionality
#'
#' @export
gBuffer<-function(spgeom, width, quadsegs=5, byid=FALSE){
    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    # st_buffer seems to be always applied on a per-id basis
    newgeom = suppressMessages(st_buffer(geom, dist=width, quadsegs=quadsegs))
    # If we didn't want byID results, then merge
    if(!byid){
        newgeom = suppressMessages(st_union(newgeom))
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}


#' limited replacement for rgeos::gCentroid using sf functionality
#'
#' @export
gCentroid<-function(spgeom, byid=FALSE){
    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    if(byid){
        newgeom = suppressMessages(st_centroid(geom))
    }else{
       #https://gis.stackexchange.com/questions/451041/equivalent-of-gcentroidx-byid-false-in-sf-system 
        newgeom = suppressMessages(st_union(geom))
        newgeom = suppressMessages(st_centroid(newgeom))
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}

#' limited replacement for rgeos::gDistance using sf functionality
#'
#' @export
gDistance<-function(spgeom1, spgeom2=NULL, byid=FALSE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }
    st_crs(geom1) = NA
    st_crs(geom2) = NA

    # st_distance always has byid=TRUE
    result = suppressMessages(st_distance(geom1, geom2, which='Euclidean'))
    if(!byid){
        result = min(result)
    }else{
        result = t(result)
    }
    return(result)
}


#' limited replacement for rgeos::gIntersection using sf functionality
#'
#' @export
gIntersection<-function(spgeom1, spgeom2, byid=FALSE, drop_lower_td=FALSE){
    # Note drop_lower_td is not used here. 
    # From the documentation it looks like st_intersection always has the
    # equivalent of drop_lower_td=TRUE.

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    newgeom = suppressMessages(st_intersection(geom1, geom2))
    if(nrow(newgeom) == 0){
        return(NULL)
    }else{
        newgeom = st_geometry(newgeom)
    }

    if(!byid){
        newgeom = suppressMessages(st_union(newgeom))
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
    
}

#' limited replacement for rgeos::gIntersects using sf functionality
#'
#' @export
gIntersects<-function(spgeom1, spgeom2 = NULL, byid = FALSE, returnDense=TRUE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    result = suppressMessages(st_intersects(geom1, geom2, sparse=(!returnDense)))
    result = t(result)

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gArea using sf functionality
#'
#' @export
gArea<-function(spgeom, byid=FALSE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    # Rgeos did not correct for the projection. Can force that in st_area
    # by throwing away the CRS
    st_crs(geom) = NA

    area = suppressMessages(st_area(geom))

    if(!byid) area = sum(area)

    return(area)
}

#' limited replacement for rgeos::gContains using sf functionality
#'
#' @export
gContains<-function(spgeom1, spgeom2, byid=FALSE, prepared=TRUE, returnDense=TRUE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    result = suppressMessages(st_contains(geom1, geom2, prepared=prepared, sparse=(!returnDense)))
    result = t(result)

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gCovers using sf functionality
#'
#' @export
gCovers<-function(spgeom1, spgeom2 = NULL, byid=FALSE, returnDense=TRUE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    result = suppressMessages(st_covers(geom1, geom2, sparse=(!returnDense)))
    result = t(result)

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gUnaryUnion using sf functionality
#'
#' @export
gUnaryUnion<-function(spgeom, id=NULL){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    # for each value of id, c(st_union(geom[k1, ]), st_union(geom[k2, ]), ...)   
    if(!is.null(id)){
        tmp = suppressMessages(aggregate(geom, by=list(id), FUN=head, dissolve=TRUE))
        tmp = suppressMessages(st_geometry(tmp)) # For consistency with gUnaryUnion
    }else{
        tmp = suppressMessages(st_union(geom))
    }

    # Convert geometry from sf to sp
    outgeom = as(tmp, 'Spatial')
    return(outgeom)
}
