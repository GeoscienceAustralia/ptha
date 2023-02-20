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
#' Like rgeos, this treats all datasets as Cartesian. It can give small
#' differences in the results of buffered polygons compared to rgeos::gBuffer,
#' due to the internals of "sf::st_buffer".  While unimportant in most
#' circumstances, it can change results, e.g.  if it leads to a point near the
#' boundary being included/excluded. Suggest to use sf::st_buffer for non-legacy applications.
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
    # Ensure that we return geometries only
    newgeom = st_geometry(newgeom)

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    # The following is for consistency with rgeos
    if( (!byid) & (length(outgeom) == 1) ) row.names(outgeom) = 'buffer'

    return(outgeom)
}


#' limited replacement for rgeos::gCentroid using sf functionality
#'
#' Like rgeos, this treats all datasets as Cartesian. Use sf::st_centroid
#' for non-legacy applications.
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

    newgeom = st_geometry(newgeom)

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}

#' limited replacement for rgeos::gDistance using sf functionality
#'
#' Like rgeos, this treats all datasets as Cartesian. Use sf::st_distance
#' for non-legacy applications.
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
#' Beware this can lead to different ordering of the geometries as compared to
#' rgeos::gIntersection. Suggest to use sf::st_intersection for non-legacy applications.
#'
#' @export
gIntersection<-function(spgeom1, spgeom2, byid=FALSE, drop_lower_td=FALSE){

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

    geom1$TEMPLOCALVARTOREMOVE_id1 = row.names(spgeom1)
    geom2$TEMPLOCALVARTOREMOVE_id2 = row.names(spgeom2)

    newgeom = suppressMessages(st_intersection(geom1, geom2))
    if(nrow(newgeom) == 0){
        return(NULL)
    }
    # Fix rownames based on code in rgeos:::RGEOSBinTopoFunc
    newgeom_rownames = paste(newgeom$TEMPLOCALVARTOREMOVE_id1, newgeom$TEMPLOCALVARTOREMOVE_id2)
    newgeom = st_geometry(newgeom)
    
    if(drop_lower_td){
        # if drop_lower_td == TRUE, objects will be dropped from
        #  output GEOMETRYCOLLECTION objects to simplify output if their
        #  topological dinension is less than the minimum topological
        #  dinension of the input objects.
        st_dimension_flag = sapply(newgeom, function(x) st_dimension(x)) 
        k = which(st_dimension_flag == max(st_dimension_flag))
        newgeom = newgeom[k]
        newgeom_rownames = newgeom_rownames[k]
    }

    if(!byid){
        newgeom = suppressMessages(st_union(newgeom))
        # Fix rownames based on code in rgeos:::RGEOSBinTopoFunc
        newgeom_rownames = "1"
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    if(length(outgeom) == length(newgeom_rownames)){
        row.names(outgeom) = newgeom_rownames
    }

    return(outgeom)
    
}

#' limited replacement for rgeos::gIntersects using sf functionality
#'
#' Suggest to use sf::st_intersects for non-legacy applications.
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
    # Here rgeos would add row.names, like row.names(result) = row.names(spgeom2)
    # For now I don't think we need this and would need to check corner cases,
    # so not implementing for now.

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gArea using sf functionality
#'
#' Like rgeos this treats all geometries as Cartesian. Suggest to use
#' sf::st_area for non-legacy applications.
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
#' Like rgeos this assumes all coordinates are Cartesian. Suggest to use
#' sf::st_contains for non-legacy applications.
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
    # Here rgeos would add row.names, like row.names(result) = row.names(spgeom2)
    # For now I don't think we need this and would need to check corner cases,
    # so not implementing for now.

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gCovers using sf functionality
#'
#' Like rgeos this assumes all geometries are Cartesian. Suggest to use
#' sf::st_covers for non-legacy applications.
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
    # Here rgeos would add row.names( .. ) that I haven't fixed here.

    if(!byid){
        result = any(apply(result, 1, any)) 
    }

    return(result)
}

#' limited replacement for rgeos::gUnaryUnion using sf functionality
#'
#' Suggest to use sf::aggregate(..., dissolve=TRUE) and/or sf::st_union for non-legacy applications.
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
    }else{
        tmp = suppressMessages(st_union(geom))
    }
    tmp = suppressMessages(st_geometry(tmp)) # For consistency with gUnaryUnion

    # Convert geometry from sf to sp
    outgeom = as(tmp, 'Spatial')

    # Here rgeos would add row.names( .. ) that I haven't fixed here.

    return(outgeom)
}


#' limited replacement for rgeos::gUnion using sf functionality
#'
#' Suggest to use sf::st_union for non-legacy applications, perhaps
#' with prior calls to sf::aggregate(..., dissolve=TRUE) in the case
#' where byid=FALSE.
#'
#' @export
gUnion<-function(spgeom1, spgeom2, byid = FALSE){
    #, id = NULL,
    #  drop_lower_td = FALSE, 
    #  unaryUnion_if_byid_false = TRUE, 
    #  checkValidity = NULL

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = suppressMessages(sf_use_s2())
    suppressMessages(sf_use_s2(FALSE))
    on.exit(suppressMessages(sf_use_s2(using_s2)))

    if(!byid){
        # By default unaryUnion_if_byid_false = TRUE
        spgeom1 = gUnaryUnion(spgeom1)
        if(!is.null(spgeom2)){
            spgeom2 = gUnaryUnion(spgeom2)
        }
    }

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    newgeom = suppressMessages(st_union(geom1, geom2, by_feature=byid) )
    newgeom = st_geometry(newgeom)

    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}
