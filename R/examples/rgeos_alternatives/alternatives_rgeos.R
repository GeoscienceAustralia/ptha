gBuffer<-function(spgeom, width, quadsegs=5, byid=FALSE){
    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = sf_use_s2()
    sf_use_s2(FALSE)    
    on.exit(sf_use_s2(using_s2))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    # st_buffer seems to be always applied on a per-id basis
    newgeom = st_buffer(geom, dist=width, quadsegs=quadsegs)
    # If we didn't want byID results, then merge
    if(!byid){
        newgeom = st_union(newgeom) 
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}


gCentroid<-function(spgeom, byid=FALSE){
    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = sf_use_s2()
    sf_use_s2(FALSE)    
    on.exit(sf_use_s2(using_s2))

    # Convert geometry from sp to sf
    geom = st_as_sf(spgeom)

    if(byid){
        newgeom = st_centroid(geom)
    }else{
       #https://gis.stackexchange.com/questions/451041/equivalent-of-gcentroidx-byid-false-in-sf-system 
        newgeom = st_union(geom)
        newgeom = st_centroid(newgeom)
    }

    # Convert geometry from sf to sp
    outgeom = as(newgeom, 'Spatial')
    return(outgeom)
}

gDistance<-function(spgeom1, spgeom2=NULL, byid=FALSE){

    # RGEOS assumed planar coordinates. Enforce that behaviour in sf
    using_s2 = sf_use_s2()
    sf_use_s2(FALSE)    
    on.exit(sf_use_s2(using_s2))

    # Convert geometry from sp to sf
    geom1 = st_as_sf(spgeom1)
    if(is.null(spgeom2)){
        geom2 = NULL
    }else{
        geom2 = st_as_sf(spgeom2)
    }

    # st_distance always has byid=TRUE
    result = drop_units(st_distance(geom1, geom2, which='Euclidean'))
    if(!byid){
        result = min(result)
    }else{
        result = t(result)
    }
    return(result)
}
