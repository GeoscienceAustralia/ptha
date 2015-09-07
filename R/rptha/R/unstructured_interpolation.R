######################################################################
#
# Here we have functions for nearest-neighbour and 'linear' (triangular)
# interpolation of unstructured point data
#
# Author: Gareth Davies, Geoscience Australia 2014
#
#
    

#' Function for nearest neighbour interpolation
#'
#' Basically a wrapper around the SearchTrees package
#'
#' @param xy = points locations associated with values to interpolation from (typically nx2 matrix)
#' @param vals = values at xy (can be a matix with 1 or more colums, and the same number of rows as xy)
#' @param newPts = points where we want interpolated values (typically mx2 matrix)
#' @return matrix/vector with as many columns as 'vals' and as many rows as 'newPts', containing the 'vals' interpolated to 'newPts'
#' @export
#' @import FNN
#' @examples
#'    # Make a single triangle in the plane z=x+y, and interpolate from it
#'    xy = matrix(c(0,0,0,1,1,1), ncol=2, byrow=TRUE)
#'    vals = c(0, 1, 2) # z=x+y
#'    newPts = matrix(c(0.1, 0.9, 0.9, 0.9), ncol=2, byrow=TRUE)
#'
#'    out = nearest_neighbour_interpolation(xy, vals, newPts)
#'    stopifnot(all.equal(out, c(1,2)))
nearest_neighbour_interpolation<-function(xy, vals, newPts){
    
    require(FNN)
    newInds = get.knnx(xy, newPts[,1:2, drop=FALSE], k=1)[[1]]

    if(is.null(dim(vals))){
        return(vals[newInds])
    }else{
        return(vals[newInds,])
    }
}

#' Function for various types of triangular ['linear'] interpolation of unstructured data. 
#' 
#' Delanay triangulation is supported, as well as a method based on linear interpolation 
#' from the 3 nearest neighbour of the interpolation point, with limiting. \cr
#' If useNearestNeighbour = FALSE then it provides a wrapper around the delanay triangulation used in the 'geometry' package.
#' Unfortunately the look-up can be slow with this method for large point clouds. \cr
#' If useNearestNeighbour=TRUE, we find the 3 nearest xy neighbours of each point to interpolate to, and
#' interpolate using the plane defined by those 3 neighbours. Limiting is used
#' to ensure the interpolated value does not exceed the range of the xy
#' neighbours. This method is fast since it relies only an a fast nearest neighbours implementation (via FNN)
#'
#' @param xy = point locations associated with the values to interpolate from  (typically nx2 matrix)
#' @param vals = values at xy (can be a matix with 1 or more colums and the same number of rows as xy)
#' @param newPts = points where we want interpolated values (typically mx2 matrix)
#' @param useNearestNeighbour = TRUE/FALSE (effect described above)
#' @return matrix/vector with as many columns as 'vals' and as many rows as 'newPts', containing the 'vals' interpolated to 'newPts'
#' @export
#' @import FNN 
#' @import geometry
#' @examples
#'    # Make a single triangle in the plane z=x+y, and interpolate from it
#'    xy = matrix(c(0, 0, 0, 1, 1, 1), ncol=2, byrow=TRUE)
#'    vals = c(0, 1, 2) # z=x+y
#'    newPts = matrix(c(0.5, 0.5, 0.3, 0.3), ncol=2, byrow=TRUE)
#'
#'    out = triangular_interpolation(xy, vals, newPts)
#'    stopifnot(all.equal(out, c(1.0,0.6)))
#'
#'    # Re-order triangle
#'    xy = xy[3:1,]
#'    vals = vals[3:1]
#'    out = triangular_interpolation(xy, vals, newPts)
#'    stopifnot(all.equal(out,c(1.0,0.6)))
#'
#'    #another one, with formula z=0.5*x+0.2*y+7
#'    xy = matrix(c(-1, -1, 1, -0.5, 0.5, 1), ncol=2,byrow=2)
#'    vals = 0.5*xy[,1]+0.2*xy[,2]+7
#'    newPts = matrix(c(0,0, 0.5, 0.3),ncol=2,byrow=TRUE)
#'    expectedVals = 0.5*newPts[,1]+0.2*newPts[,2]+7
#'    out = triangular_interpolation(xy,vals,newPts)
#'    stopifnot(all.equal(out,expectedVals))
#'
#'    # A point outside the triangle 
#'    newPts = matrix(c(-1,0, -1, 1), ncol=2, byrow=TRUE)
#'    out = triangular_interpolation(xy, vals, newPts, useNearestNeighbour=FALSE)
#'    stopifnot(all(is.na(out)))
#'    # Failure is expected here if using approximate triangulation based on nearest neighbour methods
#'
#'    # A single point
#'    newPts = matrix(c(0,0), ncol=2)
#'    out = triangular_interpolation(xy, vals, newPts)
#'    stopifnot(out == 7)
#'
#'    # Points on the triangle
#'    newPts = xy
#'    out = triangular_interpolation(xy, vals, newPts)
#'    stopifnot(all(out == vals))
#'
#'    # Point on an edge
#'    newPts = matrix(0.5*(xy[1,]+xy[2,]), ncol=2)
#'    out = triangular_interpolation(xy, vals, newPts)
#'    stopifnot(all(out == 0.5*(vals[1]+vals[2])))
triangular_interpolation<-function(xy, vals, newPts, useNearestNeighbour=TRUE){

    xy = as.matrix(xy,ncol=2)
    newPts = as.matrix(newPts)

    if(!is.null(dim(vals))){
        vals = as.matrix(vals)
    }
 
    if(dim(xy)[1]<3 | is.null(dim(newPts))){
        stop('Need at least 3 input xy points, and newPts should be a matrix')
    }

    if(is.null(dim(vals))){
        if(length(vals) != length(xy[,1])){
            stop('Length of xy[,1] and vals must be the same')
        }
    }else{
        if(length(vals[,1]) != length(xy[,1])){
            stop('Length of xy[,1] and vals[,1] must be the same')
        }
        
    }
 
    # Flag to say whether we use nearest-neighbours to define the triangulation
    #nnSort=(useNearestNeighbour & length(xy[,1])>3)
    nnSort = useNearestNeighbour

    if(!nnSort){
        # Use geometry package triangulation
        require(geometry)
        # This is slow for large problems [tsearch is presently slow], but
        # arguably the best approach
        triIndices = delaunayn(xy)
        # Use barycentric coordinates from tsearch
        triOn = tsearch(xy[,1], xy[,2], triIndices, newPts[,1], newPts[,2], bary=TRUE)
        if(is.null(dim(vals))){
            # Vals is a vector
            final = vals[triIndices[triOn$idx,1]]*triOn$p[,1]+
                vals[triIndices[triOn$idx,2]]*triOn$p[,2]+
                vals[triIndices[triOn$idx,3]]*triOn$p[,3]

        }else{
            # Vals is a matrix
            final = matrix(NA,ncol=ncol(vals), nrow=nrow(newPts))
            for(i in 1:ncol(final)){
                final[,i] = vals[triIndices[triOn$idx,1],i]*triOn$p[,1]+
                    vals[triIndices[triOn$idx,2],i]*triOn$p[,2]+
                   vals[triIndices[triOn$idx,3],i]*triOn$p[,3]
            }
        }

    }else{
        # Hack to try to speed-up tsearch on large problems.
        # Instead of using t-search, find the 3 nearest neighbours and make a triangle
        #triTree = createTree(xy)
        # Lookup the nearest index on the tree
        #lookupInds = knnLookup(triTree, newx=newPts[,1],newy=newPts[,2],k=3)
        require(FNN)
        lookupInds = get.knnx(xy, newPts[,1:2, drop=FALSE], k=3)[[1]]

        ## Interpolate. 
        ## Get vertices
        p1 = matrix(xy[lookupInds[,1],], ncol=2)
        p2 = matrix(xy[lookupInds[,2],],ncol=2)
        p3 = matrix(xy[lookupInds[,3],], ncol=2)
       
        ## Get triangle gradient  
        dx31 = p3[,1]-p1[,1]
        dy31 = p3[,2]-p1[,2]
        dx21 = p2[,1]-p1[,1]
        dy21 = p2[,2]-p1[,2]

        dxN = newPts[,1]-p1[,1]
        dyN = newPts[,2]-p1[,2]
        #
        ## Compute triangle area
        area = dx21*dy31-dx31*dy21 #dy21*dx31 - dx21*dy31

        # Gradient coefficients
        a = (dy31*dxN-dx31*dyN)/area 
        b = (-dy21*dxN+dx21*dyN)/area 

        # Treat cases with degenerate triangles -- use nearest-neighbour instead
        EPS = 1.0e-06
        a[abs(area)<EPS] = 0.
        b[abs(area)<EPS] = 0.


        if(is.null(dim(vals))){
            # Find max/min 'vals' on triangle
            valsMax = pmax(vals[lookupInds[,1]], vals[lookupInds[,2]], vals[lookupInds[,3]])
            valsmin = pmin(vals[lookupInds[,1]], vals[lookupInds[,2]], vals[lookupInds[,3]])
            
            dz31 = vals[lookupInds[,3]] - vals[lookupInds[,1]]
            dz21 = vals[lookupInds[,2]] - vals[lookupInds[,1]]

            final = vals[lookupInds[,1]] + a*dz21 + b*dz31

            # Limit
            M = (final>valsMax)
            m = (final<valsmin)
            limit = pmax(M,m)
            final = final*(1-limit) + vals[lookupInds[,1]]*limit

        }else{
            # Vals is higher dimensional

            # Find max/min 'vals' on triangle
            valsMax = pmax(vals[lookupInds[,1],], vals[lookupInds[,2],], vals[lookupInds[,3],])
            valsmin = pmin(vals[lookupInds[,1],], vals[lookupInds[,2],], vals[lookupInds[,3],])

            dz31 = vals[lookupInds[,3],] - vals[lookupInds[,1],]
            dz21 = vals[lookupInds[,2],] - vals[lookupInds[,1],]

            final = vals[lookupInds[,1],] + a*dz21 + b*dz31

            # Limit
            M = (final>valsMax)
            m = (final<valsmin)
            limit = pmax(M,m)
            # If outside min/max, use nearest neighbour only
            final = final*(1-limit) + vals[lookupInds[,1],]*limit

        }
    }
    return(final)
}


