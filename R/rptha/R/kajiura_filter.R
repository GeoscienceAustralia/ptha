# Code to apply a Kajiura filter to the earth surface deformation
#
# Gareth Davies, Geoscience Australia, 2013+


#' Compute Kajiura's filter function
#'
#' The function value is the sum over n in [0, Infty] of
#' [ 1/pi*sum( (-1)^n*(2*n+1)/( (2*n+1)^2 + r^2)^(3/2)) ].
#' For a more efficient and vectorized implementation, see kajiura_g_empirical
#' which approximates this function using splines and a log transform.
#'
#' @param r numeric. Value of r in the equation
#' @param nt Integer. Number of terms. We compare the sum of the first nt terms, and the
#' next nt terms, and if they differ by less than 'recursive_stop_factor' we
#' return the combined sum. Otherwise we double nt and try again.
#' @param verbose logical. Print when recurisng
#' @param recursive_stop_factor numeric. See information on 'nt' above
#' @return The value of the function at r
#' @export
#' 
kajiura_g<-function(r, nt=1000, verbose=FALSE, recursive_stop_factor=1.0e-05){

    # Compute kajiura's filter function
    if(length(r)>1) stop('r can only be length 1')
    # 

    # The filter has an infinite series representation
    n = seq(0,nt)
    g = 1/pi*sum( (-1)^n*(2*n+1)/( (2*n+1)^2 + r^2)^(3/2))

    n2 = seq(nt+1, 2*nt)
    g2 = 1/pi*sum( (-1)^n2*(2*n2+1)/( (2*n2+1)^2 + r^2)^(3/2))
        
    if(verbose) print(c(g2,g))

    if(abs(g2) > abs(g)*recursive_stop_factor){
        # Ensure that we have enough terms
        if(verbose) print('recursing')
        g = kajiura_g(r, nt=2*nt, verbose=verbose)
    }
    return(g+g2) 
}

#' Get a function to efficiently evaluate Kajiura's G function
#'
#' This function returns a function which can be used to efficiently
#' evaluate Kajiura's G function. It is computed by tabulating the function
#' at 'n' points between r = 0 and r=rMax, and then fitting a spline to the
#' log-transformed result. The latter spline is used to lookup log(kajiura_g), 
#' and the exponential transform is applied to correct the result. Our tests
#' suggest the approach is highly accurate, and in R it is much more efficient than
#' direct computation.
#'
#' @param rMax maximum value of r at which lookup is desired
#' @param n number of points used to empirically approximate the log-function
#' values.
#' @return Function to compute kajiura_G from an input vector of R values <= rMax.
#' @export
kajiura_g_empirical<-function(rMax=9, n=81){
    # Empirical approximation to kajiura_g
    #
    # Uses splines to log-transformed data
    #
    # It's good
    #
    kg = Vectorize(kajiura_g)
    x = seq(0,rMax,len=n)
    kg_x = kg(x)
    f_kx = splinefun(x, log(kg_x))
    kgE<-function(x) exp(f_kx(x))
    return(kgE)
}



#' Implementation of a Kajiura filter 
#'
#' Implement the filter similar to that of Glimsdal et al (2013, Eqns 10-11), based on Kajiura (1963). \cr
#' \deqn{ newDeformation(x,y) = depth(x,y)^{-2} \int \int \big{[}oldDeformation(x',y')G( \sqrt{(x'-x)^2+(y'-y)^2} / depth(x, y) )\big{]} dx' dy'}
#' This is a 2D generalisation of the cosh filter, justified for a 'temporally short'
#' earthquake with ocean governed by linear wave equations in constant depth
#' water. The definition of G is explained further in the help for \code{kajiura_g}\cr
#' Essentially: \cr
#' xyDef[,3] <-- convolution of the "old deformation" and "a smoothing kernel which integrates to 1"\cr
#' Numerically, we compute the new deformation as: \cr
#' \deqn{ numerator_{nm} = sum_{i} sum_{j} oldDeformation_{ij} * G( \sqrt{ (x_{nm} - x_{ij})^2 + (y_{nm} - y_{ij})^2} / depth_{nm} ) }
#' \deqn{ denominator_{nm} = sum_{i} sum_{j} G( \sqrt{ (x_{nm} - x_{ij})^2 + (y_{nm} - y_{ij})^2} / depth_{nm} ) }
#' \deqn{ newDeformation_{nm} = numerator_{nm} / denominator_{nm} }
#' where nm and ij denote pixel coordinates on a regular (cartesian) grid. So at any point
#' nm, the new deformation is a weighted average of nearby values, with weights
#' coming from Kajiuras function.\cr
#' This is slightly different to the approach used in Glimsdal et al (2013), but is
#' equivalent for constant depth. With non-constant depth the underlying theory is
#' not exactly valid, but should provide a reasonable approximation for slowly varying depths.\cr
#' We attempt to reduce edge effects by linearly weighting original and filtered values at edges,
#' since we cannot efficiently deal with edge effects in a better way. Therefore
#' it is best to have unimportant features around the edge of the input points. \cr
#' We allow xyDef to be unstructured, and start by gridding the results
#' on a grid with spacing approximately grid_dx,grid_dy.
#' The grid spacing is not exactly grid_dx,grid_dy, because it is forced to 
#' exactly divide the maximum input depth. \cr
#' For deformations with discontinuities, there can be artefacts due to
#' regridding, and it may be numerically beneficial to rotate the x,y input
#' coordinates so that the discontinuity is aligned with one of the coordinate
#' axes. (This is done in make_tsunami_unit_source).
#' 
#' @param xyDef  3 column matrix with x,y, deformation. x and y can be unstructured.
#' The x,y coordinate system must be cartesian and in m
#' @param depth vector with the depth at each x,y point in xyDef (in m)
#' @param grid_dx Numeric (m). See grid_dy
#' @param grid_dy Numeric (m). To apply the filter, we regrid xyDef on a grid with
#' point spacing grid_dx, grid_dy, then smooth, then transform back from the grid
#' to our original xy points. This should ideally be a fraction of the depth in
#' areas of significant deformation.
#' @param edge_buffer_value Numeric. Outside the domain edges we assume this is
#' the value of xyDef[,3], when the filter is applied. Without further correction, this
#' would make edges tend towards edge_buffer_value
#' @param edge_effect_correction_scale To reduce edge effects, after
#' filtering, we compute the distance of every point to the edge of the domain
#' 'd', and then return the solution \cr
#' WT = max( 1-d/(reference_depth*edge_effect_correction_scale), 0)**0.5 \cr
#' OUTPUT = WT*(ORIGINAL OUTPUT) + (1-WT)*FILTERED OUTPUT \cr
#' @param interpolator 'linear' or 'nearest'. Linear is better, but may be slow for
#' large point clouds. Don't use nearest unless grid_x, grid_y are 'small enough'
#' @param interpolator_categories Function of xy (matrix of coordinates)
#' returning a value for each point, with distinct values corresponding to
#' distinct 'groups' of points. Interpolation will be performed separately on
#' each group. It is passed as \code{category_function} to \code{interpolation_discontinuous}
#' @param kajiuraGmax When empirically approximating kajiuraG, we fit it from
#' x=[0, kajiuraGmax]. Values above this are evaluated to zero. Be cautious about changing this.
#' @param volume_change_error_threshold If the difference in the positive or
#' negative or total volume before and after filtering, relative to the original
#' 'absolute' volume, is more than this, then throw an error. This might not indicate a mistake,
#' but it does indicate a large change in the deformation, which is worth investigating.
#' Possible causes include strong depth variation, or regridding with coarse grid_dx, grid_dy.
#' @param volume_change_MSL The mean-sea-level used to define positive/negative deformation
#' for the volume check mentioned above.
#' @param verbose Print lots of information about the fit
#' @param fortran_inner_loop Use compiled fortran for the inner loop instead of R
#' @return replacement version of xyDef, with smoothing applied to xyDef[,3]
#' @export
kajiura_filter<-function(xyDef, 
                         depth,
                         grid_dx=max(depth)/2,
                         grid_dy=max(depth)/2, 
                         edge_buffer_value=0, 
                         edge_effect_correction_scale=1.5, 
                         kajiuraGmax=9,
                         interpolator='linear',
                         interpolator_categories = function(xy){xy[,1]*0},
                         volume_change_error_threshold=0.02,
                         volume_change_MSL=0.0,
                         verbose=FALSE,
                         fortran_inner_loop=TRUE){


    reference_depth = max(depth)

    # Get fast approximation to G function
    kgE = kajiura_g_empirical(rMax=kajiuraGmax)

    # dx/dy for gridded data + filter. Ensure reference_depth is an integer
    # multiple of dx/dy
    dx = reference_depth/ceiling(reference_depth/grid_dx)
    dy = reference_depth/ceiling(reference_depth/grid_dy)

    m0 = min(xyDef[,1])
    m1 = max(xyDef[,1])
    newX = seq(m0,m1,len=round((m1-m0)/dx)+1)
    m0 = min(xyDef[,2])
    m1 = max(xyDef[,2])
    newY = seq(m0,m1,len=round((m1-m0)/dy)+1)
    lny = length(newY)
    lnx = length(newX)

    # Compute nearest-neighbour interpolation
    newPts = as.matrix(expand.grid(newX,newY))

    if(verbose) print('Unstructured interpolation number 1...')
    if(interpolator=='nearest'){
        # Nearest neighbour interpolation
        #interp1 = nearest_neighbour_interpolation(xyDef[,1:2], 
        #    cbind(xyDef[,3], depth), newPts)
        interp1 = interpolation_discontinuous(xyDef[,1:2],
            cbind(xyDef[,3], depth), newPts, 
            category_function=interpolator_categories,
            interpolator=nearest_neighbour_interpolation)

        newVals = matrix(interp1[,1], ncol=lnx,byrow=T)
        newDepth = matrix(interp1[,2], ncol=lnx,byrow=T)

        rm(interp1); gc()

    }else if(interpolator=='linear'){
        # Cheap triangulation/3-point-knn type interpolation
        #interp1 = triangular_interpolation(xyDef[,1:2], cbind(xyDef[,3], depth),
        #    newPts)
        interp1 = interpolation_discontinuous(xyDef[,1:2], 
            cbind(xyDef[,3], depth), newPts, 
            category_function=interpolator_categories,
            interpolator=triangular_interpolation, 
            useNearestNeighbour=TRUE)

        ## Remove NA values if they exist (e.g. from points outside of a
        ## delaunay triangulation, if that type of interpolation is applied).
        ## (no longer needed)
        # interp1[is.na(interp1[,1]), 1] = edge_buffer_value
        # interp1[is.na(interp1[,2]), 2] = 1.e-12 # Set depth to a small number

        newVals = matrix(interp1[,1], ncol=lnx, byrow=T)
        newDepth = matrix(interp1[,2], ncol=lnx, byrow=T)
        
        rm(interp1); gc()

    }else{
        stop('interpolator not recognized')
    }
   
    # Compute kajiura_g filter function, on a matrix varying from +- 5 reference depths
    # Glimsdal et al (2013) highlight that this only needs to be 5 reference depths long
    filter_refDepth_range = 5
    fR = filter_refDepth_range*reference_depth
    filterXs = seq(-fR,fR, by=dx)
    filterYs = seq(-fR,fR, by=dy)
    lfx = length(filterXs)
    lfy = length(filterYs)

    if((length(filterXs)%%2 != 1)|
       (length(filterYs)%%2 != 1)){
        stop('ERROR: Filter size is wrong')
    }

    # Compute 'radius' term on filter
    filterXY = expand.grid(filterXs,filterYs)
    filterXYr = matrix( (filterXY[,1]^2+filterXY[,2]^2)**0.5,
        ncol=lfx, byrow=TRUE)
    rm(filterXY); gc()

    if(verbose) print('Applying filter ...') 

    # To apply the filter, we 'pad' the data with zeros (which means the edges
    # will taper down). 
    # The matrix has length(filterYs) rows of zeros at the top an bottom, and 
    #    length(filterXs) columns of zeros at the left and right
    filVals = matrix(edge_buffer_value, ncol=lnx+2*lfx, nrow=lny+2*lfy)
    filVals[lfy + 1:lny, lfx + 1:lnx] = newVals

    old_newVals = newVals
    # Now set newVals to zero -- it will hold the filtered results
    newVals = 0.*newVals
    GtermsSum = 0.*newVals

    # The filter is of length lfx,lfy
    # Generally lfx<<lnx, lfy<<lny
    # So for efficiency, here we loop over every element of the filter when
    # computing the weighted average
    ## THIS LOOP IS THE MOST COMPUTATIONALLY INTENSIVE PART OF THIS ROUTINE 
    depth_inv = 1.0/pmax(newDepth, 1.0e-20)
    if(!fortran_inner_loop){

        for(i in 1:lfx){
            if(verbose) print(paste0(i, ' of ', lfx ))
            for(j in 1:lfy){
                # Compute r/depth for the j,i cell of the filter,  avoid division by zero
                # Store it as G_j_i to reduce memory usage
                G_j_i = filterXYr[j,i] * depth_inv
                G_j_i = G_j_i*(G_j_i < kajiuraGmax) + kajiuraGmax*(G_j_i >= kajiuraGmax)

                # Put into a matrix which aligns with newVals
                G_j_i = kgE(G_j_i) * (G_j_i < kajiuraGmax) 
                dim(G_j_i) = c(lny, lnx)

                # Numerator of the weighted average
                newVals = newVals + filVals[(j+(lfy-1)/2) + 1:lny,(i+(lfx-1)/2) + 1:lnx] * G_j_i

                # Denominator of the weighted average. 
                GtermsSum = GtermsSum + G_j_i
            }
        }

        # Compute final weighted average 
        newVals = newVals/GtermsSum

        rm(G_j_i, GtermsSum, depth_inv); gc()
    }else{
        # Does the same as above, in Fortran! This will modify 'newVals'
        kajiura_convolution(depth_inv, filVals, filterXYr, 
            kajiuraGmax, lfx, lfy, lnx, lny, newVals)

        rm(depth_inv); gc()
    }

    if(edge_effect_correction_scale > 0.){
        if(verbose) print('Reducing edge effects ...')
        # Reduce edge-effects with a weighted average of the old values there
        # Compute the distance from an edge in the x/y directions
        xEdge = matrix( pmin((1:lnx)-0.5, (lnx:1)-0.5)*dx, byrow=TRUE, ncol=lnx, nrow=lny)
        yEdge = matrix( pmin((1:lny)-0.5, (lny:1)-0.5)*dy, byrow=FALSE, ncol=lnx, nrow=lny)
        # Convert to a weight
        xEdge = pmax(1-xEdge/(fR*edge_effect_correction_scale), 0.)**0.5
        yEdge = pmax(1-yEdge/(fR*edge_effect_correction_scale), 0.)**0.5
        edgeF = pmax(xEdge,yEdge)
        # Take weighted average of the original values and the new ones
        newVals = edgeF*old_newVals + (1-edgeF)*newVals

        rm(xEdge, yEdge, edgeF); gc()
    }
   

    # Back-compute x,y,Def with nearest neighbour lookup
    if(verbose) print('Unstructured interpolation number 2...')

    if(interpolator == 'nearest'){
        # At this stage, smoothing should have made the deformation continuous
        interp2 = nearest_neighbour_interpolation(newPts, c(t(newVals)), 
            xyDef[,1:2])

        #interp2 = interpolation_discontinuous(newPts, c(t(newVals)), 
        #    xyDef[,1:2], category_function=interpolator_categories,
        #    interpolator=nearest_neighbour_interpolation)
        new_xyDef = xyDef
        new_xyDef[,3] = interp2 

        rm(interp2); gc()

    }else if(interpolator == 'linear'){
        # At this stage, smoothing should have made the deformation continuous
        interp2 = triangular_interpolation(newPts, c(t(newVals)), xyDef[,1:2])

        #interp2 = interpolation_discontinuous(newPts, c(t(newVals)), 
        #    xyDef[,1:2], category_function=interpolator_categories,
        #    interpolator=triangular_interpolation, useNearestNeighbour=TRUE)
        #interp2[is.na(interp2)] = edge_buffer_value
        new_xyDef = xyDef
        new_xyDef[,3] = interp2 

        rm(interp2); gc()
    }

    
    newValsPosSum = sum( (newVals - volume_change_MSL)*(newVals>volume_change_MSL))
    newValsNegSum = sum( (newVals - volume_change_MSL)*(newVals<volume_change_MSL))
    old_newValsPosSum = sum((old_newVals - volume_change_MSL)*(old_newVals>volume_change_MSL))
    old_newValsNegSum = sum((old_newVals - volume_change_MSL)*(old_newVals<volume_change_MSL))

    oldNewvalsSum = sum(old_newVals)
    newvalsSum = sum(newVals)

    # Total volume of displacement (both + and - areas are counted as + in this
    # volume -- because otherwise, a deformation with equal positive and negative
    # parts would have a zero total volume)
    vol_reference = old_newValsPosSum - old_newValsNegSum

    r1 = (newValsPosSum-old_newValsPosSum)/vol_reference 
    r2 = (newValsNegSum-old_newValsNegSum)/vol_reference 
    r3 = (sum(newVals)-sum(old_newVals))/vol_reference 

    if(verbose){
        print(paste('Original re-gridded volume: ', oldNewvalsSum))
        print(paste('New gridded volume: ', newvalsSum))
        print(paste('Total volume relative change: ', r3))
        print(' ')
        print(paste('Original positive re-gridded volume: ', old_newValsPosSum))
        print(paste('New positive gridded volume: ', newValsPosSum))
        print(paste('Positive volume relative change: ', r1))
        print(' ')
        print(paste('Original negative re-gridded volume: ', old_newValsNegSum))
        print(paste('New negative gridded volume: ', newValsNegSum))
        print(paste('Negative volume relative change: ', r2))
        print(' ')
    }

    # Check for gross volume conservation errors
    if(is.finite(r3)){
        if(abs(r3) > volume_change_error_threshold){
            print(paste('r3 = ', r3))
            stop('Volume change threshold exceeded')
        }
    }else{
        # Look for other measures of problems
        # Note though that smoothing can decrease both the positive
        # and negative volume
        if(is.finite(r1)){
            if(abs(r1) > volume_change_error_threshold){
                print(paste('r1 = ', r1))
                stop('(Positive vol) Volume change threshold exceeded')

            }
        }
        if(is.finite(r2)){
            if(abs(r2) > volume_change_error_threshold){
                print(paste('r2 = ', r1))
                stop('(Negative vol) Volume change threshold exceeded')

            }
        }

    }

    return(new_xyDef)
     
}
