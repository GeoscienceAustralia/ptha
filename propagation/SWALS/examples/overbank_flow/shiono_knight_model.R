#' Analytical solution for steady-flow in the cross-channel direction
#'
#' Solve for the velocity distribution in a trapezoidal channel with floodplains,
#' using the Shiono/Knight (1991) model,
#' ignoring the slope-related drag coefficient and secondary flow:
#'     ghS0 - f/8 U^2 + d/dy{ lambda h^2 (f/8)^0.5 * 0.5 * d/dy(U^2)  } = 0
#' where 
#'     f = n^2 8 g d^(-1/3)
#' and n is Manning's n.
#' The channel is symmetric about the mid-channel (where y=0)
#'
shiono_knight_solution=function(
    width, # Domain width
    fraction_floodplain, # Fraction of the width that is floodplain
    fraction_banks, # Fraction of the width that is sloping banks
    fraction_flatbed, # Fraction of the width that is flatbed
    floodplain_elev, # Elevation of the floodplain
    flatbed_elev, # Elevation of the flat beded channel (deepest region)
    free_surface, # Elevation of the free surface
    S0,  # Downstream slope of floodplain
    lambda, # Diffusion coefficient
    manning_n, # Manning coefficient
    g, # Gravity
    dx # Numerical dx spacing
    ){

    stopifnot(isTRUE(all.equal(fraction_floodplain + fraction_banks + fraction_flatbed, 1.0)))

    stopifnot(free_surface > floodplain_elev)
    stopifnot(floodplain_elev > flatbed_elev)


    ## Y coordinate (cross-channel) -- 0 in the middle of the channel
    y = seq(-width/2, width/2, len=round((width/dx) + 1))

    ## Channel geometry
    bank_gradient = (floodplain_elev - flatbed_elev)/
        (width*(1-fraction_floodplain)/2 - width*(fraction_flatbed/2))

    elev = pmax(flatbed_elev, pmin(
                   floodplain_elev, 
                   flatbed_elev + bank_gradient * (abs(y) - width*fraction_flatbed/2)))

    depth = pmax(free_surface - elev, 0)
    f = manning_n^2 * 8 * g * depth^(-1/3)

    RHS = g * depth * S0
    DIAG = -f/8

    N = length(depth)
    depth_plus = c( 0.5*(depth[1:(N-1)] + depth[2:N]), 0) # Implicit boundary condition, depth --> 0
    depth_minus = c(0, 0.5*(depth[1:(N-1)] + depth[2:N])) # Implicit boundary condition, depth --> 0
    f_plus  = manning_n^2 * 8 * g * depth_plus^(-1/3)
    f_plus[!is.finite(f_plus)] = 0
    f_minus = manning_n^2 * 8 * g * depth_minus^(-1/3)
    f_minus[!is.finite(f_minus)] = 0
    diffuse_plus  = lambda * depth_plus**2  * sqrt(f_plus/8)  * 0.5 * 1/dx^2
    diffuse_minus = lambda * depth_minus**2 * sqrt(f_minus/8) * 0.5 * 1/dx^2

    # Use a symmetric banded matrix
    library(Matrix)
    M = bandSparse(length(depth), length(depth), # dimensions
                   (0):1, # band, diagonal is number 0. As the matrix is symmetric, we ignore the lower triangle
                   list(#diffuse_minus,
                        DIAG - diffuse_minus - diffuse_plus,
                        diffuse_plus),
                   symmetric=TRUE
                   )
    Usquared = solve(M, RHS)
    
    # Get the velocity
    U = sqrt(abs(Usquared))

    return(list(y=y, U=U, depth=depth, 
                # Discharge by trapezoidal integration
                discharge=sum(U[2:(N-1)]*depth[2:(N-1)]*dx) + sum(0.5*U[c(1,N)]*depth[c(1,N)]*dx)))
}


## INPUT PARAMETERS MATCHING THE SWALS MODEL RUN

# Geometry parameters
width = 1000
fraction_floodplain = 0.4
fraction_banks = 0.4
fraction_flatbed = 0.2
floodplain_elev = 0.0
flatbed_elev = -5.0

free_surface = floodplain_elev + 0.5

# Flow parameters
S0 = 0.001
lambda = 10
manning_n = 0.02
g = 9.8

## END INPUT PARAMETERS MATCHING THE SWALS MODEL RUN

# Compute solution without diffusion -- this will have discharge that matches
# the SWALS model discharge
solution_coarse_nolambda = shiono_knight_solution(
    width, # Domain width
    fraction_floodplain, # Fraction of the width that is floodplain
    fraction_banks, # Fraction of the width that is sloping banks
    fraction_flatbed, # Fraction of the width that is flatbed
    floodplain_elev, # Elevation of the floodplain
    flatbed_elev, # Elevation of the flat beded channel (deepest region)
    free_surface, # Elevation of the free surface
    S0,  # Downstream slope of floodplain
    0, # Diffusion coefficient
    manning_n, # Manning coefficient
    g, # Gravity
    20 # Numerical dx spacing
    )

# Get the solution with a specified discharge (instead of stage), using root-finding.
get_shiono_knight_solution_matching_discharge=function(
    width, # Domain width
    fraction_floodplain, # Fraction of the width that is floodplain
    fraction_banks, # Fraction of the width that is sloping banks
    fraction_flatbed, # Fraction of the width that is flatbed
    floodplain_elev, # Elevation of the floodplain
    flatbed_elev, # Elevation of the flat beded channel (deepest region)
    free_surface_lower_upper, # Vector with two free surface values, assume the true value is within this.
    S0,  # Downstream slope of floodplain
    lambda, # Diffusion coefficient
    manning_n, # Manning coefficient
    g, # Gravity
    dx, # Numerical dx spacing
    target_discharge # The discharge we need to match
    ){

    # Define a local solution that gives the discharge
    shiono_knight_discharge = function(free_surface, return_solution=FALSE){
        
        solution = shiono_knight_solution(
            width, # Domain width
            fraction_floodplain, # Fraction of the width that is floodplain
            fraction_banks, # Fraction of the width that is sloping banks
            fraction_flatbed, # Fraction of the width that is flatbed
            floodplain_elev, # Elevation of the floodplain
            flatbed_elev, # Elevation of the flat beded channel (deepest region)
            free_surface, # Elevation of the free surface
            S0,  # Downstream slope of floodplain
            lambda, # Diffusion coefficient
            manning_n, # Manning coefficient
            g, # Gravity
            dx # Numerical dx spacing
            )

        if(return_solution){
            return(solution)
        }else{
            return(solution$discharge - target_discharge)
        }
    }

    target_stage = uniroot(shiono_knight_discharge, interval = free_surface_lower_upper, tol=1.0e-08)$root

    solution_with_desired_discharge = shiono_knight_discharge(target_stage, return_solution=TRUE)

    return(solution_with_desired_discharge)
}

# Get the coarse-grid solution that matches the discharge of the solution without diffusion.
solution_coarse_matching_discharge = get_shiono_knight_solution_matching_discharge(
    width, # Domain width
    fraction_floodplain, # Fraction of the width that is floodplain
    fraction_banks, # Fraction of the width that is sloping banks
    fraction_flatbed, # Fraction of the width that is flatbed
    floodplain_elev, # Elevation of the floodplain
    flatbed_elev, # Elevation of the flat beded channel (deepest region)
    free_surface + c(-0.4,0.4), # Range of free surface where we search for the solution
    S0,  # Downstream slope of floodplain
    lambda, # Diffusion coefficient
    manning_n, # Manning coefficient
    g, # Gravity
    20, # Numerical dx spacing
    solution_coarse_nolambda$discharge
    )

## Get a numerical solution matching the SWALS model run
#solution_coarse = shiono_knight_solution(
#    width, # Domain width
#    fraction_floodplain, # Fraction of the width that is floodplain
#    fraction_banks, # Fraction of the width that is sloping banks
#    fraction_flatbed, # Fraction of the width that is flatbed
#    floodplain_elev, # Elevation of the floodplain
#    flatbed_elev, # Elevation of the flat beded channel (deepest region)
#    free_surface, # Elevation of the free surface
#    S0,  # Downstream slope of floodplain
#    lambda, # Diffusion coefficient
#    manning_n, # Manning coefficient
#    g, # Gravity
#    20 # Numerical dx spacing
#    )

## As above with a finer grid
#solution_fine = shiono_knight_solution(
#    width, # Domain width
#    fraction_floodplain, # Fraction of the width that is floodplain
#    fraction_banks, # Fraction of the width that is sloping banks
#    fraction_flatbed, # Fraction of the width that is flatbed
#    floodplain_elev, # Elevation of the floodplain
#    flatbed_elev, # Elevation of the flat beded channel (deepest region)
#    free_surface, # Elevation of the free surface
#    S0,  # Downstream slope of floodplain
#    lambda, # Diffusion coefficient
#    manning_n, # Manning coefficient
#    g, # Gravity
#    2.5 # Numerical dx spacing
#    )
## Check convergence -- yes it's good
# plot(solution_coarse$y, solution_coarse$U, t='o')
# points(solution_fine$y, solution_fine$U, col='red', t='l')
# points(solution_coarse_nolambda$y, solution_coarse_nolambda$U, col='green', t='l')


