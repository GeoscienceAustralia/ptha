
#' Get the default sffm model parameters
#'
#' Convenience function to make the sffm_default_model_parameters list. By default
#' these parameters are used in all sffm_ functions, unless an alternative list
#' of default parameters is provided. See the help for \code{sffm_simulate} for
#' an example of using non-default parameters. Default parameters correspond to the
#' S_{NCF} model of Davies et al., 2015
#'
#' @return a list containing default parameters affecting the sffm generation \cr
#' random_phase_generating_field --  random function f(x) which returns
#' random numbers of length(x) for the phase generating field \cr
#' negative_slip_removal_function -- a function f(x) which returns x with negative
#' values made zero or positive \cr
#' spatial_slip_decay -- Type of spatial slip decay to apply. Options
#' are 'gaussian' or 'exponential' or 'none' \cr
#' recentre_slip -- Logical. Should the location of maximum slip be
#' adjusted to correspond to the maximum slip location in the input raster? \cr
#' spectral_amplitude_fun -- Function f(kx, ky, reg_par) to return the spectral amplitude given
#' kx, ky and a vector of model parameters reg_par. The first two values of
#' reg_par should  be the corner wavenumbers kcx, kcy IN NUMERICAL SPACE (i.e. the
#' physical wavenumbers multiplied by dx and dy respectively) \cr
#' 
#' @references
#' Davies et al. (2015), 
#' Tsunami inundation from heterogeneous earthquake slip distributions:
#' Evaluation of synthetic source models, J. Geophys. Res. Solid Earth, 120,
#' 6431-6451, doi:10.1002/2015JB012272. \cr
#'
#'
#' @export
#' @examples
#' sffm_default_par = sffm_get_default_model_parameters()
sffm_get_default_model_parameters<-function(){

    # random function f(x) which returns random numbers of length(x) for the
    # phase generating field
    random_phase_generating_field = rnorm

    # function f(x) which returns x with negative values removed in some way
    negative_slip_removal_function = function(x) pmax(x, 0)

    # Type of spatial slip decay to apply. Options are 'gaussian' or
    # 'exponential' or 'none'
    spatial_slip_decay = 'gaussian'

    # Should the location of maximum slip be adjusted to correspond to 
    # the maximum slip location in the input raster
    recentre_slip = TRUE

    # Function to return the spectral amplitude given kx, ky and a vector
    # of model parameters reg_par. The first two values of reg_par should 
    # be the corner wavenumbers kcx, kcy IN NUMERICAL SPACE (i.e. the physical
    # wavenumbers multiplied by dx and dy respectively)
    spectral_amplitude_fun<-function(kx, ky, reg_par){
       return(
           (1.0 + ((kx/reg_par[1])**2 + (ky/reg_par[2])**2)**2)**(-0.5)
       )
    }
    return(as.list(environment()))
}


# Hold configuration variables in a list
.sffm_default_model_parameters = sffm_get_default_model_parameters()


###############################################################################
#'
#' Return the current value of .Random.seed 
#'
#' Utility function to get the random seed -- and make it first if required
#'
#' @export
#' @examples
#'  x = get_random_seed()
#'
get_random_seed<-function(){

    if(!exists('.Random.seed', where=.GlobalEnv)){
        # Create the random seed by calling a random number generator
        x = rnorm(1)
    }

    return(.Random.seed)
}


###############################################################################
#'
#' Compute wavenumbers in numerical space for stochastic slip modelling 
#'
#' If tg_mat has dimensions N x M, then the numerical wavenumbers are: \cr
#' kx' = pmin(0:(M-1), N - (0:(M-1)))/M \cr
#' ky' = pmin(0:(N-1), N - (0:(N-1)))/N \cr
#' The function returns kx, ky as matrices with the same dimension as
#' tg_mat. kx has rows all equal kx', and ky has columns all equal ky'. \cr
#' These numerical wavenumbers can be adjusted to physical wavenumbers by
#' division by cellsize dx for kx, or dy for ky.
#'
#' @param tg_mat = matrix or RasterLayer defining the stochastic slip grid
#' @return list of wavenumbers in numerical space
#'
#' @export
sffm_get_numerical_wavenumbers<-function(tg_mat){

    tmp = dim(tg_mat)
    B1 = tmp[1]
    A1 = tmp[2]
    X2 = pmin(0:(B1-1), B1-(0:(B1-1)))/B1
    X1 = pmin(0:(A1-1), A1-(0:(A1-1)))/A1
    kx = matrix(X1,ncol=A1,nrow=B1, byrow=T)
    ky = matrix(X2, ncol=A1,nrow=B1)

    return(list(kx,ky))
}

###############################################################################
#' Synthetic finite fault model generator
#'
#' Make a random slip surface with the same dimensions as a provided template
#' RasterLayer or matrix. The random slip surface is generated using (by default)
#' the S_{NCF} algorithm in Davies et al., (2015), based on user-provided 
#' numerical corner wave-number parameters. The corner-wavenumber parameters
#' are specified in numerical space as (kcxN, kcyN) = (kcx, kcy) * (dx, dy) 
#' where kcxN, kcyN are the NUMERICAL corner wave numbers; (kcx, kcy) are the
#' physical corner wavenumbers (units of 1/distance), and (dx, dy) are the x/y
#' pixel spacing of the template raster or matrix. \cr
#' Further explanation, see the example: \cr
#'
#' @param reg_par vector passed to sffm_pars$spectral_amplitude_fun. First
#' two entries are kcxN, kcyN, as explained above and in the example
#' @param tg_mat is a 'template' raster, or matrix. The output slip distribution
#' will have these dimensions, and if tg_mat is a RasterLayer, it will have the
#' same properties (e.g. pixel size, spatial projection). 
#' @param sffm_pars list containing sffm configuration parameters. See
#' sffm_get_default_model_parameters()
#' @param sub_sample_size vector of length 2 with integers >= 1. tg_mat is 
#' refined to have rows/columns = sub_sample_size * (original rows/columns) before
#' the algorithm is applied. Values > 1 cause tg_mat to be sampled to finer
#' resolution before we simulate the sffm. The synthetic values are re-aggregated
#' into a matrix with the size of the original tg_mat prior to returning 
#' @return Output is the same class as tg_mat
#'
#' @references
#' Davies et al. (2015), 
#' Tsunami inundation from heterogeneous earthquake slip distributions:
#' Evaluation of synthetic source models, J. Geophys. Res. Solid Earth, 120,
#' 6431-6451, doi:10.1002/2015JB012272. \cr
#'
#' @export
#' @examples
#'
#' #
#' # Example simulating an SFFM
#' #
#' 
#' tg_mat = matrix(0, nrow=8, ncol=12)
#' tg_mat[3,5] = 1 # Fix peak slip location
#' xs = seq(0, 120, len=ncol(tg_mat)) # x-coordinates of tg_mat
#' ys = seq(0, 50, len=nrow(tg_mat)) # y-coordinates of tg_mat
#' dx = xs[2] - xs[1]
#' dy = ys[2] - ys[1]
#' # Make numerical corner wavenumbers c(kcxN, kcyN), corresponding to physical
#' # corner wavenumbers (1/50, 1/20) in the (x,y) directions respectively
#' reg_par = c(1/50 * dx, 1/20 * dy) 
#' random_slip_mat = sffm_simulate(reg_par, tg_mat)
#' 
#' ## Example plot
#' filled.contour(xs, ys, t(random_slip_mat), asp=1, nlevels=30, 
#'     color.palette=rainbow)
#' 
#' # Clipping should lead to patches of zero slip
#' stopifnot(min(random_slip_mat) == 0)
#'
#' #
#' # Example changing default parameters
#' #
#' 
#' new_sffm_parameters = sffm_get_default_model_parameters()
#' # Use absolute-value transformation to remove negative values, instead of clipping
#' new_sffm_parameters$negative_slip_removal_function <-function(x) abs(x)
#' # ... other changes could be made too ... #
#' random_slip_matB = sffm_simulate(reg_par, tg_mat, sffm_pars = new_sffm_parameters)
#' 
#' ## Example plot
#' filled.contour(xs, ys, t(random_slip_matB), asp=1, nlevels=30, 
#'     color.palette=rainbow)
#' 
#' # Should no longer be patches of zero values
#' stopifnot(min(random_slip_matB) > 0)
#'
sffm_simulate<-function(reg_par, tg_mat, sffm_pars = .sffm_default_model_parameters,
    sub_sample_size = c(1,1)){ 

    # Record random seed for reproducibility
    initial_seed = get_random_seed()

    # Optionally sample to a finer grid   
    if(any(sub_sample_size > 1)){
        #ensure sub_sample_size is integer
        stopifnot(all.equal(sub_sample_size, round(sub_sample_size)))
        stopifnot(length(sub_sample_size) == 2)

	# Store input values
	old_tg_mat = tg_mat
	old_reg_par = reg_par	

	# Make a finer tg_mat
        new_tg_mat = matrix(0, nrow=nrow(tg_mat)*sub_sample_size[1], 
            ncol=ncol(tg_mat)*sub_sample_size[2])

        # Give values to interpolated tg_mat 
        max_tg_orig = which(as.matrix(tg_mat) == max(as.matrix(tg_mat)), arr.ind=TRUE)
        new_max_tg = (max_tg_orig - 1) * sub_sample_size + ceiling(0.5*sub_sample_size)
        new_tg_mat[new_max_tg[1], new_max_tg[2]] = sum(as.matrix(tg_mat)) * prod(sub_sample_size)
	
        # Adjust reg_par to match new_tg_mat
        new_reg_par = reg_par
        new_reg_par[1:2] = reg_par[1:2] / sub_sample_size

        if(class(tg_mat) == 'RasterLayer'){
            new_tg_mat = raster(new_tg_max, xmn = extent(tg_mat)@xmin, 
                xmx = extent(tg_mat)@xmax, ymn = extent(tg_mat)@ymin,
                ymx = extent(tg_mat)@ymax)
        }

        tg_mat = new_tg_mat
        reg_par = new_reg_par
    }
 
    # Make wavenumber matrices
    tmp = sffm_get_numerical_wavenumbers(tg_mat)
    kx = tmp[[1]]
    ky = tmp[[2]]

    # Compute modelled amplitude spectrum
    model_fourierspec = sffm_pars$spectral_amplitude_fun(kx, ky, reg_par)

    # Note: if kx = 0,1/N,2/N, ... are the numerical wavenumbers,
    # they are equal to the ('physical' wavenumbers) x (delX)
    # In that case, 
    # reg_par[1] = kx_Numerical = kcx*delX where kcx is the physical 'corner
    # wavenumber'

    # Make the phase generating field
    r_noise = sffm_pars$random_phase_generating_field(length(tg_mat))
    dim(r_noise) = dim(tg_mat)[1:2]
    r_noise = r_noise*sign(sum(r_noise))

    if(TRUE){
        # Use the random variables to generate the phase
        # This preserves the model spectral relation
        # (prior to clipping / filtering / etc)
        # This is the technique used by Gallovic, Mai, ...
        fake_phase = Arg(fft(r_noise))
        tmp = model_fourierspec*(cos(fake_phase)+1i*sin(fake_phase))
    }else{
        # Use the random variables and multiply in fourier space
        # This is the technique used by Geist / Lavallee / ....
        # It can lead to a heavier skew in the slip (so truncation is normally
        # applied for stable distribution)
        tmp = model_fourierspec*fft(r_noise)
    }

    # Initial SFFM
    fake_data = Re(fft(tmp, inverse=T))/prod(dim(tg_mat))

    # Deal with negative values
    fake_data_clip = sffm_pars$negative_slip_removal_function(fake_data)

    # Recentre the slip so that location of maxima in fake_data_clip is the
    # same as that in tg_mat.
    # This must happen after negative value removal, in case e.g. abs() leads
    # to a new location for the maxima
    if(sffm_pars$recentre_slip){ 
        fake_data_clip = sffm_recentre_slip(fake_data_clip, tg_mat)
    }
    
    # Apply exponential slip decay. This must happen after recentering the slip
    if(sffm_pars$spatial_slip_decay %in% c('exponential', 'gaussian')){
        # Spatial filtering
        # Find distance from maxima
        maxInds = which(fake_data_clip == max(fake_data_clip), arr.ind=T)
        # Find distance from random point
        tmp = dim(fake_data_clip)[1:2]
        nc_fdc = tmp[2]
        nr_fdc = tmp[1]
        xMat = matrix(1:nc_fdc, ncol=nc_fdc, nrow=nr_fdc, byrow=TRUE)    
        yMat = matrix(1:nr_fdc, ncol=nc_fdc, nrow=nr_fdc, byrow=FALSE)    
        #
        # reg_par[1] = kcx*delx, so to get distance*kcx
        # [=distance/length_scale], we multiply the cell-count-distance by
        # reg_par[1].
        distMat2 = ( ((xMat-maxInds[2])*reg_par[1])**2 + 
                    ((yMat-maxInds[1])*reg_par[2])**2 )

        if(sffm_pars$spatial_slip_decay == 'exponential'){
            fake_data_clip = fake_data_clip*exp(-(distMat2**0.5))
        }else if(sffm_pars$spatial_slip_decay == 'gaussian'){
            fake_data_clip = fake_data_clip*exp(-(distMat2))
        }


    }else{
        # Just make sure sffm_pars has a valid value
        if(!(sffm_pars$spatial_slip_decay %in% c('none'))){
            stop('sffm_pars$spatial_slip_decay not recognized')
        }
    }

    # If required, reaggregate the slip values
    if(any(sub_sample_size > 1)){
        agg_fake_data_clip = as.matrix(old_tg_mat) * 0
        sr = seq(1, nrow(fake_data_clip), by=sub_sample_size[1])
        sc = seq(1, ncol(fake_data_clip), by=sub_sample_size[2])
        for(ic in 1:sub_sample_size[2]){
            for(ir in 1:sub_sample_size[1]){
                agg_fake_data_clip = agg_fake_data_clip + 
                    fake_data_clip[sr + (ir-1), sc + (ic-1)]
            }
        }
	# Redefine fake_data_clip and tg_mat to be the pre-interpolation
	# matrix, so that summation is correct later on
        tg_mat = old_tg_mat 
        fake_data_clip = agg_fake_data_clip
    }

    # Ensure final mean = data mean 
    if(class(tg_mat)=='RasterLayer'){
        fake_data_clip = fake_data_clip/sum(fake_data_clip)*sum(as.matrix(tg_mat))
    }else{
        fake_data_clip = fake_data_clip/sum(fake_data_clip)*sum(tg_mat)
    }

    if(class(tg_mat)=='RasterLayer' ){
        final_rast = raster(tg_mat)
        final_rast = setValues(final_rast, fake_data_clip)
    }else{
        final_rast = fake_data_clip
    }

    # Record the random number
    attr(final_rast, 'initial_seed') = initial_seed

    return(final_rast)        
}

#######################################################################
#'
#' Fix location of peak slip in sffm
#'
#' Move asperities closer to the centre of the rupture,
#' or to the location of asperities in another slip raster (tg) \cr
#' If tg is not NULL, then move the max of m1 to the same location as
#' the max of tg, by re-ordering rows and columns \cr
#' Otherwise, compute the max along each row -- then re-order the rows so that
#' the row with the smallest max is at the top of the rupture
#' (unless is is already at the top or bottom). \cr
#' Then do the same with the columns, trying to put the smallest max
#' col on the left (unless it is already at the left or right) \cr
#'
#' @param m1 rasterLayer or matrix containing a slip surface to be adjusted
#' @param tg rasterLayer or matrix containing a reference slip surface
#' @return recentred version of m1
#' @export
#'
sffm_recentre_slip<-function(m1, tg=NULL){

    if(!is.matrix(m1)){
        m1_mat=as.matrix(m1)
    }else{
        m1_mat = m1
    }

    if(is.null(tg)){
        # Find row with lowest max slip -- we want this on the bottom or top
        # edge
        row_to_bottom=which.min(apply(m1_mat,1,max))
        # Find col with lowest max slip -- we want this on the left or right
        # edge
        col_to_left=which.min(apply(m1_mat,2,max))
        nr=dim(m1_mat)[1]
        nc=dim(m1_mat)[2]
        if(!(row_to_bottom%in%c(1,nr))){
            #row_to_bottom is not on the bottom
            m1_mat=m1_mat[c( (row_to_bottom):nr, 1:(row_to_bottom-1)), ]
        }
        if(!(col_to_left%in%c(1,nc))){
            m1_mat=m1_mat[, c( (col_to_left):nc, 1:(col_to_left-1))]
        }

    }else{
        # Move the max to the location of the max of tg
        if(!is.matrix(tg)){
            tg_mat = as.matrix(tg)
        }else{
            tg_mat = tg
        }
        #new_row_max = which.max(apply(tg_mat,1,max)) 
        #new_col_max = which.max(apply(tg_mat,2,max)) 
        tmp = which(tg_mat == max(tg_mat), arr.ind=TRUE)
        new_row_max = tmp[1]
        new_col_max = tmp[2]

        #old_row_max = which.max(apply(m1_mat,1,max)) 
        #old_col_max = which.max(apply(m1_mat,2,max)) 
        tmp = which(m1_mat == max(m1_mat), arr.ind=TRUE)
        old_row_max = tmp[1]
        old_col_max = tmp[2]

        nr = dim(m1_mat)[1]
        nc = dim(m1_mat)[2]
        newRows = (((1:nr) - (new_row_max - old_row_max))%%nr)
        newRows[which(newRows==0)] = nr
        m1_mat = m1_mat[newRows, , drop=FALSE]
        newCols = (((1:nc) - (new_col_max-old_col_max))%%nc)
        newCols[which(newCols==0)] = nc 
        m1_mat = m1_mat[, newCols, drop=FALSE]
    }
    if(class(m1) == 'RasterLayer'){
        output = raster(m1)
        output = setValues(output,m1_mat)
    }else{
        # Work directly with matrices too
        output = m1_mat
    }
    return(output)
}


###############################################################################
#'
#' Goodness of fit for SFFM parameters
#'
#' Given regression parameters, compute a goodness-of-fit statistic of the
#' model with reg_par and data (tg_mat), based on Davies et al. (2015), Equation
#' 5. Most users would not call this routine directly (see sffm_fit_parameters 
#' for parameter estimation). \cr
#'
#' @param reg_par = vector of 2 proposed regression parameters (kcxN, kcyN) **in numerical
#' space**. A 3rd parameter may be accepted in some cases, depending on the
#' values of reg_par allowed in sffm_pars$spectral_amplitude_function. See
#' \code{?sffm_simulate} for more details on 'numerical space' and 'physical space'
#' @param tg_rast slip matrix or raster to compute the goodness of fit for
#' @param verbose TRUE/FALSE -- Verbose error messages
#' @param default_seed integer -- passed to set.seed for reproducible fitting
#'        with random fault generation (original .Random.seed is restored at the end)
#' @param NumRandSf integer. Number of slip distributions simulated to compute the
#'        goodness-of-fit of the model
#' @param sffm_pars environment containing configuration parameters
#' @return A goodness-of-fit measure -- minimising this will lead to the 'best'
#'         model fit
#' @export
#'
#' @references
#' Davies et al. (2015), 
#' Tsunami inundation from heterogeneous earthquake slip distributions:
#' Evaluation of synthetic source models, J. Geophys. Res. Solid Earth, 120,
#' 6431-6451, doi:10.1002/2015JB012272. \cr
#' 
sffm_slip_goodness_of_fit<-function(
    reg_par,
    tg_rast,
    verbose=FALSE,
    default_seed=1,   
    NumRandSf=200,
    sffm_pars = .sffm_default_model_parameters
    ){
    
    # Reproducible random seed 
    initial_seed = get_random_seed()
    # Reset the random seed to the original value when the function exits
    on.exit({.Random.seed <<- initial_seed})
    set.seed(default_seed) 

    #kcx=reg_par[1]
    #kcy=reg_par[2]
    #n=reg_par[3]

    # Treat invalid values -- we assume any reg_par[3] relates to 'alpha' in
    # the stable distribution
    if( (reg_par[1]<=0) | (reg_par[2]<=0)) return(9.0e+100)
    if( (length(reg_par)>2) && ((reg_par[3]<0)|reg_par[3] > 2)) return(9.0e+100)
   
    tg_mat = as.matrix(tg_rast)
    mean_tg_mat = mean(tg_mat)

    tmp = sffm_get_numerical_wavenumbers(tg_mat)
    kx = tmp[[1]]
    ky = tmp[[2]]


    # Store spectrum of all synthetic models here
    fake_f = array(NA,dim=c(nrow(kx), ncol(kx), NumRandSf))

    # Store fraction of fault covered with zeros here
    data_zero_frac = sum(tg_mat>0.0)/length(tg_mat)


    # Generate many random faults
    for(i in 1:NumRandSf){
        fake_data_cliprast = sffm_simulate(reg_par,tg_mat, 
            sffm_pars=sffm_pars)
        fake_f[,,i] = fake_data_cliprast
    }


    # Experiment with only keeping a zero-wavenumber component
    #trick_mat = tg_mat*0
    #trick_mat[,1] = 1
    #trick_mat[1,] = 1
    #trick_mat=1

    # Compute the 'goodness of fit' statistic
    spec_data = c(Mod(fft(tg_mat)))
    modstat = matrix(NA,ncol=length(spec_data),nrow=NumRandSf)
    for(i in 1:NumRandSf){
        modstat[i,] = c(Mod(fft(fake_f[,,i])))
    }
    mean_mod = colMeans(modstat)
    output_stat = mean((mean_mod-spec_data)**2)

    if(verbose){
        print(c(output_stat, reg_par, sum(output_stat), data_zero_frac))
    }

    return(output_stat)

}

#######################################################################
#'
#' Fit SFFM parameters
#'
#' Function to compute the optimal reg_par parameters for the stochastic slip
#' model. Uses the stochastic optimization method of  Davies et al. (2015)
#'
#' @param m1 matrix or RasterLayer with slip values
#' @param default_seed integer. Force the random seed value to this, in order
#' to make the fit reproducible
#' @param NumRandSf integer. number of synthetic slip simulations to use in
#' goodness-of-fit computation
#' @param sffm_pars list containing parameters passed to \code{sffm_simulate}. See
#' \code{?sffm_simulate} and \code{?sffm_get_default_model_pars} for more info.
#' @param reg_par_start vector of starting values for parameters \code{reg_par} in
#' \code{sffm_pars$spectral_amplitude_fun} which are optimized by the current
#' function. See \code{?sffm_simulate} and \code{?sffm_get_default_model_pars}
#' for more info.
#' @return an object from 'optim' with the fit
#' 
#' @export
#'
#' @references
#' Davies et al. (2015), 
#' Tsunami inundation from heterogeneous earthquake slip distributions:
#' Evaluation of synthetic source models, J. Geophys. Res. Solid Earth, 120,
#' 6431-6451, doi:10.1002/2015JB012272. \cr
#'
#' @examples
#' 
#' # Make a random slip model with known parameters, for testing the fit
#' tg_mat = matrix(0, nrow=8, ncol=12)
#' tg_mat[3,5] = 1 # Fix peak slip location
#' xs = seq(0, 120, len=ncol(tg_mat)) # x-coordinates of tg_mat
#' ys = seq(0, 50, len=nrow(tg_mat)) # y-coordinates of tg_mat
#' dx = xs[2] - xs[1]
#' dy = ys[2] - ys[1]
#' # Make numerical corner wavenumbers c(kcxN, kcyN), corresponding to physical
#' # corner wavenumbers 1/50, 1/20
#' reg_par = c(1/50 * dx, 1/20 * dy) 
#' # Make simulation reproducible 
#' set.seed(1)
#' random_slip_mat = sffm_simulate(reg_par, tg_mat)
#' 
#' # Now try to back-estimate the reg_par parameters for random_slip_mat, using
#' # the stochastic optimization method
#' fitted_sffm_par = sffm_fit_parameters(random_slip_mat)
#' # Check we are within 20% of the correct value (this will not always be true, but
#' # is for this case)
#' stopifnot(abs(reg_par - fitted_sffm_par$par)/reg_par < 0.2)
#' 
#' # The fit relies on stochastic simulations, so we get slightly different fits
#' # using different 'default_seed' values
#' fitted_sffm_parB = sffm_fit_parameters(random_slip_mat, default_seed=2)
#' stopifnot(abs(fitted_sffm_parB$par - fitted_sffm_par$par)/fitted_sffm_par$par < 0.05)
#' 
#' # #
#' # # Example of a more comprehensive test of the fitting algorithm.
#' # #
#' # # Test the fitting algorithm with simulation (run in parallel, since it can
#' # # take some time).
#' # # Idea: Generate many sffm with known kcx/kcy values, and then estimate the
#' # # latter values with sffm_fit_parameters. Ideally all fitted values would
#' # # equal the known values. In practice the fits are only approximate, and
#' # # we graphically illustrate this below 
#' # #
#' # library(parallel)
#' # nfits = 48 # Number of fits
#' # # Spread work over a number of shared memory cores -- my machine has 12,
#' # # and the following code takes a few minutes
#' # mc_cores = detectCores() 
#' # RNGkind("L'Ecuyer-CMRG") # ensure parallel random number generation
#' # random_slips = mclapply(1:nfits, f<-function(x) sffm_simulate(reg_par, tg_mat), 
#' #     mc.cores = mc_cores)
#' # random_slip_fits = mclapply(random_slips, f<-function(x) sffm_fit_parameters(x), 
#' #     mc.cores = mc_cores)
#' # 
#' # # Convert all the fitted parameters to a 2 column matrix
#' # all_fitted_par = matrix( unlist(lapply(random_slip_fits, f<-function(x) x$par)), 
#' #     ncol=2, byrow=TRUE)
#' # par(mfrow=c(1,2))
#' # hist(all_fitted_par[,1] / dx, main = 'Fitted kcx')
#' # abline(v=reg_par[1] / dx, col='red')
#' # hist(all_fitted_par[,2] / dy, main = 'Fitted kcy')
#' # abline(v=reg_par[2] / dy, col='red')
#' 
sffm_fit_parameters<-function(
    m1,
    default_seed = 1,   
    NumRandSf = 400,
    sffm_pars = .sffm_default_model_parameters,
    reg_par_start = c(2/dim(m1)[1], 2/dim(m1)[2])
    ){

    m1_matrix = as.matrix(m1)

    output=try(
           optim(par=reg_par_start, 
                 fn=sffm_slip_goodness_of_fit,
                 tg_rast=m1_matrix,
                 default_seed = default_seed,
                 NumRandSf = NumRandSf,
                 sffm_pars = sffm_pars,
                 #hessian=TRUE,
                 control=list(maxit=2000))
            )

    if(class(output)=='try-error'){
        #print(output)
        output=list(convergence=-999, output=output, 
            reg_par_start=reg_par_start, m1=m1)

    }else if(output$convergence!=0){

        # Run again to make convergence more likely
        # (e.g. for Nelder/Mead, there might not have been enough iterations)
        reg_par_start=output$par
        output=optim(par=reg_par_start, 
                     fn=sffm_slip_goodness_of_fit,
                     tg_rast=m1_matrix,
                     default_seed = default_seed,
                     NumRandSf = NumRandSf,
                     sffm_pars = sffm_pars,
                     #hessian=TRUE,
                     control=list(maxit=2000))
    }

    return(output)
}


