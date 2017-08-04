
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


#' Create a function to simulate random length/width/kcx/kcy for synthetic
#' finite fault models.
#'
#' The default parameter values are based on the Strasser et al. (2010)
#' length/width vs Mw scaling laws, and the kcx/kcy parameters from the S_{NCF}
#' stochastic slip model of Davies et al (2015). Beware that kcx/kcy are 'physical'
#' corner wavenumbers (units of distance^{-1}), whereas the function \code{simulate_sffm} requires these
#' to be transformed to dimensionless 'numerical' corner wavenumbers, by
#' division by the sffm cell size.
#'
#' @param log10kcx_regression_par vector of length 3 with the gradient,
#' intercept and residual standard deviation for the regression of log10(kcx) vs
#' earthquake-magnitude
#' @param log10kcy_regression_par vector of length 3 with the gradient,
#' intercept and residual standard deviation for the regression of log10(kcy) vs
#' earthquake-magnitude
#' @param cor_kcx_kcy_residual Pearson correlation between the log10(kcx) and
#' log10(kcy) regression residuals
#' @return a function f(Mw) which returns a vector or matrix with random
#' length, width, kcx, kcy values to be used to simulate an sffm with a given magnitude.
#' Units are km, km, km^{-1}, km^{-1} respectively.
#' @export
#' @examples
#' # Here we use default parameter values
#' lwkc = sffm_make_random_lwkc_function()
#' print(lwkc(7.5)) # Random parameters for Mw = 7.5
#' print(lwkc(7.5)) # More random parameters
#' print(lwkc(9.2)) # Larger L, W, and smaller kc
#' # The function can also take a vector, in which case it returns a matrix with one row for each Mw.
#' print(lwkc(c(7.5, 8.0, 8.5, 9.0)))
#'
sffm_make_random_lwkc_function<-function(
    log10kcx_regression_par=c(-0.54, 2.03, 0.22),
    log10kcy_regression_par=c(-0.41, 1.18, 0.19),
    cor_kcx_kcy_residual = 0.68){

    library(rptha)

    # Make a function to simulate the parameters, with regression residuals
    # having the desired correlations
    simulate_L_W_kcx_kcy<-function(Mw){

        # Get regression coefficients for log10(L), log10(W). Note the value
        # of 7.5 for Mw is arbitrary
        AWL_sigmas = Mw_2_rupture_size(7.5, detailed=TRUE)$log10_sigmas

        L_Mw = sapply(Mw, f<-function(x) Mw_2_rupture_size(x)['length'])
        W_Mw = sapply(Mw, f<-function(x) Mw_2_rupture_size(x)['width'])

        N = length(Mw)
        new_L = 10**( log10(L_Mw) + AWL_sigmas[3] * rnorm(N))

        new_W = 10**(log10(W_Mw) + AWL_sigmas[2] * rnorm(N))

        make_correlated_random_err<-function(N, correlation_coef){
            err_kcx = rnorm(N)
            err_kcy = correlation_coef * err_kcx + sqrt(1-correlation_coef**2)*rnorm(N)
            return(cbind(err_kcx, err_kcy))
        }
        err_kc = make_correlated_random_err(N, cor_kcx_kcy_residual)
        err_kcx = err_kc[,1]
        err_kcy = err_kc[,2]
        
        a = log10kcx_regression_par
        b = log10kcy_regression_par

        physical_corner_wavenumbers = 10**cbind(a[1]*Mw + a[2] + a[3]*err_kcx, 
            b[1]*Mw + b[2] + b[3]*err_kcy)

        output = cbind(new_L, new_W, physical_corner_wavenumbers)
        colnames(output) = c('L', 'W', 'kcx', 'kcy')
        rownames(output) = NULL

        return(output)
    }

    return(simulate_L_W_kcx_kcy)
}



#' Make a rectangle with a given size and target location on a grid
#'
#' Suppose we have a 2D grid (i.e. a logically rectangular set of cells),
#' and would like to define a rectangular sub-region consisting of L, W cells
#' with a given target x,y centre location. This problem arises when trying to define
#' regions of non-zero slip over unit sources for synthetic finite fault models. 
#' In many situations this is straightforward, but in some situations it is
#' ambiguous (e.g. if L and/or W is even, then no single grid cell is in the
#' 'centre'), or impossible (e.g. if a given target centre is too close to the
#' grid boundaries, it  may be impossible to make a rectangle with L/W having the target_centre).
#' This code solves the problem, taking care of the ambiguous situations (with randomness)
#' and the impossible situations (by putting the centre somewhere else), in a manner
#' which seems satisfactory for synthetic finite fault model simulation.
#'
#' @param grid_LW vector of length 2 giving the number of rows, columns on the 2D grid
#' @param num_LW vector of length 2 giving the desired number of rows,columns
#' on the rectangular sub-region
#' @param target_centre vector of length 2 giving the row,column indices on the
#' 2D grid where we would like the centre of the sub-region to be.
#' @return integer vector of length 4. The first 2 entries are the min/max row
#' indices of the sub-region. The 3rd and 4th entries are the min/max column
#' indices of the sub-region.
#'
#' @export
#' @examples
#' #
#' # Simple test case with deterministic answer
#' # 
#'
#' v1 = rectangle_on_grid(c(4, 10), c(3, 3), c(2, 5))
#' stopifnot(all(v1 == c(1, 3, 4, 6)))
#'
#' # Push a boundary
#' v1 = rectangle_on_grid(c(4, 10), c(3, 3), c(1, 5))
#' stopifnot(all(v1 == c(1, 3, 4, 6)))
#'
#' # Push two boundaries
#' v1 = rectangle_on_grid(c(4, 10), c(3, 3), c(1, 1))
#' stopifnot(all(v1 == c(1, 3, 1, 3)))
#'
#' v1 = rectangle_on_grid(c(4, 10), c(3, 3), c(4, 10))
#' stopifnot(all(v1 == c(2, 4, 8, 10)))
#'
#' #
#' #  Case with non-deterministic answer
#' #
#' v1 = rectangle_on_grid(c(4, 10), c(2, 2), c(2, 5))
#' stopifnot(v1[2]-v1[1]+1 == 2)
#' stopifnot(v1[4]-v1[3]+1 == 2)
#'
#' stopifnot(v1[1]%in%c(1,2))
#' stopifnot(v1[3]%in%c(4,5))
#'
#' # Deterministic because of the boundary
#' v1 = rectangle_on_grid(c(4, 10), c(2, 3), c(1, 5))
#' stopifnot(v1[2]-v1[1]+1 == 2)
#' stopifnot(v1[4]-v1[3]+1 == 3)
#'
#' stopifnot(v1[1]%in%c(1))
#' stopifnot(v1[3]%in%c(4))
#'
#' # Overly large rupture because of the boundary
#' v1 = rectangle_on_grid(c(4, 10), c(5, 3), c(1, 5))
#' stopifnot(v1[4]-v1[3]+1 == 3)
#'
#' stopifnot(v1[1:2] == c(1,4))
#' stopifnot(v1[3]%in%c(4))
#'
#' print('PASS')
#'
rectangle_on_grid<-function(grid_LW, num_LW, target_centre){        

        num_L = num_LW[1]
        num_W = num_LW[2]

        # Non-zero slip cover along-strike indices sL:eL, and down-dip indices
        # sW:eW. We try to make the peak slip location be the middle of the rupture,
        # but there are practical difficulties
        if(num_L%%2 == 0){
            # Peak slip location cannot be in the middle exactly, so we
            # randomly choose one location.
            sL = target_centre[1] - num_L/2 + sample(c(0,1), size=1)
        }else{
            sL = target_centre[1] - floor(num_L/2)
        }
        # If there are not enough unit sources, further constraints are needed
        sL = min(sL, grid_LW[1] - num_L + 1)
        sL = max(sL, 1)
        eL = min(sL + num_L - 1, grid_LW[1])

        if(num_W%%2 == 0){
            # Peak slip location cannot be in the middle exactly, so we
            # randomly choose one location.
            sW = target_centre[2] - num_W/2 + sample(c(0,1), size=1)
        }else{
            sW = target_centre[2] - floor(num_W/2)
        }
        # If there are not enough unit sources, further constraints are needed
        sW = min(sW, grid_LW[2] - num_W + 1)
        sW = max(sW, 1)
        eW = min(sW + num_W - 1, grid_LW[2])

        return(c(sL, eL, sW, eW))
}


#' Make synthetic finite fault models on a discretized source
#'
#' Currently only the S_{NCF} model from Davies et al. (2015) is implemented.
#' The corner wave-number parameters are generated stochastically based on
#' the regression relations provided therein. The peak slip location varies
#' about the target_location, within a region of approximately (L/2, W/2) where
#' L, W are the length/width of earthquakes of this magnitude according to 
#' Strasser's scaling relations 
#' 
#' @param discretized_source_statistics data.frame returned from e.g. 
#' \code{discretized_source_summary_statistics} or
#' \code{discretized_source_approximate_summary_statistics}
#' @param target_location numeric vector c(lon,lat) giving the desired peak slip
#' location. By default the peak slip will be distributed stochastically around this
#' point within a region of approximately (L/2, W/2), where L, W are the
#' length/width of uniform-slip earthquakes of this magnitude based on
#' Strasser's scaling relations. See vary_peak_slip_location
#' @param target_event_mw desired magnitude of the earthquake
#' @param num_events number of stochastic slip scenarios to produce
#' @param vary_peak_slip_location If TRUE, vary the peak slip location in a window
#' around the target location, as described above. If FALSE, use the same peak
#' slip location for every synthetic event
#' @param zero_low_slip_cells_fraction real between 0 and 1. If > 0,
#' then we zero all of the smallest slip values, which contribute < 'zero_low_slip_cells_fraction'
#' to the sum of all slip values on the rupture. A conservative value is 0, 
#' but small values (e.g. 0.05) might substantially reduce the number of unit-sources
#' involved in the rupture, with little distortion of the event.
#' @param sourcename Name of source (will be included in output list, can be useful for book-keeping)
#' @param sffm_sub_sample_size vector of 2 integers. The slip raster is
#' re-sampled to a higher resolution [sub_sample_size cells in the x/y directions
#' for each original cell] before generating the sffm, and then re-averaged before
#' returning the output.
#' @param mu Shear modulus (Pascals)
#' @param return_slip_raster logical. If TRUE, include the slip as a raster object in
#' the output list. If FALSE, set the latter to NULL on output.
#' @param uniform_slip logical. If TRUE, then the slip will be uniform on each rupture.
#' This enables the rupture size to vary while still having uniform slip. By default it
#' is FALSE, and stochastic slip is used.
#' @param expand_length_if_width_limited logical or character ('random'). This
#' governs what happens to the rupture length if the desired rupture width was
#' larger than the available width on the discrete source-zone. For instance, by
#' Strasser's scaling law, an Mw 9.5 has width which is log-normally
#' distributed around a median of 283 km. Many source-zones will not have enough
#' down-dip width to accommodate this. In that case, the width must be equal to
#' the maximum source-zone width, while the length can either remain the same
#' (if FALSE), or increase in proportion to the (width/desired_width) (if TRUE),
#' or vary randomly between the behaviours with probability 0.5 (if 'random').
#' @return A list with length = num_events. Each element of the list is a list
#' containing the entries slip_matrix, slip_raster, initial_moment, peak_slip_ind,
#' numerical_corner_wavenumbers, which should be self-explanatory if you are
#' familiary with \code{sffm_simulate}. Note the slip raster is returned for plotting
#' convenience only - one limitation is that all raster cells are the same size, whereas in
#' reality the unit-source sizes vary slightly.
#'
#' @export
#' @examples
#'
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'     layer='puysegur_downdip')
#' # Make discretized_source
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'     desired_subfault_length=50,
#'     desired_subfault_width=50,
#'     downdip_lines=puysegur_downdip)
#' # Get unit source summary statistics
#' puysegur_summary_statistics = discretized_source_summary_statistics(
#'     puysegur_discretized_source,
#'     approx_dx=5000,
#'     approx_dy=5000)
#' # Make stochastic slip events, with location/magnitude based on 2009/07/15
#' # puysegur earthquake
#' puysegur_sffm_event1 = sffm_make_events_on_discretized_source(
#'     puysegur_summary_statistics,
#'     target_location = c(166.568, -45.745),
#'     target_event_mw = 7.8,
#'     num_events = 10)
#'
#' # Plot the slip raster for the first 6 synthetic events
#' par(mfrow=c(3,2))
#' for(i in 1:6) plot(puysegur_sffm_event1[[i]]$slip_raster)

sffm_make_events_on_discretized_source<-function(
    discretized_source_statistics, 
    target_location,
    target_event_mw, 
    num_events = 1,
    vary_peak_slip_location=TRUE,
    zero_low_slip_cells_fraction=0.0,
    sourcename="",
    sffm_sub_sample_size = c(1,1),
    mu=3e+10,
    return_slip_raster=TRUE,
    uniform_slip = FALSE,
    expand_length_if_width_limited = 'random'){

    nx = max(discretized_source_statistics$alongstrike_number)
    ny = max(discretized_source_statistics$downdip_number)

    # Get 'typical' rupture dimensions, and allow the peak slip location
    # to be within L/2, W/2 of the CMT location
    rs = Mw_2_rupture_size(target_event_mw)
    if(vary_peak_slip_location){
        mean_us_width = mean(discretized_source_statistics$width)
        mean_us_length = mean(discretized_source_statistics$length)
        peak_slip_unit_source_window = c(
            ceiling(0.5*rs['width']/mean_us_width),
            ceiling(0.5*rs['length']/mean_us_length))
    }else{
        # Constant peak slip location
        peak_slip_unit_source_window = c(0, 0)
    }
        

    # Make the random kcx/kcy values appropriate for the magnitude. Simulatneously
    # make the length/width over which non-zero slip can occur
    random_LWkc_function = sffm_make_random_lwkc_function()
    LWkc = random_LWkc_function(rep(target_event_mw, length=num_events))
    physical_corner_wavenumbers = LWkc[,3:4]
    nonzero_slip_LW = LWkc[,1:2]

    # Find the along-strike/down-dip unit source index near the target_location
    target_unit_source_index = which.min(distHaversine(
        discretized_source_statistics[,c('lon_c', 'lat_c')], 
        matrix(target_location, ncol=2, nrow=length(nx), byrow=TRUE)))
    target_alongstrike = discretized_source_statistics$alongstrike_number[target_unit_source_index]
    target_downdip = discretized_source_statistics$downdip_number[target_unit_source_index]

    # Get average dx/dy for unit sources, where dx is along-strike and dy is
    # down-dip
    mean_dx = mean(discretized_source_statistics$length)
    mean_dy = sum(discretized_source_statistics$width * 
        discretized_source_statistics$length) / 
        sum(discretized_source_statistics$length)

    # Record full dx/dy in matrices
    dx = matrix(NA, ncol=nx, nrow=ny)
    dy = matrix(NA, ncol=nx, nrow=ny)
    for(i in 1:nrow(discretized_source_statistics)){
        nr = discretized_source_statistics$downdip_number[i] 
        nc = discretized_source_statistics$alongstrike_number[i]
        dx[nr, nc] = discretized_source_statistics$length[i]
        dy[nr, nc] = discretized_source_statistics$width[i]
    }

    desired_M0 = M0_2_Mw(target_event_mw, inverse=TRUE)

    #source_info = list()
    slip_generator_fun<-function(j){

        # Define the allowed ranges of the peak slip location
        target_downdip_range_min = max(1, target_downdip - peak_slip_unit_source_window[1])
        target_downdip_range_max = min(ny, target_downdip + peak_slip_unit_source_window[1])
        target_alongstrike_range_min = max(1, target_alongstrike - peak_slip_unit_source_window[2])
        target_alongstrike_range_max = min(nx, target_alongstrike + peak_slip_unit_source_window[2])

        # Randomly sample the peak slip location
        peak_slip_row = sample(target_downdip_range_min:target_downdip_range_max, size=1)
        peak_slip_col = sample(target_alongstrike_range_min:target_alongstrike_range_max, size=1)

        # Choose random 'numerical' corner wavenumbers
        numerical_corner_wavenumbers = physical_corner_wavenumbers[j,1:2] * 
            c(dx[peak_slip_row, peak_slip_col], 
              dy[peak_slip_row, peak_slip_col])

        # Find indices where we allow non-zero slip.
        slip_LW = nonzero_slip_LW[j,1:2]
        # Must include at least 1x1 unit source, though could have as many as
        # the entire source zone allows.
        num_L = round(slip_LW[1] / dx[peak_slip_row, peak_slip_col])
        num_L = min(ncol(dx), max(num_L, 1))
        num_W = round(slip_LW[2] / dy[peak_slip_row, peak_slip_col])
        num_W = min(nrow(dx), max(num_W, 1))

        #print('')
        #print(paste0('num_L: ', num_L, ' num_W: ', num_W))

        if(expand_length_if_width_limited != FALSE){
            # Check if the requested width is larger than the available width
            # If so, we might increase the length, depending on the value of
            # expand_length_if_width_limited
            max_available_width = nrow(dx) * dy[peak_slip_row, peak_slip_col]
            if(slip_LW[2] > max_available_width){
                # The width is less than desired

                width_deficit = num_W * dy[peak_slip_row, peak_slip_col] / slip_LW[2]

                #print(paste0('width_deficit: ', width_deficit))

                if(width_deficit < 1){
                    alternative_length = slip_LW[1] / width_deficit
                    alternative_num_L = round(alternative_length / dx[peak_slip_row, peak_slip_col])
                    alternative_num_L = min(ncol(dx), max(alternative_num_L, 1))
                   
                    # Number to support the 'random' option  
                    p = runif(1, min=0, max=1)
                    if( (expand_length_if_width_limited == TRUE) |
                        ((expand_length_if_width_limited == 'random') & (p>= 0.5))){
                        #print('expanding')
                        num_L = alternative_num_L
                    }
                }
            }
        }

        #print(paste0('num_L: ', num_L, ' num_W: ', num_W))

        # Non-zero slip cover along-strike indices sL:eL, and down-dip indices
        # sW:eW. We try to make the peak slip location be the middle of the rupture,
        # but there are practical difficulties
        bbox = rectangle_on_grid(dim(dx), c(num_W, num_L), c(peak_slip_row, peak_slip_col))
        sL = bbox[3]
        eL = bbox[4]
        sW = bbox[1]
        eW = bbox[2]

        slip_matrix = dx * 0

        if(!uniform_slip){
            # Ensure peak slip occurs in desired location
            template_slip_matrix = dx * 0
            template_slip_matrix[peak_slip_row, peak_slip_col] = 1

            repeater = TRUE
            while(repeater){
                # Simulate (possibly on sub-sampled grid) 
                # Note we only sample on the 'non-zero slip area', i.e. sW:eW, sL:eL
                slip_matrix[sW:eW, sL:eL] = sffm_simulate(
                    numerical_corner_wavenumbers, 
                    template_slip_matrix[sW:eW, sL:eL, drop=FALSE], 
                    sub_sample_size=sffm_sub_sample_size)
                # Ensure we have not accidently set part of the length/width to be fully zero
                if(any(slip_matrix[sW,sL:eL] > 0) & any(slip_matrix[eW, sL:eL] > 0) &
                    any(slip_matrix[sW:eW, sL] > 0) & any(slip_matrix[sW:eW, eL] > 0)){
                    repeater = FALSE
                }
            }

            rm(template_slip_matrix)

        }else{
            # Use uniform slip -- will be rescaled later
            slip_matrix[sW:eW, sL:eL] = 1
        }

        # There will probably be many small but nonzero slip values
        # Set some to zero, so for efficiency later
        threshold_level = zero_low_slip_cells_fraction 
        slip_sorted = sort(slip_matrix, decreasing=FALSE)
        cumulative_slip_sorted = cumsum(slip_sorted)
        if(threshold_level == 0){
            slip_threshold = 0
        }else{
            slip_threshold_ind = max(
                which(cumulative_slip_sorted < (threshold_level*max(cumulative_slip_sorted))))
            if(is.finite(slip_threshold_ind)){
                slip_threshold = slip_sorted[slip_threshold_ind]
            }else{
                slip_threshold = 0
            }
        }
        slip_matrix = slip_matrix * (slip_matrix > slip_threshold)
        rm(slip_sorted, cumulative_slip_sorted)

        # Ensure M0 is correct
        # We need slip * dx * dy * mu = M0
        #mu = 3e+10 # Now an input argument
        initial_moment = sum(slip_matrix * dx * dy * 1e+06 * mu)
        slip_matrix = slip_matrix/initial_moment * desired_M0
        stopifnot(abs(sum(slip_matrix * dx * dy * 1e+06 * mu) - desired_M0) < 
            (1.0e-06 * desired_M0))

        if(return_slip_raster){
            # Make a raster for nice output plots
            slip_raster = raster(slip_matrix, xmn=0, xmx=nx*mean_dx, ymx=0, 
                ymn=-ny*mean_dy)
        }else{
            slip_raster = NULL
        }

        output_list = list(
            slip_matrix = slip_matrix, 
            slip_raster = slip_raster, 
            initial_moment = initial_moment,
            peak_slip_ind = c(peak_slip_row, peak_slip_col),
            numerical_corner_wavenumbers = numerical_corner_wavenumbers,
            physical_corner_wavenumbers = physical_corner_wavenumbers,
            desired_LW = slip_LW[1:2],
            target_event_mw = target_event_mw,
            target_location = target_location,
            sourcename = sourcename)


        # Garbage-collect now-and-again
        if( j%%50 == 0 ) gc()

        return(output_list)
    }

    all_sffm_events = lapply(as.list(1:num_events), slip_generator_fun)

    return(all_sffm_events)
}


#' Convert the sffm events information to a table 
#'
#' Convert the sffm events to a data.frame to ease export to csv or netcdf. Note
#' the tabular format is not ideal, since the data is relatively unstructured.
#' For instance, it is necessary to store the indices with non-zero slip in a
#' character string.  The same issue arises for the slip values themselves.
#' 
#' @param all_sffm_events output of \code{sffm_make_events_on_discretized_source}
#' @param slip_significant_figures integer or NULL. If not NULL, then truncate slip
#' to this many significant figures when concatenating to string
#' @return data.frame with: columns event_index_string and event_slip_string, having
#' all the unit-source-indices and their slip values in a character string, separated
#' by '-' and '_' respectively; and other metadata about the event
#'
#' @export
#'
sffm_events_to_table<-function(all_sffm_events, slip_significant_figures=NULL){

    # Collapse the indices of unit sources with non-zero slip to a character
    # The string is like e.g. '20-22-23-37-', i.e. numbers separated by '-'
    # This is a clumsy way of integrating an array with irregular length into
    # each row of a data.frame
    event_index_string = unlist(lapply(all_sffm_events, 
        f<-function(x) paste0(which(c(x$slip_matrix) > 0), sep="-", collapse="")))

    # Collapse the slip on unit sources with non-zero slip to a character, possibly
    # with a reduction in the number of significant figures
    #
    # The string is like e.g. '0.234_1.23543_23.7_37_', i.e. numbers separated by '_'
    #
    # This is a clumsy way of integrating an array with irregular length into
    # each row of a data.frame
    if(is.null(slip_significant_figures)){
        event_slip_string = unlist(lapply(all_sffm_events, 
            f<-function(x){
                paste0(c(x$slip_matrix[x$slip_matrix > 0]), sep="_", collapse="")
            }
        ))
    }else{
        event_slip_string = unlist(lapply(all_sffm_events, 
            f<-function(x){
                paste0(c(signif(x$slip_matrix[x$slip_matrix > 0],slip_significant_figures)), 
                    sep="_", collapse="")
            }))
    }

    # Make the final output dataframe -- put various potentially useful metadata in too
    output_data = data.frame(
        event_index_string = event_index_string, 
        event_slip_string = event_slip_string, 
        Mw = unlist(lapply(all_sffm_events, f<-function(x) x$target_event_mw)),
        target_lon = unlist(lapply(all_sffm_events, f<-function(x) x$target_location[1])),
        target_lat = unlist(lapply(all_sffm_events, f<-function(x) x$target_location[2])),
        peak_slip_downdip_ind = unlist(lapply(all_sffm_events, f<-function(x) x$peak_slip_ind[1])),
        peak_slip_alongstrike_ind = unlist(lapply(all_sffm_events, f<-function(x) x$peak_slip_ind[2])),
        physical_corner_wavenumber_x= unlist(lapply(all_sffm_events, f<-function(x) x$physical_corner_wavenumbers[1])),
        physical_corner_wavenumber_y= unlist(lapply(all_sffm_events, f<-function(x) x$physical_corner_wavenumbers[2])),
        sourcename = unlist(lapply(all_sffm_events, f<-function(x) x$sourcename)),
        stringsAsFactors=FALSE)

    return(output_data)
}
