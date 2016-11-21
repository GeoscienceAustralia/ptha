###############################################################################
#'
#' Code to fit/simulate spectral models of earthquake slip surfaces
#' , and wrappers to the Okada function
#'
#' FIXME: Integrate into rptha
#'
#' AUTHORS
#' Gareth Davies, Geoscience Australia 2013-2015, gareth.davies.ga.code@gmail.com
#'
#'

#' Convenience function to make the default_sffm_pars environment
get_default_sffm_pars<-function(){

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
    return(environment())
}

# Hold configuration variables in an environment
default_sffm_pars = get_default_sffm_pars()


###############################################################################
#'
#' Utility function to get the random seed -- and make it first if required
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
get_numerical_wavenumbers<-function(tg_mat){

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
#' Make a random slip surface from a template raster, given some regression
#' parameters (specified in NUMERICAL SPACE independent of the pixel size) and 
#' SFFM definitions in sffm_pars \cr
#' If kx = 0,1/N,2/N, ... are the numerical wavenumbers,
#' they are equal to the ('physical' wavenumbers) x (dx) \cr
#' where dx is the raster x cell size. \cr
#' In that case, if reg_par[1] and reg_par[2] are the numerical corner 
#' wavenumbers, then: \cr
#' reg_par[1] = kcx*dx  \cr
#' [where kcx is the physical 'corner wavenumber'], and similarly for
#' reg_par[2] with kcy/dy instead of kcx/dx
#'

#' @param reg_par parameter passed to sffm_pars$spectral_amplitude_fun. First
#' two entries are kcxN, kcyN, in NUMERICAL SPACE as explained above.
#' @param tg_mat is a 'template' raster, or matrix
#' @param sffm_pars environment containing configuration parameters. See
#' default_sffm_pars
#' @return Output is the same class as tg_mat
#'
#' @export
simulate_sffm<-function(reg_par, tg_mat, sffm_pars = default_sffm_pars){ 

    # Record random seed for reproducibility
    initial_seed = get_random_seed()
    
    # Make wavenumber matrices
    tmp = get_numerical_wavenumbers(tg_mat)
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
        fake_data_clip = recentre_slip(fake_data_clip, tg_mat)
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
#' Move asperities closer to the centre of the rupture,
#' or to the location of asperities in another slip raster (tg)
#'
#' If tg is not NULL, then move the max of m1 to the same location as
#' the max of tg, by re-ordering rows and columns
#' Otherwise, compute the max along each row -- then re-order the rows so that
#' the row with the smallest max is at the top of the rupture
#' (unless is is already at the top or bottom).
#' Then do the same with the columns, trying to put the smallest max
#' col on the left (unless it is already at the left or right)
#'
#' @param m1 rasterLayer or matrix containing a slip surface to be adjusted
#' @param tg rasterLayer or matrix containing a reference slip surface
#' @return recentred version of m1    
#'
recentre_slip<-function(m1, tg=NULL){

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
        m1_mat = m1_mat[newRows,]
        newCols = (((1:nc) - (new_col_max-old_col_max))%%nc)
        newCols[which(newCols==0)] = nc 
        m1_mat = m1_mat[,newCols]
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
#' Given regression parameters, compute a goodness-of-fit statistic of the
#' model with reg_par and data (tg_mat). This is useful for fitting
#' statistical models within optimization routines
#' 
#' The statistical model is 'largely' fit in 'numerical' space, not physical
#' space.
#' E.g. the smallest non-zero wavenumber resolvable by a DFT is 1/N,
#' and the nyquist wavenumber is 0.5
#' So be careful with units of kcx, kcy etc 
#'
#' @param reg_par = vector of 2 regression parameters  [kcx,kcy in numerical
#' space]. A 3rd parameter may be accepted in some cases, depending on the
#' values of reg_par allowed in sffm_pars$spectral_amplitude_function
#' 
#' @param tg_rast = slip grid raster to fit
#' @param verbose = TRUE/FALSE -- Verbose error messages
#' @param default_seed= integer -- passed to set.seed for reproducible fitting
#'        with random fault generation (original .Random.seed is restored at the end)
#' @param NumRandSf = Number of slip distributions simulated to compute the
#'        goodness-of-fit of the model
#' @param sffm_pars environment containing configuration parameters
#'
#' @return A goodness-of-fit measure -- minimising this will lead to the 'best'
#'         model fit
slip_goodness_of_fit<-function(
    reg_par,
    tg_rast,
    verbose=FALSE,
    default_seed=1,   
    NumRandSf=200,
    sffm_pars = default_sffm_pars
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

    tmp = get_numerical_wavenumbers(tg_mat)
    kx = tmp[[1]]
    ky = tmp[[2]]


    # Store spectrum of all synthetic models here
    fake_f = array(NA,dim=c(nrow(kx), ncol(kx), NumRandSf))

    # Store fraction of fault covered with zeros here
    data_zero_frac = sum(tg_mat>0.0)/length(tg_mat)


    # Generate many random faults
    for(i in 1:NumRandSf){
        fake_data_cliprast = simulate_sffm(reg_par,tg_mat, 
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
#' Function to compute the optimal reg_par parameters for the stochastic slip
#' model
#'
#' @param m1 matrix with slip values
#' @param default_seed force the random seed value
#' @param NumRandSf number of synthetic slip simulations to use in
#' goodness-of-fit computation
#' @param sffm_pars environment containing configuration parameters
#' @param reg_par_start starting values for parameters reg_par in
#' sffm_pars$spectral_amplitude_fun
#' @return an object from 'optim' with the fit
#' 
#' @export
fit_slip_parameters<-function(
    m1,
    default_seed=1,   
    NumRandSf=400,
    sffm_pars = default_sffm_pars,
    reg_par_start = c(2/dim(m1)[1], 2/dim(m1)[2])
    ){

    m1_matrix = as.matrix(m1)

    output=try(
           optim(par=reg_par_start, 
                 fn=slip_goodness_of_fit,
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
                     fn=slip_goodness_of_fit,
                     tg_rast=m1_matrix,
                     default_seed = default_seed,
                     NumRandSf = NumRandSf,
                     sffm_pars = sffm_pars,
                     #hessian=TRUE,
                     control=list(maxit=2000))
    }

    return(output)
}


