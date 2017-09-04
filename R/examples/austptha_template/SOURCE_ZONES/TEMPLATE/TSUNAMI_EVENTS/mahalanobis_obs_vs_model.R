require(nlshrink)
require(copula)

#
# Compute the squared mahalanobis distance of
# inverse-standard-normal-empirical-quantiles of rbind(model matrix and obs
# vector). The covariance matrix can be estimated with a shrinkage estimator
# (required if number of observations < dimension of problem, since otherwise 
# the covariance matrix is mathematically non-invertible), and optionally
# the obs vector can be ignored in estimating the covariance matrix
#
# Note there are statistical results for the distribution of this statistic. 
# However, they probably don't apply for our modified covariance matrix.
mahalanobis_obs_vs_models<-function(model_matrix, obs_vector, 
    use_nl_shrink = FALSE, non_parametric=TRUE, 
    ignore_obs_for_covariance_estimation = FALSE){


    if(is.null(dim(obs_vector))) dim(obs_vector) = c(1, length(obs_vector))

    full_matrix = rbind(model_matrix, obs_vector)

    p = ncol(full_matrix)
    
    # Apply transformations to data if required
    if(non_parametric){
        # Convert to values in (0-1) by empirical transform
        # Take inverse normal
        full_matrix_p_invnorm = qnorm(pobs(full_matrix), mean=0, sd=1)
    }else{
        # Assume the data is already unit-normal
        full_matrix_p_invnorm = full_matrix
    }
    
    # Create data for covariance matrix estimation
    if(!ignore_obs_for_covariance_estimation){

        data_for_cov = full_matrix_p_invnorm

    }else{

        if(non_parametric){
            data_for_cov = qnorm(pobs(model_matrix), mean=0, sd=1)
        }else{
            data_for_cov = model_matrix
        }

    }

    # Estimate covariance matrix
    if(!use_nl_shrink){
        covmat = cov(data_for_cov)
    }else{
        # Use nlshrink for covariance estimation, and suppress printing
        tmpvar = capture.output({covmat = nlshrink_cov(data_for_cov, k=1)})
    }

    dist = mahalanobis(full_matrix_p_invnorm, center=rep(0, p), cov=covmat)

    # Find whether the last value of 'dist' is large
    return(dist)
}

#' Suppose model_matrix contains rows with multivariate model realisations, and
#' obs_vector contains a row of observed multivariate data. 
#'
#' We are interested in the mahalanobis distance between the observations and
#' the data, after both are given an inverse_gaussian(empirical_cdf) transform.
#' The rank transformation is applied to make the comparison non-parametric. The
#' inverse Gaussian transform then makes the distribution to align with the idea
#' underlying the mahalanobis distance.
#'
#' A straightforward approach is to pool the model_matrix and observations, then
#' inverse_gaussian(empirical_cdf) transform, estimate the covariance matrix, and
#' then compare the MH distances for the observed data with the MH distances for
#' the model points. The idea is that if the observed data is unusual, then the
#' MH distance for the observation will be 'unusually large' compared with the
#' model points.
#'
#' However, we might be concerned that if the observation really is an outlier,
#' then its inclusion will make the covariance matrix estimator perform badly.
#' We could avoid this by not using the observed data in the covariance matrix
#' estimation -- however that could create bias in the comparison (if all the
#' other model_matrix points were included).
#'
#' An alternative is to:
#' foreach (row of model_matrix){
#'     *remove the row
#'     *compute mahalanobis distance for the observation AND the removed row,
#'        where the covariance matrix is only computed from the model_matrix
#'        with removed row.
#' }
#' Then compare the distribution of the observation MH distances, with the 
#'    removed-row MH distances.
#'
#'
bootstrap_mahalanobis_obs_vs_models<-function(model_matrix, obs_vector, 
    use_nl_shrink=FALSE, non_parametric=TRUE){

    dists_obs = rep(0, nrow(model_matrix))
    dists_leftout_model = rep(0, nrow(model_matrix))

    non_parametric = non_parametric
    use_nl_shrink = use_nl_shrink

    for(i in 1:length(dists_obs)){
        # Remove one observation from the model_matrix
        local_dat = model_matrix[-i,]
        # We will evaluate the removed observation against the remaining
        # observations.
        local_obs = model_matrix[i,]
        local_dists = mahalanobis_obs_vs_models(
            local_dat, 
            rbind(local_obs, obs_vector),
            use_nl_shrink=use_nl_shrink, 
            non_parametric=non_parametric,
            ignore_obs_for_covariance_estimation=TRUE)
        #
        dists_obs[i] = local_dists[length(local_dists)]
        dists_leftout_model[i] = local_dists[length(local_dists)-1]
    }
    return(list(
        dists_obs = dists_obs, 
        dists_leftout_model = dists_leftout_model)
        )
}

#' 
#' Code used to test it out
#'
.examples<-function(){
    library(MASS)
    # Make correlated normal random variables
    m1 = mvrnorm(1000, mu=c(0,0,0), 
        Sigma=matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), ncol=3))

    # Add an 'extremely typical' point!
    x = mahalanobis_obs_vs_models(m1, c(0,0,0))
    mean(x > x[length(x)])
    # Add a 'kind-of-unusual' point -- however, note it's an outlier in the 'usual'
    # way for this data [i.e. positive correlations between all variables]
    x = mahalanobis_obs_vs_models(m1, c(2,2,2))
    # Unusual? A bit, but not too much -- we expect positive correlations in
    # extremes
    mean(x > x[length(x)])
    # Consider a more unusual point -- it has high values of 2 var, but a very low
    # value of the others -- so the pattern does not follow the correlations in the
    # data.
    x = mahalanobis_obs_vs_models(m1, c(2,2,-2))
    # Unusual? Quite a bit!
    mean(x > x[length(x)])


    #
    # High dimensional example -- more dimensions than observations. Common situation
    # in e.g. image processing, microarrays, other 'big data' applications.
    #
    N = 20
    p = 28
    mu=rep(0, p)
    Sigma2=0.5 + diag(0.5, nrow=p, ncol=p) # Covariance matrix -- well correlated variables.
    m1 = mvrnorm(N, mu=mu, Sigma=Sigma2)
    # Test a 'very usual' point. Here, we have to use shrinkage methods to
    # get an invertible variance-covariance matrix
    x = mahalanobis_obs_vs_models(m1, rep(0, p), use_nl_shrink=TRUE)
    mean(x > x[N+1])
    # Test a seriously unusual point
    x = mahalanobis_obs_vs_models(m1, rep(c(-2,2), length.out=p), use_nl_shrink=TRUE, 
        ignore_obs_for_covariance_estimation=TRUE)
    mean(x > x[N+1])
    #
    # An 'outlier point' can mess up the covariance matrix estimate. To work around
    # that, we can compute the distance for each point, with the covariance matrix
    # estimated by ignoring itself.
    #
    #
    # Something with clear outliers
    #
    x_boot = bootstrap_mahalanobis_obs_vs_models(m1, 
        rep(c(-2,-2,-2,2),length.out=p), use_nl_shrink=TRUE)

    # Interestingly enough, a sample with all values at -2 is not 'unusual' for
    # this distribution -- it is actually 'probably more common' than any other
    # point we randomly sample
    x_boot = bootstrap_mahalanobis_obs_vs_models(m1, rep(c(-2),length.out=p), 
        use_nl_shrink=TRUE)
    # We can confirm that in a different way, by comparing the multivariate normal
    # density of points in m1 with the point [-2,-2,...]. Turns out the latter has
    # fairly high density, given our correlation matrix, as compared with the density 
    # of points in m1
    library(mvtnorm)
    dm1 = dmvnorm(m1, mean=mu, sigma=Sigma2)
    d_neg_2 = dmvnorm(rep(-2,length=p), mean=mu, sigma=Sigma2)


    #
    # High dimensional example with 'perfectly correlated' samples
    # 
    m1 = matrix(1:N, nrow=N, ncol=p)
    # This observation is extreme compared with the others, but is in keeping
    # with the idea of perfectly correlated ranks
    test_obs = rep(N+1, length=p)
    ## This fails! The covariance matrix is constant
    #x = mahalanobis_obs_vs_models(m1, test_obs, use_nl_shrink=TRUE)
}


if(FALSE){

    load('gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata')

    stochastic_S = matrix(NA, ncol=length(stochastic_slip_stats), 
        nrow=length(stochastic_slip_stats[[1]]))
    uniform_S = matrix(NA, ncol=length(uniform_slip_stats), 
        nrow=length(uniform_slip_stats[[1]]))
    data_S = rep(NA, length=length(uniform_slip_stats))

    for(j in 1:28){
        # Extract stage range data
        for(i in 1:170){
            stochastic_S[i,j] = diff(stochastic_slip_stats[[j]][[i]]$model_range)
        }
        for(i in 1:17){
            uniform_S[i,j] = diff(uniform_slip_stats[[j]][[i]]$model_range)
        }
        data_S[j] = diff(stochastic_slip_stats[[j]][[1]]$data_range)
    }

    # Compare model distance with data distance
    mh_S = mahalanobis_obs_vs_models(stochastic_S, data_S, use_nl_shrink=TRUE, 
        ignore_obs_for_covariance_estimation=TRUE)
    mh_U = mahalanobis_obs_vs_models(uniform_S, data_S, use_nl_shrink=TRUE, 
        ignore_obs_for_covariance_estimation=TRUE)

    # Compare model distance with data distance, when data points are left out
    # of model 1-by-1
    mh_S_boot = bootstrap_mahalanobis_obs_vs_models(stochastic_S, data_S, 
        use_nl_shrink=TRUE)

    mh_U_boot = bootstrap_mahalanobis_obs_vs_models(uniform_S, data_S, 
        use_nl_shrink=TRUE)

}
