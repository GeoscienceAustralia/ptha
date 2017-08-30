context('test_GR')

test_that('test_GR', {

    set.seed(1234)

    negloglik_GR<-function(b, x, mw_min){
        -sum(log(dGR(x, b, mw_min)))
    }

    negloglik_GR_binned<-function(b, x, mw_min, delta=0.1){
        -sum(log(pGR(x+delta/2, b, mw_min) - pGR(x-delta/2, b, mw_min)))
    }
    
    parlist<-list(
        list(b=1, mw_min=5),
        list(b=0.7, mw_min=6),
        list(b=1.3, mw_min=4))
    
    for(par in parlist){
        b = par$b
        mw_min = par$mw_min
        # Quick check that pGR/qGR are inverse functions
        p = runif(1000)
        expect_that(all(abs(range(
            pGR(
                qGR(p, b=b, mw_min=mw_min), b=b, mw_min=mw_min) - p)
            ) < 1.0e-12), is_true())
        #print('PASS')

        # Quick check that integral of dGR is equal to pGR
        qs =  seq(mw_min, 8.5, len=100)
        p_direct = pGR(qs, b=b, mw_min=mw_min)
        p_indirect = p_direct*0
        for(i in 1:100){
            p_indirect[i] = integrate(dGR, mw_min, qs[i], 
                b=b, mw_min=mw_min)$value
        }
        expect_that(all(abs(range(p_indirect - p_direct)) < 1.0e-12), is_true())
        #print('PASS')

        large_random_sample = rGR(1e+05, b=b, mw_min=mw_min)
        local_fit = optimize(negloglik_GR, c(0.6, 1.5), 
            x=large_random_sample, mw_min=mw_min)
        expect_that(abs(local_fit$minimum - b) < 0.02, is_true())
        #print('PASS')
    }

    #
    # Truncated GR tests
    #

    negloglik_truncGR<-function(b, x, mw_min, mw_max=Inf){
        -sum(log(dtruncGR(x, b, mw_min, mw_max)))
    }

    negloglik_truncGR_binned<-function(b, x, mw_min, mw_max=Inf, delta=0.1){
        # Make sure we do not pass mw < mw_min or mw > mw_max
        -sum(log(ptruncGR(x+delta/2, b, mw_min, mw_max) - 
            ptruncGR(x-delta/2, b, mw_min, mw_max)))
    }


    parlist=list(
        list(b=1, mw_min=5, mw_max=8.5),
        list(b=0.7, mw_min=5, mw_max=8.5),
        list(b=1.2, mw_min=5, mw_max=8.5),
        list(b=1.2, mw_min=3, mw_max=Inf),
        list(b=1., mw_min=3, mw_max=Inf)
    )

    for(par in parlist){

        b = par$b
        mw_min = par$mw_min
        mw_max = par$mw_max

        # Quick check that ptruncGR/qtruncGR are inverse functions
        p = runif(1000)
        expect_that(all(abs(range(
            ptruncGR(
                qtruncGR(p, b=b, mw_min=mw_min, mw_max=mw_max), 
                b=b, mw_min=mw_min, mw_max=mw_max) - p)
            ) < 1.0e-12), is_true())
        #print('PASS')

        # Quick check that integral of dtruncGR is equal to ptruncGR
        qs =  seq(5., 8.5, len=100)
        p_direct = ptruncGR(qs, b=b, mw_min=mw_min, mw_max=mw_max)
        p_indirect = p_direct*0
        for(i in 1:100){
            p_indirect[i] = integrate(dtruncGR, mw_min, qs[i], b=b, 
                mw_min=mw_min, mw_max=mw_max)$value
        }
        expect_that(all(abs(range(p_indirect - p_direct)) < 1.0e-12), is_true())
        #print('PASS')

        # Test that we can fit it.
        large_random_sample = rtruncGR(1e+05, b=b, mw_min=mw_min, mw_max=mw_max)
        local_fit = optimize(negloglik_truncGR, c(0.6, 1.5), x=large_random_sample, 
            mw_min=mw_min, mw_max=mw_max)
        expect_that(abs(local_fit$minimum - b) < 0.02, is_true())
        #print('PASS')
    }
})
