library(rptha)
all_Rdata = Sys.glob('*.Rdata')

for(RdataFile in all_Rdata){

    #load(RdataFile)
    attach(RdataFile)

    for(stat_measure in 1:3){

        # Get the 'best' stochastic and uniform events
        if(stat_measure == 1){
            # Time-domain similarity statistic
            stat_name = 'time'
            stoc_score = do.call(cbind, similar_s_time)
            unif_score = do.call(cbind, similar_u_time)
            vu_score = do.call(cbind, similar_vu_time)
        }else if(stat_measure == 2){
            # Spectral domain similarity statistic
            stat_name = 'spec'
            stoc_score = do.call(cbind, similar_s_spec)
            unif_score = do.call(cbind, similar_u_spec)
            vu_score = do.call(cbind, similar_vu_spec)
        }else if(stat_measure == 3){
            # Combined Time-domain and spectral domain statistic
            stat_name = 'spectime'
            stoc_score = 0.8*do.call(cbind, similar_s_time) + 0.2*do.call(cbind, similar_s_spec)
            unif_score = 0.8*do.call(cbind, similar_u_time) + 0.2*do.call(cbind, similar_u_spec)
            vu_score = 0.8*do.call(cbind, similar_vu_time) + 0.2*do.call(cbind, similar_vu_spec)
        }

        stoc_score_median = apply(stoc_score, 1, median)
        unif_score_median = apply(unif_score, 1, median)
        vu_score_median = apply(vu_score, 1, median)
        # Find the 'best' modelled events, using the 
        # minimum(median of the goodness-of-fit metric at all
        # gauges) as the definition of 'best'
        si = which.min(stoc_score_median)
        ui = which.min(unif_score_median)
        vui = which.min(vu_score_median)
        
        output_png = paste0(output_name_base, '_allmodel_goodness_fit_', stat_name, '_plot.png')
        png(output_png, width=13, height = 5, units='in', res=300)
        par(mfrow=c(1,3))
        hist(stoc_score_median, freq=FALSE, 
            main=paste0('Stochastic slip model median GoF: ', stat_name), 
            xlim=c(-0.1,1.1), col='red', density=20)
        hist(unif_score_median, freq=FALSE, 
            main=paste0('Uniform slip model median GoF: ', stat_name), 
            xlim=c(-0.1,1.1), col='blue', density=20)
        hist(vu_score_median, freq=FALSE, 
            main=paste0('Variable uniform slip model median GoF: ', stat_name), 
            xlim=c(-0.1,1.1), col='green', density=20)
        dev.off()

        #
        # Plot all the gauges and model results
        #
        n_gauges = length(uniform_slip_stats)
        ncol = (n_gauges > 7) + 1
        nrow = ceiling((n_gauges+1)/ncol)

        # Extract the gauge names from the files
        gauge_names = unlist(lapply(strsplit(basename(names(uniform_slip_stats)), split='_'), f<-function(x) substr(x[2], 1, 5)))

        output_png = paste0(output_name_base, '_gauges_best_fit_', stat_name, '_plot.png')
        png(output_png, width=7.5*ncol, height = max(nrow,5), units='in', res=300)
        par(mfrow=c(nrow, ncol))
        par(mar=c(0.3,1.5,0.5,1.5))
        par(omd=c(0.05, 0.95, 0,1))
        for(i in 1:n_gauges){

            uss = uniform_slip_stats[[i]][[ui]]
            sss = stochastic_slip_stats[[i]][[si]]
            svu = variable_uniform_slip_stats[[i]][[vui]]

            #xmin = min(c(uss$data_t, uss$model_t, sss$data_t, sss$model_t))
            xmin = min(c(sss$data_t, sss$model_t))
            xlim = c(xmin, xmin+4.0*3600)

            ymin = min(c(uss$data_s, uss$model_s, sss$data_s, sss$model_s, svu$model_s, svu$data_s))
            ymax = max(c(uss$data_s, uss$model_s, sss$data_s, sss$model_s, svu$model_s, svu$data_s))
            ylim = c(ymin, ymax)

            plot(uss$data_t, uss$data_s, t='o', col='black', pch=19, cex=0.4, xlim=xlim, ylim=ylim, axes=FALSE)
            points(sss$data_t, sss$data_s, t='o', col='black', pch=19, cex=0.4)
            points(svu$data_t, svu$data_s, t='o', col='black', pch=19, cex=0.4)
            points(uss$model_t - uss$model_time_offset, uss$model_s, t='l', col='blue', lwd=2, lty='longdash')
            points(sss$model_t - sss$model_time_offset, sss$model_s, t='l', col='red', lwd=1)
            points(svu$model_t - svu$model_time_offset, svu$model_s, t='l', col='green', lwd=1)
            abline(v = seq(xmin + 1, xmin + 7*3600, len=7), col='tan2', xpd=TRUE)
            abline(h=0, col='tan2')

            tmp = pretty(ylim) 
            axis(side=4, at = c(min(tmp[tmp>ylim[1]]), max(tmp[tmp<=ylim[2]])), las=1, line=-3, cex.axis=2.0)
            text(xlim[2]-1800, ylim[2]*0.8, labels=gauge_names[i], cex=3)
        }

        plot(c(0,1), c(0,1), axes=FALSE, col=0)
        legend('center', 
            c('Data', paste0('Stoch. slip (', si, ')'), paste0('Unif. slip (', ui, ')'), paste0('VarU. slip (', vui, ')')), 
            bty='n',
            lty=c('solid', 'solid', 'longdash', 'solid'),
            col=c('black', 'red', 'blue', 'green'), pch=c(19, NA, NA, NA), cex= 2.0 + (nrow > 2),
            pt.cex=1,
            ncol=2) 

        dev.off()

    }
    # Remove nearly all objects, except the ones controlling the outer loop
    rm(list=setdiff(ls(all=TRUE), c('RdataFile', 'all_Rdata')))

}
