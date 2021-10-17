library(rptha)

ptha18_rate_curve_session = '../../../ptha_access/compute_rates_all_sources_session.RData'
if(!file.exists(ptha18_rate_curve_session)) stop('You need download compute_rates_all_sources_session.RData in the ptha_access directory. This can be done by sourcing "get_detailed_PTHA18_source_zone_info.R" in that folder. See the README in that folder for details')

load(ptha18_rate_curve_session)

#
# 'Nice-ish' source-zone specific plot
#

for(site in c('kermadectonga2', 'kermadectonga2_tonga', 'kermadectonga2_kermadec', 'kermadectonga2_hikurangi')){

    if(site == 'kermadectonga2'){

        site_title = 'Full source-zone unsegmented \n (50% weight on unsegmented)'
        title_cex = 1.8 # 2 #1.8
        lab_cex = 1.5
        axis_cex = 1.4
        legend_cex=1.4
        ylim = c(1.0e-04, 10)
        with_posterior_mean_CI=TRUE
        add_legend=TRUE

    }else if(site == 'kermadectonga2_tonga'){

        site_title = 'Tonga segment \n (50% weight on union of segments)'
        title_cex = 2.2 #1.8
        lab_cex = 1.5
        axis_cex = 1.4
        legend_cex=1.4
        ylim = c(1.0e-04, 10)
        with_posterior_mean_CI=TRUE
        add_legend=FALSE

    }else if(site == 'kermadectonga2_kermadec'){

        site_title = 'Kermadec segment \n (50% weight on union of segments)'
        title_cex = 2.2 # 2 #1.8
        lab_cex = 1.5
        axis_cex = 1.4
        legend_cex=1.4
        ylim = c(1.0e-04, 10)
        with_posterior_mean_CI=TRUE
        add_legend=FALSE

    }else if(site == 'kermadectonga2_hikurangi'){

        site_title = 'Hikurangi segment \n (50% weight on union of segments)'
        title_cex = 2.2 # 2 #1.8
        lab_cex = 1.5
        axis_cex = 1.4
        legend_cex=1.4
        ylim = c(1.0e-04, 10)
        with_posterior_mean_CI=TRUE
        add_legend=FALSE

    }

    se = source_envs[[site]]

    all_rate_curves = se$mw_rate_function(NA, return_all_logic_tree_branches=TRUE)

    outfile = paste0('rate_curve_plot_', site, '.png')
    png(outfile, width=7, height=5, units='in', res=200)

    options(scipen=5)

    plot(mw, se$mw_rate_function(mw), t='o', ylim=ylim, log='y', xlab="", 
        ylab='Exceedance Rate (events/year)', col='white', cex.lab=1.5, cex.axis=axis_cex)
    title(xlab=expression(M[w]), cex.lab=lab_cex*1.4, line=2)
    grid(col='orange')
    for(i in 1:nrow(all_rate_curves$all_rate_matrix)){
        points(all_rate_curves$Mw_seq, pmax(all_rate_curves$all_rate_matrix[i,], 0e-100), 
               t='l', col='grey', lwd=0.3)
    }
    add_log_axis_ticks(side=2)

    points(mw, se$mw_rate_function(mw), t='o', col='red', lwd=1, pch=19, cex=0.5)
    points(mw, pmax(se$mw_rate_function(mw, quantiles=0.975), 1.0e-100), t='l', col='purple', lty='dashed', lwd=2)
    points(mw, pmax(se$mw_rate_function(mw, quantiles=0.025), 1.0e-100), t='l', col='purple', lty='dashed', lwd=2)
    points(mw, pmax(se$mw_rate_function(mw, quantiles=0.84), 1.0e-100), t='l', col='blue', lty='dashed', lwd=2)
    points(mw, pmax(se$mw_rate_function(mw, quantiles=0.16), 1.0e-100), t='l', col='blue', lty='dashed', lwd=2)

    if(add_legend){
        legend('topright', 
            c('Logic-tree mean', 
              'Posterior 16/84 %', 'Posterior 2.5/97.5 %'),
            lty=c('solid', 'dashed', 'dashed'), 
            pch=c(19, NA, NA), bg=rgb(1, 1, 1, alpha=0.5), 
            cex=legend_cex, lwd = c(1, 2, 2), box.col=rgb(1,1,1,alpha=0.3),
            col=c('red', 'blue', 'purple'))

        legend('bottomleft', 
            c('Data', 'Logic-tree \nbranches'),
            lty=c('solid','solid'), pch=c(17, NA), pt.cex=c(1.5, NA),
            col=c('darkgreen', 'grey'), bg=rgb(1,1,1,alpha=0.3), 
            lwd = c(1, 2),
            cex=legend_cex, box.col=rgb(1,1,1,alpha=0.5))
    }

    title(paste0(site_title), cex.main=title_cex)

    # Add empirical Mw-vs-rate for GCMT data
    gcmt_data = se$gcmt_data
    if(nrow(gcmt_data) > 0){
        rnk = rank(gcmt_data$Mw)
        N = nrow(gcmt_data)

        # Empirical rate 
        aep = (N+1 - rnk) / gcmt_access$cmt_duration_years

        ordr = order(gcmt_data$Mw)

        points(gcmt_data$Mw[ordr], aep[ordr], col='darkgreen', pch=17, cex=1.5, t='o')

    }

    dev.off()
}
