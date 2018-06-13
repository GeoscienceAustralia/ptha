#
# Single station stage exceedance rate computation
#
# This is actually a slow computational method ('quick' means 'quick to code'),
# but is useful to cross check the other results
#


quick_source_deagg<-function(lon, lat){
    #lon = 151.41
    #lat = -34.08

    tsunami_files = Sys.glob(
        '../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_*.nc')
    earthquake_only_files = gsub('_tsunami', '', tsunami_files)
    gauge_file = 'tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc'

    library(rptha)

    # Get hazard points -- faster to not use the '_tsunami' file
    fid = nc_open(gauge_file, readunlim=FALSE)
    n = length(fid$dim$station$vals)
    hp = data.frame(
        lon     = ncvar_get(fid, 'lon'),
        lat     = ncvar_get(fid, 'lat'),
        elev    = ncvar_get(fid, 'elev'),
        gaugeID = ncvar_get(fid, 'gaugeID')
        )
    nc_close(fid)

    # Find index of point nearest to lon/lat
    ni = lonlat_nearest_neighbours(cbind(lon, lat), cbind(hp$lon, hp$lat))

    # Get stage and rates, for each source

    stage_rate = vector(mode='list', length=length(tsunami_files))
    names(stage_rate) = basename(dirname(dirname(tsunami_files)))
    for(i in 1:length(tsunami_files)){
        print(basename(tsunami_files[i]))
        fid = nc_open(tsunami_files[i], readunlim=FALSE)
        fid_rates = nc_open(earthquake_only_files[i], readunlim=FALSE)

        event_rate = ncvar_get(fid_rates, 'variable_mu_rate_annual')
        event_rate_upper = ncvar_get(fid_rates, 'variable_mu_rate_annual_upper_ci')
        event_rate_lower = ncvar_get(fid_rates, 'variable_mu_rate_annual_lower_ci')
        event_Mw = round(ncvar_get(fid_rates, 'Mw'), 3) # Deal with floating point imperfections in netcdf
        event_Mw_vary_mu = round(ncvar_get(fid_rates, 'variable_mu_Mw'))
        peak_stage = ncvar_get(fid, 'max_stage', start=c(1, ni), count=c(-1,1))
        site = rep(basename(dirname(dirname(tsunami_files[i]))), length=length(event_rate))
        row_index = 1:length(event_rate)

        stage_rate[[i]] = data.frame(
            event_rate = event_rate,
            event_rate_upper = event_rate_upper,
            event_rate_lower = event_rate_lower,
            peak_stage = peak_stage,
            site = site,
            row_index=row_index,
            event_Mw = event_Mw,
            event_Mw_vary_mu = event_Mw_vary_mu)

        nc_close(fid)
        nc_close(fid_rates)
    }

    # Back-calculate the stage-vs-rate curves
    stage_rate_all = do.call(rbind, stage_rate)

    odr = rev(order(stage_rate_all$peak_stage))

    stg = stage_rate_all$peak_stage[odr]
    er = cumsum(stage_rate_all$event_rate[odr])
    er_up = cumsum(stage_rate_all$event_rate_upper[odr])
    er_lo = cumsum(stage_rate_all$event_rate_lower[odr])
    
    # Convergence check
    s1 = seq(1, length(odr), by=2)
    s2 = seq(2, length(odr), by=2)
    stg_conv_1 = stage_rate_all$peak_stage[odr[s1]]
    stg_conv_2 = stage_rate_all$peak_stage[odr[s2]]
    er_conv_1 = cumsum(stage_rate_all$event_rate[odr[s1]]) * 2
    er_conv_2 = cumsum(stage_rate_all$event_rate[odr[s2]]) * 2

    # We will compare with the values in the file
    fid = nc_open('tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc', 
        readunlim=FALSE)
    stages = fid$dim$stage$vals
    ers = ncvar_get(fid, 'variable_mu_stochastic_slip_rate'            , start=c(1,ni), count=c(-1,1))
    ers_up = ncvar_get(fid, 'variable_mu_stochastic_slip_rate_upper_ci', start=c(1,ni), count=c(-1,1))
    ers_lo = ncvar_get(fid, 'variable_mu_stochastic_slip_rate_lower_ci', start=c(1,ni), count=c(-1,1))
    nc_close(fid)


    #
    # Plot the data
    #
    plot_stage_vs_rate<-function(){
        xmax = max(stg*(er>0), na.rm=TRUE)
        ylim = c(1e-05, max(max(er), 1))
        plot(stg, er, t='l', log='xy', 
            xlim=c(min(0.01, max(xmax-0.005, 1.0e-05)), xmax),
            ylim=ylim,
            xlab='Stage (m)', ylab= 'Exceedance rate (events/year)')
        grid(col='brown')
        points(stg, er_up, t='l', col='red')
        points(stg, er_lo, t='l', col='red')
        title(paste0('Stage-vs-exceedance-rate @ (lon=', round(hp$lon[ni],3), ', lat=', 
            round(hp$lat[ni], 2), ', elev=', round(hp$elev[ni],2), ', ID=', round(hp$gaugeID[ni], 2),
            ') \n (Lines and points should overlap everywhere -- or there is a database mistake!)'))

        points(stages, ers, pch=19, cex=1.0, col='brown')
        points(stages, ers_up, pch=19, cex=1.0, col='pink')
        points(stages, ers_lo, pch=19, cex=1.0, col='pink')

        legend('topright', 
            c('Peak stage exceedance rate (mean over all logic-tree branches)',
              '95% credible interval'),
            col=c('brown', 'pink'),
            pch=c(19, 19), bg='white')
   
        # Add convergence check information 
        
        plot(stg, er, t='l', log='xy', 
            xlim=c(min(0.01, max(xmax-0.005, 1.0e-05)), xmax),
            ylim=ylim,
            xlab='Stage (m)', ylab= 'Exceedance rate (events/year)')
        grid(col='brown')
        points(stg_conv_1, er_conv_1, t='l', col='orange')
        points(stg_conv_2, er_conv_2, t='l', col='red')
        title(paste0(
            'The convergence_check1 and convergence_check2 curves are made using half the data each. They should agree fairly well \n except for rare events. Where they disagree significantly, the rates should be considered unreliable (i.e. avoid use)'))

        legend('topright', 
            c('Peak stage exceedance rate (mean over all logic-tree branches)',
              'convergence_check1', 'convergence_check2'),
            col=c('black', 'orange', 'red'),
            lty=c(1, 1, 1), bg='white')
    }
    

    #
    # Function to examine the distribution of earthquake magnitudes
    #
    peak_stage_magnitude_summary<-function(stage_threshold, source_zone){

        is_site = as.character(stage_rate_all$site) == source_zone

        mw_max = max(stage_rate_all$event_Mw[is_site & stage_rate_all$event_rate > 0])
        mw_min = min(stage_rate_all$event_Mw[is_site & stage_rate_all$event_rate > 0])

        if(is.finite(mw_max) & is.finite(mw_min)){
            k = which( (stage_rate_all$peak_stage > stage_threshold) & is_site & (stage_rate_all$event_rate>0) &
                stage_rate_all$event_Mw >= mw_min & stage_rate_all$event_Mw <= mw_max)
        }else{
            k = c()
        }

        if(length(k) == 0){
            output = NA
        }else{
            output = aggregate(stage_rate_all$event_rate[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x), length(x)))
            output_upper = aggregate(stage_rate_all$event_rate_upper[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x)))
            output_lower = aggregate(stage_rate_all$event_rate_lower[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x)))


            output = data.frame(Mw=output[,1], rate_exceeding=output$x[,1], 
                n=output$x[,2], rate_exceeding_lower=output_lower$x, 
                rate_exceeding_upper=output_upper$x)
        
            # Compute overall numbers of Mw events
            fracs = rep(NA, length(output[,1]))
            for(i in 1:length(output[,1])){
                n = sum(is_site & (stage_rate_all$event_Mw == output[i,1]))
                fracs[i] = output[i,3]/n
            }
        
            output = cbind(output, data.frame(fraction_events = fracs))
        }

        return(output)
    }

    #
    # Plot
    #
    plot_deaggregation_summary<-function(stage_threshold){

        k = which( (stage_rate_all$peak_stage > stage_threshold) & (stage_rate_all$event_rate > 0))
        if(length(k) == 0){
            plot(c(0, 1), c(0, 1), 
                main=paste0('No events exceeding stage_threshold = ', stage_threshold, 'm'))
        }else{

            par(mfrow=c(2,2))

            rate_by_source = aggregate(stage_rate_all$event_rate[k], 
                by=list(source_zone=as.character(stage_rate_all$site[k])), sum)

            # Sort from highest to lowest
            m1 = order(rate_by_source$x, decreasing=TRUE)
            if(length(m1) > 10){
                # Plot at most 10 source-zones
                rate_by_source = rate_by_source[m1[1:10],]
                m1 = order(rate_by_source$x, decreasing=TRUE)
            }

            # Color the top 3 source-zones differently
            colz = c('grey', 'red')[ 1 + (( (1:length(m1)) %in% m1[1:3]) & (rate_by_source$x > 0))]

            #dotchart(rate_by_source$x, labels=rate_by_source[,1], 
            #    main=paste0('All source-zones: Rate of events with \n stage exceeding ', stage_threshold),
            #    xlab='Rate (events/year)', 
            #    color=colz)
            oldmar = par('mar')
            par('mar' = oldmar + c(0, 8, 0, 0)) # new margins to allow source-zone names to fit on plot
            barplot(rate_by_source$x[m1], 
                names.arg=as.character(rate_by_source[m1,1]),
                col=colz[m1], density=100, horiz=TRUE, las=1, 
                xlab='Rate (events/year)', 
                main=paste0('Top ', length(m1), 
                    ' source-zones: Rate of events with \n peak_stage > ', stage_threshold, 'm', 
                    ' @ (', round(hp$lon[ni],3), ', ', round(hp$lat[ni],3), ')')
                )
            par('mar' = oldmar) # Back to old margins

            # Also plot rate-vs-Mw for the 3 largest contributors 
            for(i in 1:min(length(m1), 3)){
                # Name of source-zone
                sz = rate_by_source[m1[i], 1]
                rate_by_Mw = peak_stage_magnitude_summary(stage_threshold, sz) 
                dotchart(rate_by_Mw$rate_exceeding,
                    labels=paste0(rate_by_Mw$Mw, ' (', round(rate_by_Mw$fraction_events*100, 1), ')'),
                    xlab='Rate > stage_threshold (events/year) with 95% CI', ylab='', 
                    pch=19, xlim=c(0, max(rate_by_Mw$rate_exceeding_upper)))
                points(rate_by_Mw$rate_exceeding_upper, 1:nrow(rate_by_Mw), col='red')
                points(rate_by_Mw$rate_exceeding_lower, 1:nrow(rate_by_Mw), col='red')
                mtext(side=2, 'Magnitude, fixed mu, (% of scenarios > stage_threshold)', line=2.3, cex=0.5)
                title(paste0(sz, ': Rate of events of each magnitude (fixed_mu) \n that are > stage_threshold (', signif(stage_threshold, 2), 
                    ')'))
            }
        }
    }

    pdf(paste0('station_summary_', lon, '_', lat, '.pdf'), width=12, height=7)
    plot_stage_vs_rate()
    plot_deaggregation_summary(0.3)
    plot_deaggregation_summary(1.0)
    plot_deaggregation_summary(2.0)
    dev.off()

}
