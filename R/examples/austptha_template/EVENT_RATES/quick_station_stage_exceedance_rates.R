#
# Single station stage exceedance rate computation
#
# This is actually a slow computational method ('quick' means 'quick to code'),
# but is useful to cross check the other results
#

.preamble_title = paste0('2018 Australian Probabilistic Tsunami Hazard Assessment single station summary')
.preamble_text = paste0('\n', 
                       "This file gives a summary of the Geoscience Australia's 2018 PTHA results at a single station. See the README on:\n",
                       '    https://github.com/GeoscienceAustralia/ptha/tree/master/ptha_access \n',
                       'for further information on accessing the results and associated reports. Users should understand this material before using the results.\n',
                       '\n',
                       'The plots are: \n',
                       '\n',
                       '    1) A hazard curve, containing the peak-tsunami-stage vs the exceedance-rate at the specified hazard point location, \n',
                       '       with 95% credible intervals describing the uncertainty. Note the peak-tsunami-stage is the defined as the maximum water level\n ',
                       '       attained by the tsunami (above mean-sea-level, ignorning tides)\n',
                       '\n',
                       '    2) A "convergence check" of the above hazard curve. The PTHA hazard result rely on simulating a large number of random\n',
                       '       earthquake-tsunami scenarios. At sufficiently rare exceedance rates (less frequent than some particular value, e.g. 1/10000 years)\n', 
                       '       there will be few random model scenarios, so exceedance rates become sensitive to the "random details" of those scenarios (i.e. less\n', 
                       '       reliable). To help users judge when this happens we compare 2 hazard curves, each derived from half the scenarios. \n',
                       '\n',
                       '    3-4-5) Information on which source-zones dominate the hazard (i.e. a hazard-deaggregation plot) for peak-stage thresholds of 0.3m,\n', 
                       '       1m, and 2m. In each case, for the top 3 source-zones we show rates separated by magnitude ("constant shear modulus magnitude"),\n',
                       '       to highlight the model scenarios most likely to cause tsunami above the peak-stage threshold.\n',
                       '       These plots are useful for determining which sources and scenarios to focus on, when conducting for tsunami hazard studies.\n',
                       '       Note that if the true shear modulus varies with depth, then the true magnitudes may differ from the "constant shear modulus \n',
                       '       magnitude" reported on this figure. Our event files also report on the "variable_mu_Mw" which is the magnitude assuming a \n',
                       '       depth-dependent shear modulus on subduction zones, following Bilek and Lay (1999). See the report for further information. \n',
                       '\n',
                       '\n',
                       'A copy of the script used to make this plot can be found at:\n',
                       'https://github.com/GeoscienceAustralia/ptha/blob/master/R/examples/austptha_template/EVENT_RATES/quick_station_stage_exceedance_rates.R \n',
                       'Other codes used to develop the PTHA can be found in the same git repository. \n',
                       '\n',
                       '\n',
                       'References\n',
                       '\n',
                       '    Bilek, S.L. and Lay, T. (1999) Rigidity variations with depth along interplate megathrust faults in subduction zones. Nature, 400, 443-446 \n'
                        
	)

.plot_preamble<-function(){

    plot(c(0,1), c(0,1), col='white', frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    title(.preamble_title, cex=1.5)
    text(0, 0.7, .preamble_text, pos=4)

}


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

    # Reduce the size of some lines that occur in the plot
    # This reduces the file size, important if we distribute many
    stages_interp = approx(stages, n=5*length(stages))$y
    stg_small = stages_interp
    er_small = approx(stg, er, xout=stages_interp)$y
    er_up_small = approx(stg, er_up, xout=stages_interp)$y
    er_lo_small = approx(stg, er_lo, xout=stages_interp)$y
    stg_conv_1_small = stages_interp
    stg_conv_2_small = stages_interp
    er_conv_1_small = approx(stg_conv_1, er_conv_1, xout=stages_interp)$y
    er_conv_2_small = approx(stg_conv_2, er_conv_2, xout=stages_interp)$y

    #
    # Plot the data
    #
    plot_stage_vs_rate<-function(){
        xmax = max(stg*(er>0), na.rm=TRUE)
        ylim = c(1e-05, max(max(er), 1))
        plot(stg_small, er_small, t='l', log='xy', 
            xlim=c(min(0.02, max(xmax-0.005, 1.0e-05)), xmax),
            ylim=ylim,
            xlab='Stage (m)', ylab= 'Exceedance rate (events/year)')
        abline(h=10**(seq(-5, 0)), lty='dotted', col='brown')
        abline(v=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), lty='dotted', col='brown')
        points(stg_small, er_up_small, t='l', col='red')
        points(stg_small, er_lo_small, t='l', col='red')
        title(paste0('Stage-vs-exceedance-rate @ (lon=', round(hp$lon[ni],3), ', lat=', 
            round(hp$lat[ni], 2), ', elev=', round(hp$elev[ni],2), ', ID=', round(hp$gaugeID[ni], 2),
            ') \n (Lines and points should overlap -- if they differ, contact the PTHA maintainer!)'))

        points(stages, ers, pch=19, cex=1.0, col='brown')
        points(stages, ers_up, pch=19, cex=1.0, col='pink')
        points(stages, ers_lo, pch=19, cex=1.0, col='pink')

        legend('topright', 
            c('Peak stage exceedance rate (mean over all logic-tree branches)',
              '95% credible interval'),
            col=c('brown', 'pink'),
            pch=c(19, 19), bg='white')
   
        # Add convergence check information 
        
        plot(stg_small, er_small, t='l', log='xy', 
            xlim=c(min(0.02, max(xmax-0.005, 1.0e-05)), xmax),
            ylim=ylim,
            xlab='Stage (m)', ylab= 'Exceedance rate (events/year)')
        abline(h=10**(seq(-5, 0)), lty='dotted', col='brown')
        abline(v=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), lty='dotted', col='brown')
        points(stg_conv_1_small, er_conv_1_small, t='l', col='orange')
        points(stg_conv_2_small, er_conv_2_small, t='l', col='red')
        title(paste0(
            'The two "convergence check" curves are made using half the model scenarios each. \n They should agree fairly well except for rare events.'))

        legend('topright', 
            c('Peak stage exceedance rate (mean over all logic-tree branches)',
              'convergence check 1', 'convergence check 2'),
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

            # Make an outer-margin
            par(oma=c(5,0,0,0))
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
                mtext(side=2, paste0('Magnitude (assumes constant shear modulus)'), line=2.3, cex=0.7)
                title(paste0(sz, ': Rate of events in each magnitude category \n with peak stage > ', 
                    signif(stage_threshold,2), 'm'))
            }

            peak_stage_text = paste0('stage=',signif(stage_threshold, 2), 'm')
            mtext(text=paste0('In the plots containing rate vs magnitude, the black dots give the rate of events for each individual magnitude bin. The red dots give the 95% credible intervals.\n',
                              'The number in parenthesis on the vertical axis (beside the magnitude) gives the percentage of scenarios with that magnitude that exceed ', peak_stage_text, '. If the \n',
                              'latter percentage is reasonably high (e.g. > 20%), then it means that "fairly typical" modelled tsunamis with the specified magnitude can exceed ', peak_stage_text, '\n', 
                              'However, if the percentage is low, it means that only "extreme" modelled tsunamis with that magnitude are exceeding ', peak_stage_text, '. We suggest avoiding the latter\n',
                              "case if possible. In general, we should be more skeptical about the model's representation of unusual or extreme events, as this is more difficult to test.\n"),
                side=1, outer=TRUE, padj=0.85, adj=0.05, cex=0.9)

            # Back to default outer-margins
            par(oma=c(0,0,0,0))
        }
    }

    pdf(paste0('station_summary_', lon, '_', lat, '.pdf'), width=12, height=7)
    .plot_preamble()
    plot_stage_vs_rate()
    plot_deaggregation_summary(0.3)
    plot_deaggregation_summary(1.0)
    plot_deaggregation_summary(2.0)
    dev.off()

}
