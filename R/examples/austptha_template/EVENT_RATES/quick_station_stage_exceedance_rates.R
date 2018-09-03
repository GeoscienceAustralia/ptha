#
# Single station stage exceedance-rate computation
#

library(rptha)
options(scipen=5) # Suppress scientific notation (e.g. 0.0001 rather than like 1e-04)

.preamble_title = paste0('2018 Australian Probabilistic Tsunami Hazard Assessment single station summary')
.preamble_text = paste0( 
                       "This file gives a summary of the Geoscience Australia's 2018 PTHA model results at a single station. See the README on:\n",
                       '    https://github.com/GeoscienceAustralia/ptha/tree/master/ptha_access\n',
                       'for further information on accessing the results, source code, and associated reports describing the methodology and its limitations.\n',
                       '\n',
                       'This product is provided under a Creative Commons 4.0 International Licence (http://creativecommons.org/licenses/by/4.0/legalcode)\n',
                       'Geoscience Australia has tried to make the information in this product as accurate as possible. However, it does not guarantee that the\n',
                       'information is totally accurate or complete. Therefore, you should not solely rely on this information when making a commercial decision.\n',
                       '\n',
                       'If using points far from Australia, note we ignore earthquake source-zones that are not considered relevant for Australian tsunami hazard.\n'
)

.preamble_title1 = 'Information on the plots'
.preamble_text1 = paste0( 
                       'The plots are: \n',
                       '\n',
                       '    1) A hazard curve, containing the peak-tsunami-stage vs the exceedance-rate at the specified hazard point location,\n',
                       '       with logic-tree mean and quantiles illustrating the uncertainty. Note the peak-tsunami-stage is the maximum water level\n',
                       '       attained by the tsunami (above mean-sea-level=0, ignorning tides) at the site.\n',
                       '\n',
                       '    2) A "convergence check" of the above hazard curve. The PTHA hazard results rely on simulating a large number of random\n',
                       '       earthquake-tsunami scenarios. At sufficiently rare exceedance-rates (less frequent than some site specific value, e.g. 1/10000 years)\n', 
                       '       there will be few random model scenarios, so exceedance-rates become sensitive to the "random details" of those scenarios (i.e. less\n', 
                       '       reliable). To help users judge when this happens we compare 2 hazard curves, each derived from half the scenarios.\n',
                       '\n',
                       '    3-8) Information on which source-zones dominate the hazard (i.e. hazard-deaggregation plots) for mean exceedance rates of 1/100,\n', 
                       '       1/500, 1/2500. In each case, for the top 3 source-zones we show rates separated by magnitude (assuming "constant shear modulus"),\n',
                       '       to highlight the model scenarios most likely to cause tsunami above the threshold peak-stage.\n',
                       '       The "spatial hazard deaggregation" plot gives an idea of where earthquakes exceeding the stage threshold might occur. For every\n',
                       '       unit-source, it shows:\n',
                       '           SUM( (event-slip-on-the-unit-source / sum-of-the-event-slip-on-all-the-unit-sources) X individual-mean-event-rate )\n',
                       '       , where the SUM includes events that exceed the peak-stage threshold. Results are normalised to [0-1].\n',
                       '       Please note that the appearance of the "spatial hazard deaggregation" plot is significantly affected by the choice of colour scheme.\n',
                       '       Interpretations should always be cross-checked with the "Top 10 source-zones" bar plot. A common mistake is to focus on\n',
                       '       "small-areas with high values", without considering that "large-areas with lower values" may contribute more to the overall hazard.\n',
                       '       Furthermore, details of the hazard deaggregation may be affected by convergence issues analogous to those mentioned in point 2 above.\n',
                       '\n',
                       '    9-10) The peak-stage at the station of interest resulting from 1m of slip on each unit source. This may help to understand and\n',
                       '       cross-check the hazard deaggregation.\n'
)

.preamble_title2 = 'Information on the plots'
.preamble_text2 = paste0('\n',
                       'An important caveat regarding the interpretation of earthquake magnitude:\n', 
                       '  Earthquake magnitudes in these plots are derived from the event slip and area, assuming a "constant shear modulus" of 30 GPa (for thrust\n',
                       '  fault sources) or 60 GPa (for normal fault sources). However, in reality the shear modulus may vary with depth (e.g. Bilek and Lay, 1999).\n',
                       '  In that case the real earthquake magnitude may differ from the "constant shear modulus magnitude" reported in the Figures.\n', 
                       '  For this reason, our downloadable event metadata also reports on the "variable_mu_Mw" which is the magnitude assuming a depth-dependent\n',
                       '  shear modulus for thrust source-zones only (based on Bilek and Lay, 1999). See the report for further information.\n',
                       '  You should consider such details if comparing events of a particular magnitude with historical data, in a situation where the shear modulus\n', 
                       '  might differ significantly from 30 GPa (for thrust fault source-zones only). On normal fault source-zones the current study does not account\n',
                       '  for variation of the shear modulus with depth.\n',
                       '\n',
                       'An important caveat regarding the rate-by-magnitude plots:\n',
                       '  For the plots which show "rate-by-magnitude" for the top 3 source-zones, the rates (mean, 16%, 84%) are derived by differentiating the \n',
                       '  Mw-exceedance-rate curve at the logic-tree (mean, 16%, 84%), and distributing among scenarios. This ensures a tight relation between \n',
                       '  Mw-exceedance-rate curves and the scenario rates. However, this procedure does NOT ensure any particular ordering on the (mean, 16%, 84%) rates\n',
                       '  at the scenario level. Although the Mw-exceedance-rate curves are ordered, their derivatives (used to derive scenario rates) might not be. \n',
                       '  Irrespective the "rate-by-magnitude" plots give a "heuristic indication" of the magnitudes most likely to lead to the stage-exceedance.'
                       '\n',
                       'A copy of the script used to make this plot can be found at:\n',
                       '    https://github.com/GeoscienceAustralia/ptha/tree/master/R/examples/austptha_template/EVENT_RATES\n',
                       'in the file named "quick_station_stage_exceedance_rates.R"\n',
                       'Other codes used to develop the PTHA can be found in the same git repository.\n',
                       '\n',
                       'References\n',
                       '    Bilek, S.L. and Lay, T. (1999) Rigidity variations with depth along interplate megathrust faults in subduction zones. Nature, 400, 443-446\n'
)

.plot_preamble<-function(){

    plot(c(0,1), c(0,1), col='white', frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    title(.preamble_title, cex=1.5)
    text(0, 1, .preamble_text, adj=c(0,1), cex=0.95)

    plot(c(0,1), c(0,1), col='white', frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    title(.preamble_title1, cex=1.5)
    text(0, 1, .preamble_text1, adj=c(0,1), cex=0.95)

    plot(c(0,1), c(0,1), col='white', frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    title(.preamble_title2, cex=1.5)
    text(0, 1, .preamble_text2, adj=c(0,1), cex=0.95)
}

#
# In one of the figures we use a coarse background raster to show land/water
# This function reads it (and creates it if it does not exist)
#
get_background_raster<-function(){

    bg_raster_file = 'background_raster/land_water_raster.tif'

    if(!file.exists(bg_raster_file)){
        # Make the background raster
        input_raster = '../DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
        modelling_rast = raster(input_raster)
        modelling_rast_lw = (modelling_rast > 0)
        basic_lw = raster(extent(modelling_rast_lw), nrow=nrow(modelling_rast_lw)/25, 
            ncol=ncol(modelling_rast_lw)/25)
        basic_lw = resample(modelling_rast_lw, basic_lw)
        dir.create('background_raster', showWarnings=FALSE)
        writeRaster(basic_lw, bg_raster_file, options=c('COMPRESS=DEFLATE'))
    }

    bg_raster = raster(bg_raster_file)
   
    return(bg_raster) 
}
background_raster = get_background_raster()


#
# For putting the spatial hazard deaggregation results on the map, these
# functions are helpful
#
plot_hazard_curves_utilities = new.env()
source('plot_hazard_curves_utilities.R', local=plot_hazard_curves_utilities)

#' Make a single-station summary plot of the tsunami hazard results
#'
#' @param lon longitude near the gauge to check
#' @param lat latitude near the gauge to check
#' @param output_dir directory for the pdf plot
#'
quick_source_deagg<-function(lon, lat, output_dir='.'){

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

    #
    # Get stage and rates, for each source zone, from the event files.
    #
    stage_rate = vector(mode='list', length=length(tsunami_files))
    names(stage_rate) = basename(dirname(dirname(tsunami_files)))
    for(i in 1:length(tsunami_files)){

        print(basename(tsunami_files[i]))

        fid = nc_open(tsunami_files[i], readunlim=FALSE)
        fid_rates = nc_open(earthquake_only_files[i], readunlim=FALSE)

        event_rate = ncvar_get(fid_rates, 'variable_mu_rate_annual')
        event_rate_upper = ncvar_get(fid_rates, 'variable_mu_rate_annual_upper_ci')
        event_rate_lower = ncvar_get(fid_rates, 'variable_mu_rate_annual_lower_ci')
        event_rate_median = ncvar_get(fid_rates, 'variable_mu_rate_annual_median')
        event_rate_16pc = ncvar_get(fid_rates, 'variable_mu_rate_annual_16pc')
        event_rate_84pc = ncvar_get(fid_rates, 'variable_mu_rate_annual_84pc')

        event_Mw = round(ncvar_get(fid_rates, 'Mw'), 3) # Deal with floating point imperfections in netcdf
        #event_Mw_vary_mu = round(ncvar_get(fid_rates, 'variable_mu_Mw'))
        peak_stage = ncvar_get(fid, 'max_stage', start=c(1, ni), count=c(-1,1))

        site = rep(basename(dirname(dirname(tsunami_files[i]))), length=length(event_rate))
        #row_index = 1:length(event_rate)

        stage_rate[[i]] = data.frame(
            event_rate = event_rate,
            event_rate_upper = event_rate_upper,
            event_rate_lower = event_rate_lower,
            event_rate_median = event_rate_median,
            event_rate_16pc = event_rate_16pc,
            event_rate_84pc = event_rate_84pc,
            peak_stage = peak_stage,
            site = site,
            #row_index=row_index,
            event_Mw = event_Mw)
            #event_Mw_vary_mu = event_Mw_vary_mu)

        nc_close(fid)
        nc_close(fid_rates)
        rm(event_rate, event_rate_upper, event_rate_lower, event_Mw, 
            event_rate_median, event_rate_16pc, event_rate_84pc,
            #event_Mw_vary_mu, row_index,
            peak_stage, site)
        gc()
    }

    # Back-calculate the stage-vs-rate curves
    stage_rate_all = do.call(rbind, stage_rate)
    rm(stage_rate); gc()

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
    ers_median = ncvar_get(fid, 'variable_mu_stochastic_slip_rate_median', start=c(1,ni), count=c(-1,1))
    ers_16pc = ncvar_get(fid, 'variable_mu_stochastic_slip_rate_16pc', start=c(1,ni), count=c(-1,1))
    ers_84pc = ncvar_get(fid, 'variable_mu_stochastic_slip_rate_84pc', start=c(1,ni), count=c(-1,1))
    nc_close(fid)

    # Get stages at various exceedance rates FOR THE MEDIAN CURVE
    #ex_rates = c(1/100, 1/500, 1/2500)
    #stages_at_exceedance_rates_median_curve = approx(ers_median, stages, xout=ex_rates, ties='min')$y

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
    # Plot the stage-vs-rate curves
    #
    plot_stage_vs_rate<-function(exceedance_rates_median_curve=NULL){

        xmax = max(stg*(er>0), na.rm=TRUE)
        ylim = c(1e-05, max(max(er), 1))
        xlim = c(min(0.02, max(xmax-0.005, 1.0e-05)), xmax)


        # The stage-vs-rate info from the event files
        plot(stg_small, er_small, t='l', log='xy', xlim=xlim, ylim=ylim, 
            xlab='Stage (m)', ylab= 'Exceedance rate (events/year)', col='brown')
        points(stg_small, er_up_small, t='l', col='red')
        points(stg_small, er_lo_small, t='l', col='red')

        points(stages, ers, pch=4, cex=1.0, col='brown')
        points(stages, ers_up, pch=19, cex=0.5, col='red')
        points(stages, ers_lo, pch=19, cex=0.5, col='red')
        points(stages, ers_median, pch=19, t='o', col='black')
        points(stages, ers_16pc, pch=19, t='o', col='orange', cex=0.5)
        points(stages, ers_84pc, pch=19, t='o', col='orange', cex=0.5)

        # Note that the lines/points for median/16/84 pc will automatically overlap.
        # However, consistency in the files is required for 95% & mean overlap, so this
        # check is still useful.
        main_title_extra = '\n (Lines and points should overlap for all curves)'
        #main_title_extra = ''

        main_title = paste0('Stage-vs-exceedance-rate @ (lon=', 
            round(hp$lon[ni],3), ', lat=', round(hp$lat[ni], 2), ', elev=', 
            round(hp$elev[ni],2), ', ID=', round(hp$gaugeID[ni], 2),') ', 
            main_title_extra)

        title(main_title)

        abline(h=10**(seq(-5, 0)), lty='dotted', col='brown')
        #abline(h=2*10**(seq(-5, 0)), lty='dotted', col='brown')
        #abline(h=1/100, 1/500, 1/2500, col='yellow')
        abline(v=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), 
            lty='dotted', col='brown')

        legend('topright', 
            c('Peak-stage exceedance-rate (Median from logic-tree)',
               'Mean', '68% credible interval', '95% credible interval'),
            col=c('black', 'brown', 'orange', 'red'),
            pch=c(19, 4, 19, 19), bg='white', 
            pt.cex=c(1, 1, 0.5, 0.5))
  
        # Add convergence check information 
        plot(stg_small, er_small, t='l', log='xy', 
            xlim=c(min(0.02, max(xmax-0.005, 1.0e-05)), xmax),
            ylim=ylim,
            xlab='Stage (m)', ylab= ' Mean exceedance rate (events/year)')
        abline(h=10**(seq(-5, 0)), lty='dotted', col='brown')
        abline(v=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), 
            lty='dotted', col='brown')
        points(stg_conv_1_small, er_conv_1_small, t='l', col='orange')
        points(stg_conv_2_small, er_conv_2_small, t='l', col='red')
        title(paste0(
            'The two "convergence check" curves are made using half the model scenarios each.\n',
            'They should agree fairly well except for rare events.'))

        legend('topright', 
            c('Peak-stage exceedance-rate (mean over all logic-tree branches)',
              'convergence check 1', 'convergence check 2'),
            col=c('black', 'orange', 'red'),
            lty=c(1, 1, 1), bg='white')
    }
    

    #
    # Summarise the distribution of earthquake magnitudes by source-zone
    #
    peak_stage_magnitude_summary<-function(stage_threshold, source_zone){

        is_site = as.character(stage_rate_all$site) == source_zone

        mw_max = max(stage_rate_all$event_Mw[is_site & stage_rate_all$event_rate > 0])
        mw_min = min(stage_rate_all$event_Mw[is_site & stage_rate_all$event_rate > 0])

        if(is.finite(mw_max) & is.finite(mw_min)){
            k = which( (stage_rate_all$peak_stage > stage_threshold) & 
                        is_site & 
                        (stage_rate_all$event_rate>0) &
                        (stage_rate_all$event_Mw >= mw_min) & 
                        (stage_rate_all$event_Mw <= mw_max) )
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
            output_median = aggregate(stage_rate_all$event_rate_median[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x)))
            output_16pc = aggregate(stage_rate_all$event_rate_16pc[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x)))
            output_84pc = aggregate(stage_rate_all$event_rate_84pc[k], 
                by=list(stage_rate_all$event_Mw[k]), 
                f<-function(x) c(sum(x)))


            output = data.frame(Mw=output[,1], 
                rate_exceeding=output$x[,1], 
                n=output$x[,2], 
                rate_exceeding_lower=output_lower$x, 
                rate_exceeding_upper=output_upper$x, 
                rate_exceeding_median=output_median$x,
                rate_exceeding_16pc=output_16pc$x,
                rate_exceeding_84pc=output_84pc$x)
        
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
    # Plot the distribution of earthquake magnitudes by source-zone
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

            rate_by_source_mean = aggregate(stage_rate_all$event_rate[k], 
                by=list(source_zone=as.character(stage_rate_all$site[k])), sum)

            # Sort from highest to lowest
            m1 = order(rate_by_source_mean$x, decreasing=TRUE)
            if(length(m1) > 10){
                # Plot at most 10 source-zones
                rate_by_source_mean = rate_by_source_mean[m1[1:10],]
                m1 = order(rate_by_source_mean$x, decreasing=TRUE)
            }

            # Color the top 3 source-zones differently
            colz = c('grey', 'red')[ 1 + (( (1:length(m1)) %in% m1[1:3]) & (rate_by_source_mean$x > 0))]

            oldmar = par('mar')
            par('mar' = oldmar + c(0, 8, 0, 0)) # new margins to allow source-zone names to fit on plot
            barplot(rate_by_source_mean$x[m1], 
                names.arg=as.character(rate_by_source_mean[m1,1]),
                col=colz[m1], density=100, horiz=TRUE, las=1, 
                xlab='Logic-tree Mean rate (events/year)', 
                main=paste0('Top ', length(m1), 
                    ' sources: Mean rate of events with \n peak_stage > ', signif(stage_threshold,3), 'm', 
                    ' @ (', round(hp$lon[ni],3), ', ', round(hp$lat[ni],3), ')')
                )
            par('mar' = oldmar) # Back to old margins

            # Also plot rate-vs-Mw for the 3 largest contributors 
            for(i in 1:min(length(m1), 3)){
                # Name of source-zone
                sz = rate_by_source_mean[m1[i], 1]
                rate_by_Mw = peak_stage_magnitude_summary(stage_threshold, sz)
                dotchart(rate_by_Mw$rate_exceeding,
                    labels=paste0(rate_by_Mw$Mw, ' (', round(rate_by_Mw$fraction_events*100, 1), ')'),
                    xlab='Rate > stage_threshold (events/year) ', ylab='', 
                    pch=4, col='black', xlim=c(0, max(rate_by_Mw$rate_exceeding_84pc)))
                # Overwrite the black (needed to make labels black)
                points(rate_by_Mw$rate_exceeding, 1:nrow(rate_by_Mw), col='brown', pch=4) 
                # Other summaries
                #points(rate_by_Mw$rate_exceeding_upper, 1:nrow(rate_by_Mw), col='red')
                #points(rate_by_Mw$rate_exceeding_lower, 1:nrow(rate_by_Mw), col='red')
                #points(rate_by_Mw$rate_exceeding_median, 1:nrow(rate_by_Mw), col='black', pch=19)
                points(rate_by_Mw$rate_exceeding_16pc, 1:nrow(rate_by_Mw), col='orange')
                points(rate_by_Mw$rate_exceeding_84pc, 1:nrow(rate_by_Mw), col='orange')
                mtext(side=2, paste0('Magnitude (assumes constant shear modulus)'), line=2.3, cex=0.7)
                title(paste0(sz, ': Rates with peak-stage > ', signif(stage_threshold,3), 'm \n Split by magnitude category'))
            }
            ## Because the 'mean', 'median', etc refer to the "Exceedance Rate" curve, but the scenario weights
            ## are related to the derivative of this, it is possible that the scenario weights are not ordered
            ## as one would naturally expect (e.g. sometimes it is not true that scenarios 16% < mean < 84%). 
            ## I expect this will be confusing, so it is better not to show the median.
            #legend('bottomright', c('Median', 'Mean', '16/84%'), col=c('black', 'brown', 'orange'), pch=c(19, 4, 1))
            legend('bottomright', c('Mean', '16/84%'), col=c('brown', 'orange'), pch=c(4, 1))

            peak_stage_text = paste0('stage=',signif(stage_threshold, 3), 'm')
            mtext(text=paste0(" The rate vs magnitude plots give an indication of which magnitudes are most likely to generate tsunamis exceeding ", peak_stage_text, ". They are derived by partitioning each\n", 
                              "source's magnitude-exceedance rate curves (mean, 16%, 84%) into individual scenario rates, and then summing by magnitude for events that exceed ", peak_stage_text, ".\n",
                              'The number in parenthesis on the vertical axis (beside the magnitude) gives the percentage of scenarios with that magnitude that exceed ', peak_stage_text, '. High values\n',
                              'suggest that typical modelled tsunamis with that magnitude can exceed ', peak_stage_text, ', while low values indicate that unusual events dominate. \n'),
                side=1, outer=TRUE, padj=0.85, adj=0.05, cex=0.9)

            # Back to default outer-margins
            par(oma=c(0,0,0,0))
        }


    }

    pdf(paste0(output_dir, '/station_summary_', lon, '_', lat, '.pdf'), width=12, height=7)
    .plot_preamble()

    plot_stage_vs_rate()
    gc()

    # Plot at a few exceedance rates, median curve
    ex_rates = c(1/100, 1/500, 1/2500)
    stages_ex_rates_mean = approx(ers, stages, xout=ex_rates, ties='min')$y
    median_rates_for_stages   = approx(stages, ers_median, xout=stages_ex_rates_mean)$y
    for(sv in 1:length(ex_rates)){

        #site_deagg = plot_hazard_curves_utilities$get_station_deaggregated_hazard(lon, lat, 
        #    slip_type='stochastic', exceedance_rate=ex_rates[sv], shear_modulus_type='variable_mu_')
        site_deagg = plot_hazard_curves_utilities$get_station_deaggregated_hazard(lon, lat, 
            slip_type='stochastic', stage=stages_ex_rates_mean[sv], shear_modulus_type='variable_mu_')

        stage_level = signif(site_deagg[[1]]$stage_exceed, 3)

        # Spatial hazard plot
        par(mfrow=c(1,1))
        plot_hazard_curves_utilities$plot_station_deaggregated_hazard(site_deagg, scale=0.01,
            background_raster=background_raster, 
            main=paste0('Spatial hazard deaggregation, peak-stage exceeding ', stage_level,  'm',
                ' \n Logic-tree-mean-exceedance-rate = 1/', (1/ex_rates[sv]), 
                '; median-exceedance-rate = 1/', round(1/median_rates_for_stages[sv])))
        rm(site_deagg); gc()

        # Bar charts
        plot_deaggregation_summary(stage_level)
    }

    par(mfrow=c(1,1))

    # Add the peak stage from all unit sources
    site_index = plot_hazard_curves_utilities$plot_unit_source_wave_heights_at_station(
        lon, lat,
        background_raster=background_raster, rake_range = c(89, 91), 
        main='Peak stage from each thrust (rake=90) unit-source with 1m slip')
    gc()
    # Add the peak stage from all unit sources. Use the site index from the
    # previous plot to speed it up
    plot_hazard_curves_utilities$plot_unit_source_wave_heights_at_station(
        lon, lat,
        site_index=site_index,
        background_raster=background_raster, rake_range = c(-91, -89), 
        main='Peak stage from each normal (rake=-90) unit-source with 1m slip')
    gc()

    dev.off()

}



