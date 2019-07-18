#
# Make plots for benchmark problem 06 in the NTHMP suite.
# Conical Island experiment
#
source('../../../plot.R')

# Case A is optional in the current NTHMP -- and we had
# issues with the runup in that one only (recall in teh NTHMP 2011 report,
# some modellers mention scaling the gauge for case A -- we would have to
# do this too). Here just do B, C.
forcing_cases = c(2, 3)
forcing_cases_name = c('B', 'C')

# Allowed relative error in "max(runup anywhere around island)" between model and measured
ERR_TOL_MAX_RUNUP = 0.2
# Allowed mean of relative error in gauge maxima
ERR_TOL_MAX_GAUGE_MEAN = 0.1
# Allowed mean of relative error in island runup 
ERR_TOL_MAX_RUNUP_MEAN = 0.3

for(model_run in 1:length(forcing_cases)){

    forcing_case = forcing_cases[model_run]
    forcing_case_name = forcing_cases_name[model_run]

    # Run the model
    system(paste0('./BP06 ', forcing_case, ' > outfile.log'))

    #
    # Get the multidomain
    #
    multidomain_dir = rev(Sys.glob('OUTPUTS/RUN*'))[1]
    sink('tmp_sink')
    x = lapply(Sys.glob(paste0(multidomain_dir, '/RUN*')), f<-function(y) get_all_recent_results(y, quiet=TRUE))
    sink()

    #
    # Get gauges data
    #
    if(forcing_case_name == 'A'){
        gauges_data_files = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2a.txt'
        gauge_plot_ylim = c(-1, 1)*0.03
        runup_data_file = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/run2a.txt'

    }else if(forcing_case_name == 'B'){
        gauges_data_files = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2b.txt'
        gauge_plot_ylim = c(-1, 1)*0.06
        runup_data_file = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/run2b.txt'

    }else if(forcing_case_name == 'C'){
        gauges_data_files = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2cnew1.txt'
        gauge_plot_ylim = c(-1, 1)*0.1
        runup_data_file = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/run2c.txt'

    }else{
        stop('unrecognized forcing case')
    }

    gauges_data = read.table(gauges_data_files[1], header=FALSE, skip=8)
    # The wavemaker forcing begins at 22s, corresponding to a model start time of 0
    obs_time = gauges_data[,1] - 22.0
    obs_gauges = gauges_data[,2:9]


    #
    # Plot gauges
    #

    # For each observed gauge (there are 8), store the corresponding model domain
    # index and gaugeID. Note the location of gauges 1-4 changes with the forcing_case,
    # leading to a different set of gauges stored in the inner model domain in each case.
    #gauge_mappings = list(c(1,1), c(1,2), c(1,3), c(1,4), 
    #                      c(2,6), c(2,9), c(2,16), c(2,22))
    gauge_IDs = c(1, 2, 3, 4, 6, 9, 16, 22)

    png(paste0('Gauges_plot_', forcing_case_name, '.png'), width=10, height=9, units='in', res=300)
    par(mfrow=c(4,2))
    par(mar = c(4.5,4.5,2,1))
    err_store = rep(NA, length(gauge_IDs))
    md_gauges = merge_multidomain_gauges(x)
    for(i in 1:length(gauge_IDs)){
        plot(obs_time, obs_gauges[,i], t='p', pch=19, cex=0.2, 
             ylim=gauge_plot_ylim, xlab='Time (s)', ylab='Stage (m)',
             cex.axis=1.4, cex.lab=1.6, xlim=c(0, 40))

        # Get the gauge coordinate
        gi = which(round(md_gauges$gaugeID) %in% gauge_IDs[i])

        # Check the domain for which it has priority. This can change because
        # the gauge locations shift
        points(md_gauges$time, md_gauges$time_var$stage[gi,],t='l', col='red')
        title(paste0('Gauge ', gauge_IDs[i]), cex.main=1.8)
        grid()

        err_store[i] = abs(max(md_gauges$time_var$stage[gi,]) - max(obs_gauges[,i]))/diff(range(obs_gauges[,i]))
    }
    dev.off()


    if(mean(err_store) < ERR_TOL_MAX_GAUGE_MEAN){
        print(c('PASS', mean(err_store)))
    }else{
        print(c('FAIL', mean(err_store)))
    }
    #print(err_store)


    #
    # Plot multidomain elevation
    #
    png(paste0('Flume_plot_', forcing_case_name, '.png'), width=7.2, height=8, units='in', res=300)
    multidomain_image(dirname(x[[1]]$output_folder), variable='elev', time_index=1, 
        xlim=range(x[[1]]$xs), ylim=range(x[[1]]$ys), zlim=c(-0.35, 0.65),
        cols=colorRampPalette(c('lightblue', 'green', 'yellow', 'orange', 'red', 'purple', 'black'))(255))
    points(md_gauges$lon, md_gauges$lat, pch=20)
    text(md_gauges$lon, md_gauges$lat, round(md_gauges$gaugeID), col='purple', 
        pos=c(rep(4,4), 2, rep(4,3)), cex=1.4)
    dev.off()


    #
    # Plot runup around Island
    #
    runup_data = read.table(runup_data_file, skip=10, header=FALSE)
    island_centre = c(12.96, 13.80)

    png(paste0('Runup_plot_', forcing_case_name, '.png'), width=6, height=5, units='in', res=300)

    plot(runup_data[,2], runup_data[,3]/100, t='p', xlab='Degrees around island', ylab='Runup (m)', 
         ylim=c(0, max(runup_data[,3]/100)*1.5), cex.lab=1.5)

    # Get max-stage from around island
    # Find wet cells near a wet/dry boundary
    D_ID = 2
    wet_or_dry = (x[[D_ID]]$elev0 < x[[D_ID]]$maxQ - 1.0e-06) # True if wet
    n = dim(wet_or_dry)
    wet_or_dry_smooth = wet_or_dry
    # Single jacobi iteration.
    wet_or_dry_smooth[2:(n[1]-1), 2:(n[2]-1)] = 0.25 * (
        wet_or_dry_smooth[1:(n[1]-2), 2:(n[2]-1)] + 
        wet_or_dry_smooth[3:(n[1]-0), 2:(n[2]-1)] + 
        wet_or_dry_smooth[2:(n[1]-1), 1:(n[2]-2)] + 
        wet_or_dry_smooth[2:(n[1]-1), 3:(n[2]-0)] )
    near_wd = (wet_or_dry > 0.5)&(wet_or_dry_smooth < 0.8) # Wet cell with at least 1 dry neighbour
    near_wd_keep = which(near_wd, arr.ind=TRUE)
    # x, y, maximum-stage of wet/dry boundary 
    near_wd_x = x[[D_ID]]$xs[near_wd_keep[,1]]
    near_wd_y = x[[D_ID]]$ys[near_wd_keep[,2]]
    near_wd_stage = apply(near_wd_keep, 1, f<-function(y) x[[D_ID]]$maxQ[y[1], y[2]])
    # Use a degree coordinate, like the data
    deg_loc = atan2(near_wd_x - island_centre[1], near_wd_y - island_centre[2])/pi*180
    deg_loc[deg_loc<0] = deg_loc[deg_loc < 0] + 360
    deg_loc_order = order(deg_loc)
    points(deg_loc[deg_loc_order], near_wd_stage[deg_loc_order], t='l', col='red', lwd=2)
    title(paste0('Runup around island, Case', forcing_case_name), cex.main=1.8)
    grid()
    dev.off()

    model_at_obs_site = approx(deg_loc[deg_loc_order], near_wd_stage[deg_loc_order], xout=runup_data[,2], rule=2)
    err_stat = mean(abs(model_at_obs_site$y*100 - runup_data[,3])/runup_data[,3])
    #
    # The NTHMP reports seem to use max(model) - max(obs) for this problem, which means 
    # models that get the max runup quite wrong, but have a comparably large
    # runup elsewhere, still do ok.
    if(err_stat < ERR_TOL_MAX_RUNUP_MEAN){
        print(c('PASS', err_stat))
    }else{
        print(c('FAIL', err_stat))
    }

    err_statB = (max(model_at_obs_site$y*100) - max(runup_data[,3]))/max(runup_data[,3])
    if(err_statB < ERR_TOL_MAX_RUNUP){
        print(c('PASS', err_statB))
    }else{
        print(c('FAIL', err_statB))
    }


}
