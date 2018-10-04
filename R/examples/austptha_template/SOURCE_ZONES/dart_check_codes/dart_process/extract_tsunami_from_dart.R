#
# This code extracts tsunami from manually chosen DART buoys during manually
# chosen time-frames. This allows us to narrow our focus to a few relevant
# gauges, and also to use gauges which are missing from the metadata [e.g. 
# one of the gauges that recorded the 2009 Puysegur earthquake is missing
# from the dart metadata I downloaded, but we have the file and can manually
# deal with that].
#

# Function to read a single .txt file with a year of DART data
parse_dart_year<-function(dart_file){
    fileinfo = read.table(dart_file, comment.char='#', skip=1, header=FALSE, 
        colClasses=c(rep('character', 6), rep('numeric', 2)))

    names(fileinfo) = c('year', 'month', 'day', 'hour', 'min', 'sec', 'T', 
        'height')

    dart_timestring = paste0(fileinfo$year, '-', fileinfo$month, '-', 
        fileinfo$day, ' ', fileinfo$hour, ':', fileinfo$min, ':', fileinfo$sec)
    # Convert time to numeric
    dart_time = strptime(dart_timestring,
        format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

    hts = as.numeric(fileinfo$height)
    tokeep = which(hts != 9999.0)

    output = data.frame(time=dart_time[tokeep], height=hts[tokeep])

    return(output)
}

#'
#' Function to extract data from a single dart buoy, between start_date and
#' end_date
#'
#' @param dart_ID character giving that part of the dart filename that
#' comes before 't', e.g. '55012'
#' @param start_date character giving the start Y-m-d to extract, e.g.
#' '2009-07-14'
#' @param end_date character giving the end Y-m-d to extract, e.g. '2009-07-18'
#' @param remove_time_repetitions I discovered a dataset which seems to
#' repeat a day, with almost the same data. So time does not move forward
#' exclusively! As a work around I introduced this argument. Switched off
#' by default because a large amount of prior work was done without this -
#' those results are obviously OK - and I don't want to by chance mess
#' with them (e.g. due to an un-noticed single point time reversal)
extract_tsunami_from_dart_series<-function(dart_ID, start_date, end_date,
    remove_time_repetitions=FALSE){

    start_end_date = as.POSIXct(strptime(c(start_date, end_date), 
        format='%Y-%m-%d', tz='Etc/GMT'))

    year_start = format(start_end_date[1], '%Y')
    year_end = format(start_end_date[2], '%Y')

    if(year_start != year_end){
        stop(paste0(
            'Code currently needs start_date and end_date in the same year, \n', 
            'although a little new code could fix this'))
    }

    dart_file_name = paste0('dart_historical/', dart_ID, 't', year_end, '.txt')    
    if(!file.exists(dart_file_name)){
        stop('Could not find dart file')
    }

    dart_data = parse_dart_year(dart_file_name)

    xlim = start_end_date

    kk = which(dart_data$time > xlim[1] & dart_data$time < xlim[2])
    if(length(kk) == 0) stop('No dart data found in these time limits')

    if(remove_time_repetitions){
        #browser()
        m1 = as.numeric(diff(dart_data$time[kk]))
        if(any(m1 < 0)){
            l = min(which(m1 < 0))
            if(l > 1){
                kk = kk[1:(l-1)]
            }else{
                stop('Time reversal at the start of the time series! Need to write code to deal with this case!')
            }
        }
    }

    data = data.frame(time=dart_data$time[kk], height=dart_data$height[kk], 
        id=rep(dart_ID, length(kk)))

    return(data)
}

#' Remove the tide from a dart record with loess smoothing 
#'
#' Requires a subjective choice of the 'loess span' smoothing parameter, 
#' see ?loess for guidelines.
#'
#' @param dart_data result of extract_tsunami_from_dart_series
#' @param loess_span span parameter for loess smoothing. This is chosen
#' subjectively, which is why the code produces plots by default.
#' @param plot logical Make a plot?
#' @param zoom_plot_inds indices of the dart series to make a zoomed plot of
#' the de-tided tsunami
#' @param hours_before_after vector of length 2, giving the number of hours
#' before/after the maximum absolute residual to include in the zoom plot
#' @param interpolate_over list with 2 vectors 'start' and 'end' giving the 
#' numeric julian times of one or more patches of data which should be removed prior
#' to interpolation. This lets us 'jump' spikes in the data
#' @param start_limit numeric julian time, before which we do not wish to treat
#' the series as a tsunami. Useful to e.g. dodge seismic waves. The output incudes
#' a vector 'allowed' which is TRUE if the time is after the start_limit.
#' @return a data.frame with the original series, the smoothed data, and the 
#' residual (the latter being the 'tsunami' if it is filtered well)
remove_tide_from_dart_series<-function(dart_data, loess_span = 0.05,
    plot=TRUE, zoom_plot_inds = NULL, hours_before_after = NULL, interpolate_over=NULL,
    start_limit=NULL, ...){

    dart_time = dart_data[,1]
    dart_t = as.numeric(julian(dart_data$time))
    dart_h = dart_data$height
    dart_id = dart_data$id

    # 'Patch' regions by removing them before the interpolation
    if(!is.null(interpolate_over)){

        if(!is.null(zoom_plot_inds)) print('Beware: zoom_plot_inds may interact badly with interpolate_over')

        to_remove = c()
        for(i in 1:length(interpolate_over$start)){
            to_remove = c(to_remove, 
                which(dart_t > interpolate_over$start[i] & dart_t < interpolate_over$end[i]))
        }
        #browser()
        if(length(to_remove) == 0) stop('interpolate_over was provided but no points match')
        dart_t = dart_t[-to_remove]
        dart_h = dart_h[-to_remove]
        dart_id = dart_id[-to_remove]
        dart_time = dart_time[-to_remove]
    }

    # Interpolate evenly spaced at 15s time interval
    np = round((max(dart_t) - min(dart_t))*(3600 * 24 / 15)) + 1
    dart_interp = spline(dart_t, dart_h, n = np)

    loess_fit = loess(y ~ x, data=dart_interp, span=loess_span)
    loess_fitted = predict(loess_fit, data.frame(x = dart_t))

    if(!is.null(start_limit)){
        allowed = (dart_t >= start_limit)
    }else{
        allowed = rep(TRUE, length(dart_t))
    }

    output = data.frame(
        time = dart_time, 
        juliant = dart_t, 
        id = dart_id,
        height=dart_h, loess = loess_fitted, 
        resid = dart_h - loess_fitted, 
        zoom = rep(0, length(dart_t)),
        allowed = allowed)

    if(plot){

        if(is.null(zoom_plot_inds)){
            if(is.null(hours_before_after)) hours_before_after = c(3,6)
            before_hours = as.difftime(hours_before_after[1], format='%H', 
                units='hours')
            after_hours = as.difftime(hours_before_after[2], format='%H', 
                units='hours')
            one_day = as.difftime('24', format='%H', units='hours')

            ll = length(output[,1])
            kk = which.max( abs(output$resid) * 
                (output$time - output$time[1] > one_day) * 
                (output$time[ll] - output$time > one_day) )

            zoom_plot_inds = which(
                (output$time[kk] - output$time < before_hours) &
                (output$time - output$time[kk] < after_hours))
        }

        output$zoom[zoom_plot_inds] = 1

        par(mfrow=c(4,1))
        plot(output$time, output$height, t='l', 
            main=paste0('DART series and loess smooth, ', dart_data$id[1]))
        points(output$time, output$loess, t='l', col='red')
        abline(v=output$time[zoom_plot_inds[1]], col='orange')
        abline(v=output$time[zoom_plot_inds[length(zoom_plot_inds)]], 
            col='orange')

        plot(output$time[zoom_plot_inds], output$height[zoom_plot_inds], 
            t='l', main='DART series_and_loess_smooth_zoom')
        points(output$time[zoom_plot_inds], output$loess[zoom_plot_inds], 
            t='l', col='red')
        grid()
        abline(v=output$time[zoom_plot_inds[1]], col='orange')
        abline(v=output$time[zoom_plot_inds[length(zoom_plot_inds)]], 
            col='orange')

        plot(output$time[zoom_plot_inds], output$resid[zoom_plot_inds], 
            main='DART residual zoom', ...)
        grid()
        abline(v=output$time[zoom_plot_inds[1]], col='orange')
        abline(v=output$time[zoom_plot_inds[length(zoom_plot_inds)]], 
            col='orange')

        kk = min(which(output$allowed))
        abline(v=output$time[kk], col='blue', lty='dashed')

        hourly_seq = as.difftime(seq(-200, 200), format='%H', units='hours')
        abline(v=output$time[zoom_plot_inds[1]] + hourly_seq, col='green')

        # This is a helpful diagnostic -- derivative of the loess -- it will
        # show up if we are smoothing too little
        plot(output$time[-1], diff(output$loess)/diff(output$juliant), t='l')
        title('Time-derivative of loess smoother')
        abline(v=output$time[zoom_plot_inds[1]], col='orange')
        abline(v=output$time[zoom_plot_inds[length(zoom_plot_inds)]], 
            col='orange')
    }

    return(output)
}


#' Convenience function combining the two above functions, with some IO
extract_and_detide_dart<-function(event_name, dart_ID, 
    start_date, end_date, 
    loess_span, topdf = FALSE,  zoom_plot_inds=NULL, 
    hours_before_after=NULL, interpolate_over = NULL,
    start_limit=NULL, remove_time_repetitions = FALSE, ...){

    dir.create(paste0('dart_extract/', event_name), showWarnings=FALSE, 
        recursive=TRUE)

    dart = extract_tsunami_from_dart_series(dart_ID, start_date, end_date,
        remove_time_repetitions=remove_time_repetitions)

    if(topdf){
        pdf(paste0('dart_extract/', event_name, '/', event_name, '_', dart_ID, 
                '.pdf'), 
            width=10, height=12.5)
    }

    #interpolate_over = interpolate_over
    #start_limit = start_limit

    dart_filt = remove_tide_from_dart_series(
        dart_data=dart, 
        loess_span=loess_span,
        zoom_plot_inds = zoom_plot_inds, 
        hours_before_after = hours_before_after,
        interpolate_over = interpolate_over,
        start_limit=start_limit,
        ...)

    if(topdf) dev.off()
   
    write.csv(dart_filt, 
        file=paste0('dart_extract/', event_name, '/', event_name, '_', 
            dart_ID, '.csv'),
        row.names=FALSE)

    return(invisible(dart_filt))
} 

###############################################################################
# Main code
stop('Deliberate stop')


################################################################################
# Kaikoura Mw 7.9 2016/11/13
#
## interpolate_over = NULL
## start_limit = NULL
## # Nothing at 51426, 51425, 55012, 55015, 55023
## x = extract_and_detide_dart('kaikoura_2016_11_13_Mw7.9',
##     '55023',
##     '2016-11-11', '2016-11-17',
##     topdf=FALSE, loess_span=0.015, # DO NOT SEND TO PDF
##     hours_before_after=c(2,6.),
##     interpolate_over = interpolate_over,
##     start_limit=start_limit,
##     ylim=c(-1,1)*0.01, t='l')




################################################################################
# Kuril Mw 8.3 2006/11/15
#

# This one is problematic due to data truncation an hour after the tsunami begins
# This prevents my loess method from working. 
interpolate_over = NULL
start_limit = 13467.50208333333394
x = extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '21414',
    '2006-11-13', '2006-11-18',
    topdf=FALSE, loess_span=0.01, # DO NOT SEND TO PDF HERE, THIS ONE IS UNUSUAL SO REQUIRES FIX BELOW
    hours_before_after=c(2,4),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.08, t='l')
#
# HACK HERE -- just do a linear de-tiding over the (short) tsunami record, in this problematic case.
#
inds = 299:440
gradient = -0.06/(1/24)
x$loess[inds] = x$loess[inds[1]] + gradient*(x$juliant[inds] - x$juliant[inds[1]])
x$resid[inds] = x$height[inds] - x$loess[inds]
x$zoom = 0
x$zoom[inds] = 1
plot(x$juliant, x$height - x$loess, t='l', xlim=c(13467.5, 13467.6))
write.csv(x, 
    file=paste0('dart_extract/', 'kuril_2006_11_15_Mw8.3', '/', 'kuril_2006_11_15_Mw8.3','_', 21414, '.csv'),
    row.names=FALSE)



# Another gauge

interpolate_over = NULL
start_limit = 13468.30208333333394
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '32401',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(4,8),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.04, t='l')

start_limit = 13467.596527777777737
interpolate_over = list(start=13467.656852715701461, end=13467.659138987826736)
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46402',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(1.5,10.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

start_limit = 13467.596527777777737 + 1/24
interpolate_over = list(start=13467.650591457917471, end=13467.653617492413105)
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46403',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(9.5,-0.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

start_limit = 13467.82291666666606 + 0.333/24
interpolate_over = NULL
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46404',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(1.0,5.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

start_limit = 13467.578472222221535 + 0.2/24
interpolate_over = NULL
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46408',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(1.5,4.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

start_limit = 13467.683854166667516 + 1.2/24
interpolate_over = NULL 
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46409',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,6.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

start_limit = 13467.668229166667516 + 2/24
interpolate_over = NULL 
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46410',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.035,
    hours_before_after=c(0,4.0),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')

start_limit = 13467.668229166667516 + 2.7/24
interpolate_over = list(start=13467.802761459206522, end= 13467.803921227436149)
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46411',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.035,
    hours_before_after=c(6,3.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')

start_limit = 13467.79375000000072 + 0.66/24
interpolate_over = NULL 
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46412',
    '2006-11-13', '2006-11-17',
    topdf=TRUE, loess_span=0.035,
    hours_before_after=c(6,3.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')

start_limit = 13467.79375000000072 - 6/24
interpolate_over = NULL 
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '46413',
    '2006-11-13', '2006-11-17',
    topdf=TRUE, loess_span=0.035,
    hours_before_after=c(6,7.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.07, t='l')

# Might be something in this one, but a few spikes...
start_limit = 13467.7472222 - 2/24
interpolate_over = list()
interpolate_over$start = c(13467.6887, 13467.7409241, 13467.7510853)
interpolate_over$end = c(13467.6969, 13467.7438878,  13467.7574360)
extract_and_detide_dart('kuril_2006_11_15_Mw8.3',
    '51407',
    '2006-11-13', '2006-11-18',
    topdf=TRUE, loess_span=0.025,
    hours_before_after=c(2.5,2.5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

##############################################################################
# South America 2016/04 Mw 7.83
# maybe {32411, 32413, 43413} 
interpolate_over=NULL
start_limit=16908.009 + 1.75/24
extract_and_detide_dart('southamerica_2016_04_16_Mw7.8',
    '32411',
    '2016-04-15', '2016-04-20',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,4),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')

interpolate_over=NULL
start_limit=16908.009 + 1.75/24
extract_and_detide_dart('southamerica_2016_04_16_Mw7.8',
    '32413',
    '2016-04-15', '2016-04-20',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,4),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')

###############################################################################
# Solomons 2016 Mw 7.83
#### def{ 52406, 55012, 55023}, weak{51425}
interpolate_over=NULL
start_limit=17143.7569444
extract_and_detide_dart('solomons_2016_12_08_Mw7.8',
    '55012',
    '2016-12-06', '2016-12-10',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,5),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

interpolate_over=NULL
start_limit=17143.78125
extract_and_detide_dart('solomons_2016_12_08_Mw7.8',
    '55023',
    '2016-12-06', '2016-12-10',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,9),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

interpolate_over=NULL
start_limit=17143.7590278
extract_and_detide_dart('solomons_2016_12_08_Mw7.8',
    '52406',
    '2016-12-06', '2016-12-10',
    topdf=TRUE, loess_span=0.03,
    hours_before_after=c(0,3),
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')

## This one only has a few mm -- skip for now    
## interpolate_over=NULL
## start_limit=NULL
## x = extract_and_detide_dart('solomons_2016_12_08_Mw7.8',
##     '51425',
##     '2016-12-06', '2016-12-10',
##     topdf=FALSE, loess_span=0.03,
##     hours_before_after=c(-0,8),
##     interpolate_over = interpolate_over,
##     start_limit=start_limit,
##     ylim=c(-1,1)*0.01, t='l')



################################################################
# South America (Tocopilla earthquake) 2007
interpolate_over=NULL
start_limit=13831.666941 + 1/24 * 1/6
extract_and_detide_dart('southamerica_2007_11_14_Mw7.8', 
    '32401', 
    '2007-11-13', '2007-11-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(0,4), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')

interpolate_over=NULL
start_limit=13831.666941 + 1/24 * (1 + 1/6)
extract_and_detide_dart('southamerica_2007_11_14_Mw7.8', 
    '32412', 
    '2007-11-13', '2007-11-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(0,4.2), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')


################################################################
# KermadecTonga outer rise 2011
interpolate_over=NULL
# Start at 7.25 pm. The wave shouldn't reach the dart before then, and this will remove some seismic waves
start_limit=15161 + (12 + 7 + (25/60) )/24 
extract_and_detide_dart('kermadectonga_2011_07_06_Mw7.6', 
    '54401', 
    '2011-07-04', '2011-07-10', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,2), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

#################################################################
# Solomons 2007
interpolate_over = list()
interpolate_over$start = 13605.07552
interpolate_over$end = 13605.0769
start_limit=13605 + 0.75/24
extract_and_detide_dart('solomons_2007_04_01_Mw8.1', 
    '52402', 
    '2007-03-30', '2007-04-03', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,10), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')

# Note -- here we remove repeated data in the DART series!
interpolate_over = list()
interpolate_over$start = 13605.0949653
interpolate_over$end = 13605.0963542
start_limit = 13605-1/24
extract_and_detide_dart('solomons_2007_04_01_Mw8.1', 
    '52403', 
    '2007-03-30', '2007-04-03', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,9), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    remove_time_repetitions=TRUE,
    ylim=c(-1,1)*0.01, t='l')

## NOTHING HERE!
## interpolate_over=NULL
## start_limit=NULL
## extract_and_detide_dart('solomons_2007_04_01_Mw8.1', 
##     '21413', 
##     '2007-03-30', '2007-04-03', 
##     topdf=FALSE, loess_span=0.05, 
##     hours_before_after=c(1,15), 
##     interpolate_over = interpolate_over,
##     start_limit=start_limit,
##     ylim=c(-1,1)*0.01, t='l')


######################################################################
# South America 2007 Mw 8.0
interpolate_over=list()
interpolate_over$start=13741.1828
interpolate_over$end=13741.1842
start_limit=13741.14
extract_and_detide_dart('southamerica_2007_08_15_Mw8.0', 
    '32411', 
    '2007-08-13', '2007-08-18', 
    topdf=TRUE, loess_span=0.02, 
    hours_before_after=c(1,4.5), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')

start_limit = 13741.01
interpolate_over = NULL
extract_and_detide_dart('southamerica_2007_08_15_Mw8.0', 
    '32401', 
    '2007-08-13', '2007-08-18', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,5.5), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 13741.2 + 1/24
interpolate_over=NULL
extract_and_detide_dart('southamerica_2007_08_15_Mw8.0', 
    '43412', 
    '2007-08-13', '2007-08-18', 
    topdf=TRUE, loess_span=0.02, 
    hours_before_after=c(1,4), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.02, t='l')

################################################################
# Puysegur 2009
interpolate_over=NULL
start_limit = 14440.4088
extract_and_detide_dart('puysegur_2009_07_15_Mw7.8', 
    '55013', 
    '2009-07-13', '2009-07-18', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,2), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

extract_and_detide_dart('puysegur_2009_07_15_Mw7.8', 
    '55015', 
    '2009-07-13', '2009-07-18', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,2), 
    interpolate_over = interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')


################################################################
# Tonga march 2009
start_limit = 14322.7820842
interpolate_over=NULL
extract_and_detide_dart('tonga_2009_03_19_Mw7.7', 
    '51426', 
    '2009-03-18', '2009-03-22', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,4), 
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')


################################################################
# North vanuatu October 2009
start_limit = 14525.0
interpolate_over=NULL
extract_and_detide_dart('vanuatu_north_2009_10_07_Mw7.8', 
    '51425', 
    '2009-10-06', '2009-10-9', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,5),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.03, t='l')
#
# Can't see anything but noise at 51426 or 52401 or 52404 or 21413
# (these sites looked possible in the larger figure)
#
#x = extract_and_detide_dart('vanuatu_north_2009_10_07_Mw7.8', 
#    '21413', 
#    '2009-10-06', '2009-10-9', 
#    topdf=FALSE, loess_span=0.05, 
#    hours_before_after=c(1,8),
#    interpolate_over=interpolate_over,
#    start_limit=start_limit,
#    ylim=c(-1,1)*0.01, t='l')


###############################################################
# Samoa september 2009
start_limit = 14516.780
interpolate_over=NULL
extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
    '51425', 
    '2009-09-27', '2009-10-01', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 14516.7721052
interpolate_over=NULL
extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
    '51426', 
    '2009-09-27', '2009-10-01', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 14516.8131952
interpolate_over = NULL
extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
    '54401', 
    '2009-09-27', '2009-10-01', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 14517.3171447
interpolate_over=NULL
extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
    '32401', 
    '2009-09-27', '2009-10-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 14517.2431447
interpolate_over=NULL
extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
    '32412', 
    '2009-09-27', '2009-10-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

# Sites nearer aust -- do they get a signal? Commented out if not.
#extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
#    '55023', 
#    '2009-09-27', '2009-10-02', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(1,8),
#    ylim=c(-1,1)*0.05, t='l')
#extract_and_detide_dart('samoa_2009_09_29_Mw8.1', 
#    '55015', 
#    '2009-09-27', '2009-10-02', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(1,3),
#    ylim=c(-1,1)*0.1, t='l')


###########################################################################
# Chile 2010_02_27


start_limit = 14667.314704 + 1.5/24
interpolate_over = NULL
xx = extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '32412', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.08, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 21413
start_limit = 14668.11 + 0.65/24
interpolate_over = list()
interpolate_over$start = 14668.1789931
interpolate_over$end = 14668.1809028
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '21413', 
    '2010-02-26', '2010-03-02', 
    topdf=TRUE, loess_span=0.08, 
    hours_before_after=c(2,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 43412
start_limit = 14667.64 + 0.5/24
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '43412', 
    '2010-02-26', '2010-03-02', 
    topdf=TRUE, loess_span=0.08, 
    hours_before_after=c(1,4),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 46403
start_limit = 14667.8020833 + 4.5/24
interpolate_over = list()
interpolate_over$start = 14668.0227431 - 2/(60*24)
interpolate_over$end = 14668.0227431 + 2/(60*24)
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '46403', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(5,10),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

# 46404
start_limit = 14667.314704
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '46404', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(-0.1,10),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 46407
start_limit = 14667.314704
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '46407', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,13),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 46409
start_limit = 14667.95 + 0.75/24
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '46409', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,8),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.1, t='l')

# 46410 is a mess, not using it
# 46412 ""
# 46419
start_limit = 14667.314704
interpolate_over = list()
interpolate_over$start = 14667.9331597
interpolate_over$end = 14667.9354167
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '46419', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,9),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')


start_limit = 14667.8560693
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '51425', 
    '2010-02-25', '2010-03-01', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

start_limit =14667.7419026
interpolate_over = NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '51426', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,24),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

# Good already!
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '52401', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,24),
    ylim=c(-1,1)*0.10, t='l')

# 
start_limit = 14668 + 2/24
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '52402', 
    '2010-02-28', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

# 
start_limit = 14668 + 3.5/24
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '52403', 
    '2010-02-28', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(1,12),
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

# 
start_limit = 14668.15 +2/24
x = extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '52405', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(4,12),
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

# 
start_limit = 14667.6 + 3.5/24
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '54401', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(4,12),
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

start_limit = 14667.97 + 0.33/24
interpolate_over=NULL
extract_and_detide_dart('chile_2010_02_27_Mw8.8', 
    '55012', 
    '2010-02-25', '2010-03-02', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(4,12),
    interpolate_over=interpolate_over,
    start_limit=start_limit,
    ylim=c(-1,1)*0.10, t='l')

#############################################################################
# Andaman 2009

interpolate_over=list()
interpolate_over$start = 14466.8870
interpolate_over$end = 14466.8891
start_limit = 14466.85

extract_and_detide_dart('andaman_2009_08_10_Mw7.5', 
    '23401', 
    '2009-08-08', '2009-08-12', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(4,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')



##############################################################################
# Sumatra 2010 

start_limit = 14705.9745146
interpolate_over=NULL
extract_and_detide_dart('sumatra_2010_04_06_Mw7.8', 
    '23401', 
    '2010-04-03', '2010-04-09', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(4,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')

## ON REFLECTION -- THE SIGNAL IS THE SAME SIZE AS THE NOISE! REMOVE
#start_limit = 14705.98
#interpolate_over=NULL
#extract_and_detide_dart('sumatra_2010_04_06_Mw7.8', 
#    '56001', 
#    '2010-04-03', '2010-04-09', 
#    topdf=TRUE, loess_span=0.02, 
#    hours_before_after=c(4,12),
#    interpolate_over=interpolate_over,
#    start_limit = start_limit,
#    ylim=c(-1,1)*0.01, t='l')

#################################################################################
# Mentawai 2010
start_limit = 14907.6645061
interpolate_over=NULL
extract_and_detide_dart('mentawai_2010_10_25_Mw7.9', 
    '56001', 
    '2010-10-24', '2010-10-27', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(4,15),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')

#############################################################################
# Tonga 2006 -- note: From the T2 report (BOM), it seems that higher res
# data is available.
#extract_and_detide_dart('tonga_2006_05_03_Mw8.0', 
#    '51407', 
#    '2006-05-01', '2006-05-06', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(4,12),
#    ylim=c(-1,1)*0.01, t='l')

# Here we force it to use the high res copy (for details see ./dart_51407_2005_2007)
start_limit = 13271.291493055556202 + 14/24
interpolate_over=NULL
extract_and_detide_dart('tonga_2006_05_03_Mw8.0', 
    '51407B', 
    '2006-05-01', '2006-05-06', 
    topdf=TRUE, loess_span=0.015, 
    hours_before_after=c(6,24),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.02, t='l')

##############################################################################
# Sumatra 2007. 

start_limit = 13768.4958277 + 1/24
interpolate_over=list()
# Remove the seismic waves, and also waves from a second earthquake around
# 12 hours after the first
interpolate_over$start= c( 13768.6651716, 13768.996527)
interpolate_over$end  = c( 13768.6719454, 13769.0100694)
extract_and_detide_dart('sumatra_2007_09_12_Mw8.5', 
    '23401', 
    '2007-09-10', '2007-09-15', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(10,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

##############################################################################
# Japan 2011
start_limit= 15044.2648337
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21401', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,36),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.5, t='l')

start_limit = 15044.273443
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21413', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,36),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.5, t='l')

start_limit = 15044.303443 + 1.3/24
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21414', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,36),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 15044.303443 + 1/24
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21415', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,36),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

# This one is very 'spikey' so we de-spike here
start_limit = 15044.2522684
interpolate_over=list()
interpolate_over$start = c(15044.3122096, 15044.3326772, 15044.3575309,
    15044.3750746, 15044.4576764, 15044.5015358, 15044.5863305, 
    15044.6272659, 15044.7508031, 15044.9174686, 15045.0029943)
interpolate_over$end = c(15044.3158645, 15044.3370632, 15044.3604548,
    15044.3765366, 15044.4606004, 15044.5066527, 15044.5892545,
    15044.6316518, 15044.7551890, 15044.9203926, 15045.0066493)
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21418', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*2.0, t='l')

start_limit = 15044.2707264
interpolate_over=NULL
xx = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '21419', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.7, t='l')

# 
start_limit = 15044.96 + 3/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '32401', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(5,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.74 + 4.33/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '32411', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(5,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.83 + 4.5/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '32412', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(12,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')


# 
start_limit = 15044.81 + 3.25/24
x = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '32413', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.56 + 5.25/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '43412', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.62 + 5.5/24
x = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '43413', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 15044.4 + 0.5/24
interpolate_over = list()
interpolate_over$start=15044.4700906 
interpolate_over$end=15044.4715294
xx = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46402', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,16),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.29 +4/24
x = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46403', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,16),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

#
start_limit = 15044.47 + 3/24
x = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46404', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,24),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit=15044.287 + 3/24
interpolate_over=NULL
xx = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46408', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,24),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit=15044.287 + 5/24
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46409', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,24),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.1, t='l')

# 
start_limit = 15044.49 + 1/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46410', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,24),
    start_limit = start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit=15044.287 + 7.5/24
interpolate_over=NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46411', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,35),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 15044.6 + 1.42/24
interpolate_over = list()
interpolate_over$start = c(15044.5847222-0.002, 15044.6217319, 15044.6625-0.001, 15044.6791667 - 0.0005, 15044.6845806 - 0.0005, 15044.7018715 - 0.0005, 15044.7226792 - 0.001, 15044.6661027, 15044.921378 - 0.001)
interpolate_over$end   = c(15044.5847222+0.000, 15044.6315348, 15044.6625+0.001, 15044.6791667 + 0.0005, 15044.6845806 + 0.0005, 15044.7018715 + 0.0005, 15044.7226792 + 0.001, 15044.6681660, 15044.921378 + 0.001)
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '46412', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(1,12),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.15, t='l')

# 
start_limit = 15044.32 + 5.25/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '51407', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

# 
start_limit = 15044.41 + 3.75/24
x = extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '51425', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    start_limit=start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 15044.35 + 0.5/24
interpolate_over = NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '52402', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,35),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.3, t='l')


# 
start_limit = 15044.31 + 3/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '52403', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,35),
    start_limit=start_limit,
    ylim=c(-1,1)*0.3, t='l')

start_limit = 15044.3 + 1.5/24
interpolate_over = NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '52405', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 15044.28 + 5/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '52406', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    start_limit = start_limit,
    ylim=c(-1,1)*0.15, t='l')

start_limit = 15044.53 + 1/24
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '55012', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

start_limit = 15044.3 + 7/24
interpolate_over = NULL
extract_and_detide_dart('tohoku_2011_03_11_Mw9.1', 
    '55023', 
    '2011-03-09', '2011-03-16', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,35),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')


##############################################################################
# South america 2014

start_limit = 16161.9945
interpolate_over=NULL
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '32401', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(6,35),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 16161.9957 + 0.5/24
interpolate_over=NULL
xx = extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '32402', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(2,10),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')


start_limit = 16162.04 + 0.5/24
interpolate_over=NULL
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '32412', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,10),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit = 16162.48 + 1.5/24
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '51426', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(4,10),
    start_limit = start_limit,
    ylim=c(-1,1)*0.02, t='l')

# Hard to do this one confidently?
start_limit = 16162.48 + 5.5/24
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '21414', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(5,1),
    #interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')


# Hard to do this one confidently?
start_limit = 16162.48 + 8/24
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '21418', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(3,2),
    #interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')

# Hard to do this one confidently?
start_limit = 16162.75 + 1/24
x = extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '52402', 
    '2014-03-28', '2014-04-04', 
    topdf=FALSE, loess_span=0.02, 
    hours_before_after=c(6,0.5),
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')

# Hard to be confident?
start_limit = 16162.04 + 3.5/24
interpolate_over=NULL
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '32411', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.015, 
    hours_before_after=c(7,1),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.01, t='l')

## This one is clear
interpolate_over = list()
# Get rid of seismic waves and spikes
interpolate_over$start = c(16161.9951389, 16162.2246528 -0.001)
interpolate_over$end   = c(16162.0026042, 16162.2246528 +0.001)
start_limit = 16162.0 + 3/24
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '32413', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.015, 
    hours_before_after=c(2,7),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.03, t='l')

# Clear
start_limit = 16162.65625 + 2.24/24
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '52406', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.025, 
    hours_before_after=c(4,10),
    start_limit = start_limit,
    ylim=c(-1,1)*0.02, t='l')

# Clear
start_limit = 16162.65625 + 2.24/24 - 5/24
interpolate_over=list()
interpolate_over$start = 16162.5579861
interpolate_over$end   = 16162.5583333 
extract_and_detide_dart('southamerica_2014_04_01_Mw8.2', 
    '51407', 
    '2014-03-28', '2014-04-04', 
    topdf=TRUE, loess_span=0.025, 
    hours_before_after=c(4,10),
    start_limit = start_limit,
    interpolate_over = interpolate_over,
    ylim=c(-1,1)*0.02, t='l')

##############################################################################
# south america 2015

# 
start_limit = 16695.47 + 8.5/24
x = extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '21413', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,18),
    start_limit = start_limit,
    ylim=c(-1,1)*0.04, t='l')

# 
start_limit = 16695.63 + 2/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '21414', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,18),
    start_limit = start_limit,
    ylim=c(-1,1)*0.04, t='l')

# 
start_limit = 16695.66 + 2/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '21415', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,18),
    start_limit = start_limit,
    ylim=c(-1,1)*0.04, t='l')

# 
start_limit = 16694.97
interpolate_over = NULL
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '32402', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,18),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit = 16695 + 1/24
interpolate_over = NULL
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '32412', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,18),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

# Already good
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '43412', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,9),
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.62 + 1.5/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '46408', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.64 + 1/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '46413', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.48 + 1.5/24
x = extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '51407', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(2,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.41 + 3.5/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '51425', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.61 + 4.5/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52401', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.61 + 4.5/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52402', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.67 + 4/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52403', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.64 + 7.0/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52404', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# 
start_limit = 16695.73 + 5.5/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52405', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

# Good
start_limit = 16695.52 + 4.25/24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '52406', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(6,9),
    start_limit=start_limit,
    ylim=c(-1,1)*0.05, t='l')

# Not really detected
#extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
#    '55015', 
#    '2015-09-14', '2015-09-19', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(6,9),
#    ylim=c(-1,1)*0.05, t='l')

# Not really detected
#extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
#    '55023', 
#    '2015-09-14', '2015-09-19', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(6,9),
#    ylim=c(-1,1)*0.05, t='l')

# Not really detected
#extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
#    '55042', 
#    '2015-09-14', '2015-09-19', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(6,9),
#    ylim=c(-1,1)*0.05, t='l')

start_limit = 16695.1041667 + 2/24
interpolate_over = NULL
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '32411', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(3,18),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

start_limit = 16695.52 + 3./24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '46403', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(2,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')


start_limit = 16695.52 + 3./24
extract_and_detide_dart('southamerica_2015_09_16_Mw8.3', 
    '46409', 
    '2015-09-14', '2015-09-19', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(2,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.05, t='l')

###############################################################
# northNewHebrides_2013_02_06

start_limit =15742.08
interpolate_over = NULL
extract_and_detide_dart('northNewHebrides_2013_02_06_Mw7.9', 
    '55012', 
    '2013-02-04', '2013-02-08', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,9),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.2, t='l')

start_limit =15742.08
interpolate_over = NULL
extract_and_detide_dart('northNewHebrides_2013_02_06_Mw7.9', 
    '52406', 
    '2013-02-04', '2013-02-08', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,9),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.1, t='l')

start_limit =15742.16 + 1/24
interpolate_over = NULL
extract_and_detide_dart('northNewHebrides_2013_02_06_Mw7.9', 
    '52402', 
    '2013-02-04', '2013-02-08', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(3,9),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.03, t='l')

# Good
start_limit =15742.2 + 1/24
extract_and_detide_dart('northNewHebrides_2013_02_06_Mw7.9', 
    '52401', 
    '2013-02-04', '2013-02-08', 
    topdf=TRUE, loess_span=0.05, 
    hours_before_after=c(3,9),
    start_limit = start_limit,
    ylim=c(-1,1)*0.03, t='l')

start_limit =15742.1
interpolate_over = NULL
extract_and_detide_dart('northNewHebrides_2013_02_06_Mw7.9', 
    '51425', 
    '2013-02-04', '2013-02-08', 
    topdf=TRUE, loess_span=0.03, 
    hours_before_after=c(3,9),
    interpolate_over=interpolate_over,
    start_limit = start_limit,
    ylim=c(-1,1)*0.03, t='l')

##################################################################
# Cascadia_2012_10_28_Mw7.8
#
# This does have minor tsunami in NE Pacific.
#
#extract_and_detide_dart('cascadiaNorth_2012_10_08_Mw7.9', 
#    '51425', 
#    '2012-10-06', '2012-10-10', 
#    topdf=TRUE, loess_span=0.05, 
#    hours_before_after=c(3,9),
#    ylim=c(-1,1)*0.03, t='l')
