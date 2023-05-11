# Read the dart buoy timeseries

# This can be run after 'read_dart_metadata.R', as it relies on a shapefile
# created by the latter.

library(rptha) # provides readOGR, writeOGR
dart_metadata = readOGR('dart_locations', layer='dart_locations')

dart_files = list()
for(i in 1:nrow(dart_metadata)){
    dart_name = as.character(dart_metadata$WMO.ID[i])
    dff = Sys.glob(paste0('dart_historical/', dart_name, 't*.txt'))
    if(length(dff) == 0){
        dart_files[[dart_name]] = NA
    }else{
        dart_files[[dart_name]] = dff
    }
}

stopifnot(all(names(dart_files) == as.character(dart_metadata$WMO.ID)))

cat('Manually appending higher res DART51407 in 2006..')
dart_files = c(dart_files, 'dart_historical/51407Bt2006.txt')

# Function to read a single .txt file with a year of DART data
parse_dart_year<-function(dart_file){
    fileinfo = read.table(dart_file, comment.char='#', skip=1, header=FALSE, 
        colClasses=c(rep('character', 6), rep('numeric', 2)))

    names(fileinfo) = c('year', 'month', 'day', 'hour', 'min', 'sec', 'T', 'height')

    dart_timestring = paste0(fileinfo$year, '-', fileinfo$month, '-', fileinfo$day, 
            ' ', fileinfo$hour, ':', fileinfo$min, ':', fileinfo$sec)
    # Convert time to numeric
    dart_time = strptime(dart_timestring,
        format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

    hts = as.numeric(fileinfo$height)
    tokeep = which(hts != 9999.0)

    output = data.frame(time=dart_time[tokeep], height=hts[tokeep])

    return(output)
}

# Read and concatenate all the input data
dart_series = list()
for(i in 1:length(dart_files)){
    dff = dart_files[[i]]
    if(!is.na(dff[1])){
        print(dff[1])
        series_accum = parse_dart_year(dff[1])
        if(length(dff) >= 2){
            for(j in 2:length(dff)){
                print(paste0('   ', dff[j]))
                series_accum = rbind(series_accum, parse_dart_year(dff[j]))
            }
        }
        # Add the filenames as an attribute, to make it harder to link the data
        # with the wrong station later
        attr(series_accum, 'source_files') = dff
        dart_series[[names(dart_files)[i]]] = series_accum
    }else{
        dart_series[[names(dart_files)[i]]] = list()
    }
}

saveRDS(dart_series, 'dart_historical/dart_series.RDS')

# Get a world map polygon as a SpatialPolygonsDataFrame. 
# Originally I got this from the 'maptools' package, but that is retiring as of 2023, 
# so we make a work-around.
get_wrld_simpl<-function(use_maptools = FALSE){

    if(use_maptools){
        require(maptools, quietly=TRUE)
        data(wrld_simpl)
    }else{
        require(spData)
        require(sp)
        require(sf)
        data(world)
        wrld_simpl = as(world, 'Spatial')
    }
    return(wrld_simpl)
}

# Convenience plotting routine
plot_dart_series<-function(dart_series, start_date, end_date){
    dir.create('FIG', showWarnings=FALSE)

    pdf_out = paste0('FIG/event_', start_date, '.pdf')

    wrld_simpl = get_wrld_simpl()

    pdf(pdf_out, width=15, height=12)
    for(i in 1:length(dart_series)){

        if(length(dart_series[[i]]) == 0) next

        xlim = as.POSIXct(strptime(c(start_date, end_date), format='%Y-%m-%d', tz='Etc/GMT'))

        kk = which(dart_series[[i]]$time > xlim[1] & dart_series[[i]]$time < xlim[2])
        if(length(kk) == 0) next

        par(mfrow=c(2,1))
        par(mar=c(4,3,2,1))
        plot(dart_series[[i]]$time[kk], dart_series[[i]]$height[kk], t='l')
        points(dart_series[[i]]$time[kk], dart_series[[i]]$height[kk], col='red', pch='.', cex=2)
        grid()
        title(names(dart_series)[i])
        par(mar=c(0,0,0,0))
        plot(wrld_simpl, asp=1, axes=TRUE)
        plot(dart_metadata, add=TRUE, col='blue')
        points(dart_metadata[i,], col='red', pch=19)
    }
    dev.off()
}

# A small earthquake at the south of Puysegur (can't see anything)
#plot_dart_series(dart_series, '2007-09-28', '2007-10-02')

#stop()

# Kaikoura 2016
plot_dart_series(dart_series, '2016-11-12', '2016-11-14')

# Kuril Islands 2006
plot_dart_series(dart_series, '2006-11-13', '2006-11-17')

# Kermadec 2006 (can't see anything)
# plot_dart_series(dart_series, '2006-05-01', '2006-05-05')

# South America Mw 7.8
plot_dart_series(dart_series, '2016-04-14', '2016-04-18')

# Solomons 2016
plot_dart_series(dart_series, '2016-12-06', '2016-12-10')

# South America Mw 7.75
plot_dart_series(dart_series, '2007-11-13', '2007-11-16')

# West Papua
plot_dart_series(dart_series, '2009-01-01', '2009-01-05')

# Tohoku
plot_dart_series(dart_series, '2011-03-09', '2011-03-15')

# Sumatra 2010
plot_dart_series(dart_series, '2010-04-05', '2010-04-09')

# Mentawai 2010
plot_dart_series(dart_series, '2010-10-23', '2010-10-27')

# Andaman strike-slip 2010
plot_dart_series(dart_series, '2010-06-10', '2010-06-14')

# Chile 2010
plot_dart_series(dart_series, '2010-02-26', '2010-03-01')

# Central america
plot_dart_series(dart_series, '2012-09-04', '2012-09-08')

# Tonga
plot_dart_series(dart_series, '2009-03-18', '2009-03-22')

# North Vanuatu
plot_dart_series(dart_series, '2009-10-06', '2009-10-10')

# Samoa 2009
plot_dart_series(dart_series, '2009-09-27', '2009-10-01')

# Andaman - Arakan 2009
plot_dart_series(dart_series, '2009-08-08', '2009-08-12')

# Puysegur-- one dart is missing from the spatial database and doesn't
# show up here, although we have the data
plot_dart_series(dart_series, '2009-07-14', '2009-07-17')

# Tonga 2006 -- one dart is missing from the spatial database and doesn't
# show up here, although we have the data (but apparently at lower temporal
# resolution than BOM have)
plot_dart_series(dart_series, '2006-05-01', '2006-05-06')

# Sumatra 2007
plot_dart_series(dart_series, '2007-09-10', '2007-09-14')

# South america 2014
plot_dart_series(dart_series, '2014-03-28', '2014-04-04')

# South america 2015
plot_dart_series(dart_series, '2015-09-15', '2015-09-20')

# Java tsunami eq
plot_dart_series(dart_series, '2006-07-15', '2006-07-20')

# Cascadia -- minor tsunami from this.
plot_dart_series(dart_series, '2012-10-26', '2012-10-30')

# Solomons 2007 Mw 8.1
plot_dart_series(dart_series, '2007-03-30', '2007-04-4')

# Solomons 2014 Mw 7.5
plot_dart_series(dart_series, '2014-04-10', '2014-04-14')

# KermadecTonga outer-rise Mw 7.6
plot_dart_series(dart_series, '2011-07-04', '2011-07-09')

# South-america 2007 Mw 8.0
plot_dart_series(dart_series, '2007-08-13', '2007-08-18')

