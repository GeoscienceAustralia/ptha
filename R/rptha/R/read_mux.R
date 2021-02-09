#' Read mux2 data format
#'
#' Returns a list with data from mux2 files (an output file format of the URS
#' tsunami propagation solver). This format binary format consists of \cr
#' 1) A 4byte integer witht the number of stations \cr
#' 2) A (large) table containing location and grid information for all stations.
#' This includes a flag 'in_grids' which is -1 if the gauge is not in the grid
#' (in which case its detailed timeseries output is not recorded) \cr
#' 3) A row of 1's (one for each station) \cr
#' 4) A row with -1 and then the number of timesteps for each station \cr
#' 5) A large array with each column giving the stage at each gauge, but the first
#' column giving the time. Only gauges with ig > -1 are recorded. \cr
#' 
#' @param mux2file character. Filename of a mux2 file which is to be read
#' @param inds optional integer vector of indicies of stations to keep. This can be 
#' useful if you want to work on the stations in chunks (e.g. to save memory)
#' @param return_nstations_only logical. If TRUE, then return an integer giving
#' the number of stations in the file FOR WHICH TIMESERIES DATA IS RECORDED.
#' This might not be the same as the number of stations, because of the 'in_grids'
#' treatment.
#' @param file_subtype string with either 'GAR15' or 'PTHA2008'. The former corresponds
#' to a version of the mux2 file from 2015, and the latter corresponds to a version that
#' was used in the 2008 Australian PTHA (at least at some time). The mux2 file format
#' has been varied over time, so it is possible that neither of these works
#' @return A list with entries mux2file, loc, t, wave, wave_tail which
#' hold the data from the file
#' @export
#' @examples
#' \dontrun{
#' # Hypothetical system file
#' filename = "mux2_data_file.mux2"
#' # Read into a list 'x', access with e.g. x$wave, x$t, etc
#' x = read_mux2_data(filename)
#' names(x)
#' # Example of accessing the table holding the point location information
#' summary(x$loc)
#' }
#'
read_mux2_data<-function(mux2file, inds=NULL, return_nstations_only=FALSE,
    file_subtype='GAR15'){

    # Read mux2data into mux2data object
    tgscon = file(mux2file,'rb')

    ## Relevant section of C code
    # fwrite(&nsta,sizeof(int),1,fp);     
    ## Read the number of stations
    nsta = readBin(tgscon, what='int', n=1,size=4)

    ## Read a tgsrwg struct
    # Relevant section of C code
    #struct tgsrwg
    #   {
    #   float geolat;
    #   float geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
    #   float mcolat;
    #   float mcolon;
    #   int ig;
    #   int ilon;     /* grid point location */
    #   int ilat;
    #   float z;             /* water depth at this location */
    #   float center_lat, center_lon;        /**** center of this array *****/
    #   float offset,az,baz;                 /* for arrays this is the distance in km from center of array */
    #   float dt;            /* sampling rate for this site */
    #   int nt;              /* number of points */
    #   char id[16];         /* identifier */
    #   };
    #
    #

    # Read location information
    # These variables will be in the data.frame
    geolat = rep(NA,nsta)
    geolong = rep(NA,nsta)
    mclat = rep(NA,nsta)
    mclong = rep(NA,nsta)
    ig = rep(NA,nsta)
    ilon = rep(NA,nsta)
    ilat = rep(NA,nsta)
    z = rep(NA,nsta)
    center_lat = rep(NA,nsta)
    center_lon = rep(NA,nsta)
    offset = rep(NA,nsta)
    az = rep(NA,nsta)
    baz = rep(NA,nsta)
    nt = rep(NA,nsta)
    dt = rep(NA,nsta)
    id = rep(NA,nsta)

    for(i in 1:nsta){
        mm = readBin(tgscon, what='double', n=4, size=4)
        geolat[i] = mm[1]
        geolong[i] = mm[2]
        mclat[i] = mm[3]
        mclong[i] = mm[4]

        mm = readBin(tgscon, what='int', n=3, size=4)
        ig[i] = mm[1]
        ilon[i] = mm[2]
        ilat[i] = mm[3]

        mm = readBin(tgscon, what='double', n=7, size=4)
        z[i] = mm[1]
        center_lat[i] = mm[2]
        center_lon[i] = mm[3]
        offset[i] = mm[4]
        az[i] = mm[5]
        baz[i] = mm[6]
        dt[i] = mm[7]

        mm = readBin(tgscon, what='int', n=1, size=4)
        nt[i] = mm
        mm = readBin(tgscon, what='char', n=16, size=1)
        id[i] = paste(mm, sep="", collapse="")
    }


    # At this point out treatment of the file types varies
    if(file_subtype == 'GAR15'){

        # Find points inside the grids
        in_grids = (ig >= 0)
        
        if(return_nstations_only){
            close(tgscon)
            return(sum(in_grids))
        }

        # Key location summary statistics
        loc = data.frame(
            geolat=geolat[in_grids],
            geolong=geolong[in_grids],
            mclat=mclat[in_grids],
            mclong=mclong[in_grids],
            ig=ig[in_grids],
            ilon=ilon[in_grids],
            ilat=ilat[in_grids],
            z=z[in_grids],
            center_lat=center_lat[in_grids],
            center_lon=center_lon[in_grids],
            offset=offset[in_grids],
            az=az[in_grids],
            baz=baz[in_grids],
            dt=dt[in_grids],
            id=id[in_grids])

        # Read the timeseries data
        nsta2 = sum(in_grids)
        wave = matrix(NA, ncol=nt[1], nrow=nsta2)
        wavet = rep(NA, nt[1])
        wave_tail = NA

        for(i in 1:(nt[1]+2)){
            # For files with hazard points outside the domain, this seems
            # necessary
            if(i %in% c(1,2)){
                # First 2 rows have nsta points, and basically contain header type
                # information
                mm = readBin(tgscon, what='double', n=nsta, size=4)
                next
            }else{
                # Remaining rows have nsta2+1 points (includes a time
                # value)
                mm = readBin(tgscon, what='double', n=nsta2+1, size=4)
            }

            if( (i > 2) & (length(mm) == (nsta2+1))){
                wavet[i-2] = mm[1]
                wave[,i-2] = mm[-1]
            }else{
                print('Warning: length of wave data not as expected')
                wave_tail = mm
                break
            }
        }

        close(tgscon)

        if(is.null(inds)){
            output = list(mux2file=mux2file, loc=loc, t=wavet, wave=wave, wave_tail=wave_tail)
        }else{

            if( min(inds) <= 0 ) stop('inds must all be > 0')
            if( max(inds) > length(loc[,1]) ) stop('inds must all be <= number of stations')

            output = list(mux2file=mux2file, loc=loc[inds,], t=wavet, 
                wave=wave[inds,], wave_tail=wave_tail)
        }
    
        # END of specific treatment for GAR15 file_subtype

    }else if(file_subtype == 'PTHA2008'){


        # Find points inside the grids
        #in_grids = (ig >= 0)
        in_grids = 1:length(geolat)
        
        if(return_nstations_only){
            close(tgscon)
            return(sum(in_grids))
        }

        # Key location summary statistics
        loc = data.frame(
            geolat=geolat,
            geolong=geolong,
            mclat=mclat,
            mclong=mclong,
            ig=ig,
            ilon=ilon,
            ilat=ilat,
            z=z,
            center_lat=center_lat,
            center_lon=center_lon,
            offset=offset,
            az=az,
            baz=baz,
            dt=dt,
            id=id)

        # Read the timeseries data

        # Recording for each station 'starts' and 'ends' at a particular time
        t_start = readBin(tgscon, what='int', n=nsta, size=4)
        t_end = readBin(tgscon, what='int', n=nsta, size=4)

        # Now read the remainder of the file. Here 'n' will be too large, unless every
        # gauge recorded every event
        binary_dump = readBin(tgscon, what='double', n=nsta*max(t_end), size=4)

        wave = matrix(0, ncol=max(t_end), nrow=nsta)
        wavet = rep(NA, max(t_end))

        counter = 0
        for(i in 1:length(wavet)){
            counter = counter+1
            local_t = binary_dump[counter] #readBin(tgscon, what='double', n=1, size=4)
            if(length(local_t) == 0){
                stop('READING ERROR, length(local_t) == 0')
            }
            wavet[i] = local_t
            kk = which(t_start <= i & t_end >= i)
            if(length(kk) > 0){
                local_wave = binary_dump[counter + (1:length(kk))] #readBin(tgscon, what='double', n=length(kk), size=4)
                wave[kk,i] = local_wave 
                counter = counter+length(kk)
            }
        }

        close(tgscon)

        # If we correctly parsed the file, this should be true
        stopifnot(counter == length(binary_dump))


        # Output
        if(is.null(inds)){
            output = list(mux2file=mux2file, loc=loc, t=wavet, wave=wave)
        }else{

            if( min(inds) <= 0 ) stop('inds must all be > 0')
            if( max(inds) > length(loc[,1]) ) stop('inds must all be <= number of stations')

            output = list(mux2file=mux2file, loc=loc[inds,], t=wavet, 
                wave=wave[inds,])
        }


    }

    return(output)

}

#/*** DEFINITION FOR TIDE GAUGE OUTPUT ***/
#struct tgs
#   {
#   int ista;
#   float geolat,geolon; /* location and water depth in geographic coordinates -180/180/-90/90 */
#   int   ilat,ilon;     /* grid point location */
#
#   float z;             /* water depth at this location */
#
#   float center_lat, center_lon;        /**** center of this array *****/
#   float offset,az,baz;                 /* for arrays this is the distance in km from center of array */
#
#   float dt;            /* sampling rate for this site */
#   int nt;              /* number of points */
#
#   char id[16];         /* identifier */
#   };




#' read a mux2 data file
#'
#' This function does basically the same thing as the previous read_mux2_data
#' file (for GAR15 file_subtype only). It is arguably more elegant. 
# I had hoped it would be faster but it does not seem to be
#' significantly faster. Anyway it might be useful for debugging so I keep it
#' here, but am not exporting it into the package namespace
#'
#' @param mux2file mux2 filename
#' @param inds indices to extract, see \code{read_mux2_data} for details
#' @return see \code{read_mux2_data} for details
read_mux2_data_alternative<-function(mux2file, inds=NULL){

    desired_inds = inds
    tgscon = file(mux2file,'rb')

    ## Relevant section of C code
    # fwrite(&nsta,sizeof(int),1,fp);     
    ## Read the number of stations
    nsta = readBin(tgscon, what='int', n=1, size=4)

    # There are the equivalent of 19 x 4byte columns in the table, with one row for each station
    # Read first a raw binary, and then work on each column separately
    table_binary_stream = readBin(tgscon, what = 'raw', n = nsta*19*4, size=1)

    # Trick to extract indices corresponding to 'first column data' bytes
    bigseq = seq(0, 19*4*nsta, 19*4)
    bigMat = matrix(NA, nrow=length(bigseq), ncol=4)
    for(i in 1:ncol(bigMat)) bigMat[,i] = bigseq + i
    inds = c(t(bigMat))

    # Read each column
    geolat = readBin(table_binary_stream[inds], what='double', n=nsta, size=4)
    geolong = readBin(table_binary_stream[inds+4], what='double', n=nsta, size=4)
    mclat = readBin(table_binary_stream[inds+8], what='double', n=nsta, size=4)
    mclong = readBin(table_binary_stream[inds+12], what='double', n=nsta, size=4)
    ig = readBin(table_binary_stream[inds + 16], what = 'int', n=nsta, size=4)
    ilon = readBin(table_binary_stream[inds + 20], what = 'int', n=nsta, size=4)
    ilat = readBin(table_binary_stream[inds + 24], what = 'int', n=nsta, size=4)
    z = readBin(table_binary_stream[inds + 28], what = 'double', n=nsta, size=4)
    center_lat = readBin(table_binary_stream[inds + 32], what = 'double', n=nsta, size=4)
    center_lon = readBin(table_binary_stream[inds + 36], what = 'double', n=nsta, size=4)
    offset = readBin(table_binary_stream[inds + 40], what = 'double', n=nsta, size=4)
    az = readBin(table_binary_stream[inds + 44], what = 'double', n=nsta, size=4)
    baz = readBin(table_binary_stream[inds + 48], what = 'double', n=nsta, size=4)
    dt = readBin(table_binary_stream[inds + 52], what = 'double', n=nsta, size=4)
    nt = readBin(table_binary_stream[inds + 56], what = 'int', n=nsta, size=4)
    # Final character
    id = matrix('', ncol=16,nrow=nsta)
    for(i in 1:16){
        id[,i] = readBin(table_binary_stream[bigMat[,1] + 60 + i-1], what = 'char', n=nsta, size=1)
    }
    id = apply(id, 1, function(x) paste(x, sep="", collapse=""))

    # Find points inside the grids
    in_grids = (ig >= 0)
        
    # Key location summary statistics
    loc = data.frame(
        geolat=geolat[in_grids],
        geolong=geolong[in_grids],
        mclat=mclat[in_grids],
        mclong=mclong[in_grids],
        ig=ig[in_grids],
        ilon=ilon[in_grids],
        ilat=ilat[in_grids],
        z=z[in_grids],
        center_lat=center_lat[in_grids],
        center_lon=center_lon[in_grids],
        offset=offset[in_grids],
        az=az[in_grids],
        baz=baz[in_grids],
        dt=dt[in_grids],
        id=id[in_grids])

    ## Now read the wave timeseries at each gauge ##
    nsta2 = sum(in_grids)
    wave = matrix(NA, ncol=nt[1], nrow=nsta2)
    wavet = rep(NA, nt[1])
    wave_tail = NA

    # Read a line of ones and negative ones showing which stations are included
    ones = readBin(tgscon, what = 'integer', n = nsta, size=4)
    # Read a line that contains the number of time-steps (nt[1]) in every entry, except
    # those where no data are recorded, where it has (-1)
    another_header = readBin(tgscon, what = 'integer', n = nsta, size=4)

    if(nt[1] > 1){
        wave_double_stream = readBin(tgscon, what = 'double', n = (nsta2+1)*(nt[1]), size=4)

        if(length(wave_double_stream)!= (nsta2+1)*nt[1]){
            print('Warning: Possible parsing error reading file')
        }

        wave_time_indices =  seq(0, (nsta2+1)*(nt[1]-1), by=(nsta2 + 1)) + 1

        next_int_stream = readBin(tgscon, what = 'integer', n = 2*nsta, size=4)

        wavet = wave_double_stream[wave_time_indices]
        wave[,1:nt[1]] = wave_double_stream[-wave_time_indices]
    }


    close(tgscon) 
    # Return output
    if(is.null(desired_inds)){
        output = list(mux2file=mux2file, loc=loc, t=wavet, wave=wave, wave_tail=wave_tail)

    }else{
        if( min(desired_inds) <= 0 ) stop('inds must all be > 0')
        if( max(desired_inds) > length(loc[,1]) ){
            stop('inds must all be <= number of stations')
        }

        output = list(mux2file=mux2file, loc=loc[desired_inds,], t=wavet, 
            wave=wave[desired_inds,], wave_tail=NA)
    }

    return(output)

}


# Compute the zero-crossing-period of x
#
# This function computes both the up-crossing and down-crossing periods and
# returns their average. No interpolation is used. Consider \code{gauge_statistics_simple}
# for a faster alternative that also uses interpolation (so is more accurate).
#
# @param x numeric vector containing a timeseries with mean (approx) zero
# @param dt The time between consecutive samples of x
# @return The zero crossing period of x
# 
# @examples
# x = seq(1,2000, by=2)
# period = 30
# y = sin(x*2*pi/period)
# # Check we can compute this
# y_period = zero_crossing_period(y, dt=2)
# stopifnot(abs(y_period/period - 1) < 0.1)
#
zero_crossing_period <-function(x, dt=1){
    sg_x = sign(x)
    n = length(sg_x)
    # Get 'positive' zero crossings
    up_cross = which((diff(sg_x) > 0) & (sg_x[2:n] > 0.0))
    # Get 'negative' zero crossings
    down_cross = which((diff(sg_x) < 0) & (sg_x[2:n] < 0.0))

    lu = length(up_cross)
    ld = length(down_cross) 

    if( lu < 2 | ld < 2){
        #warning('To few zero crossings for period computation')
        return(NA)
    }

    # Compute periods
    up_period = (up_cross[lu] - up_cross[1])*dt/(lu - 1)
    down_period = (down_cross[ld] - down_cross[1])*dt/(ld - 1)

    return(0.5*(up_period + down_period))
}
