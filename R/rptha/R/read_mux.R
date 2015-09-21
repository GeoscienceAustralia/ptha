#' Read mux2 data format
#'
#' S4 class to hold data from mux2 files (an output format of the URS tsunami propagation solver)
#' 
#' @param mux2file character. Filename of a mux2 file which is to be read
#' @return An object of class mux2data with slots mux2file, loc, t, wave, wave_tail which
#' hold the data in the file
#' @export
#' @examples
#' \dontrun{
#' # Hypothetical system file
#' filename = "mux2_data_file.mux2"
#' # Read into an object 'x', access with e.g. x@@wave, x@@t, etc
#' x = mux2data(filename)
#' slotNames(x)
#' # Example of accessing the table holding the point location information
#' summary(x@@loc)
#' }
#'
mux2data=setClass('mux2data', 
                  slots=c(mux2file='character', 
                          loc="data.frame", 
                          t='vector', 
                          wave="matrix",
                          wave_tail='vector'))

setMethod(
    'initialize', 'mux2data',
     function(.Object, mux2file = 'character'){        

        # Read mux2data into mux2data object
        .Object@mux2file = mux2file
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

            mm = readBin(tgscon, what='int', n=1, size=4)
            nt[i] = mm
            mm = readBin(tgscon, what='char', n=16, size=1)
            id[i] = paste(mm, sep="", collapse="")
        }

        # Find points inside the grids
        in_grids = (ig>=0)

        .Object@loc = data.frame(geolat=geolat[in_grids],
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
            id=id[in_grids])

        # Read the timeseries data
        nsta2 = sum(in_grids)
        wave = matrix(NA, ncol=nt[1]+1, nrow=nsta2)
        wavet = rep(NA, nt[1]+1)
        .Object@wave_tail = NA
        for(i in 1:(nt[1]+1)){
            # For files with hazard points outside the domain, this seems
            # necessary
            if(i %in% c(1,2)){
                # First 2 rows have nsta points?
                mm = readBin(tgscon, what='double', n=nsta, size=4)
                mm = c(NA, mm[in_grids])
            }else{
                # Remaining rows have nsta2+1 points (includes a time
                # value)
                mm = readBin(tgscon, what='double', n=nsta2+1, size=4)
            }

            if(length(mm) == (nsta2+1)){
                wave[,i] = mm[2:(nsta2+1)]
                wavet[i] = mm[1]
            }else{
                print('Warning: length of wave data not as expected')
                .Object@wave_tail = mm
                break
            }
        }
        .Object@wave = wave
        .Object@t = wavet

        close(tgscon)
        return(.Object)
    }
)


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
