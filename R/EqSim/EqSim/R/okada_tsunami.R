#'
#' Compute ground surface displacements in linear elastic half-space from known
#' subfault ruptures, using Okada's (1985) solution for rectangular subfaults
#' 
#' Complex ruptures can be treated using multiple sub-faults
#'
#'@param elon -- numeric vector with x rupture centroid location (m)
#'@param elat -- numeric vector with y of rupture centroid location (m)
#'@param edep -- numeric vector with centroid depth of rupture (km)
#'@param strk -- numeric vector with strike of sub-fault (degrees clockwise from North)
#'@param dip -- numeric vector with dip of rupture (degrees below the horizontal, dipping to the right when looking in the along-strike direction)
#'@param lnth -- numeric vector with length of sub-fault (km)
#'@param wdt -- numeric vector with width of sub-fault (km)
#'@param disl1 -- numeric vector with along-strike disloacation on the sub-fault (m)
#'@param disl2 -- numeric vector with up-dip disloacation on the sub-fault (m)
#'@param rlon -- numeric vector with x locations where output is desired (m)
#'@param rlat -- numeric vector with y locations where output is desired (m)
#'@param dstmx -- optional maximum distance at which sub-faults can cause displacement 
#'@param verbose -- TRUE/FALSE -- print info on ground deformation
#'
#'@return -- list with edsp, ndsp, zdsp giving the displacements in the
#' east,north, and vertical directions. Often we are only interested in the
#' latter.
#'@export
#'@useDynLib EqSim
#'@examples
#' # Simple example -- pure thrust fault
#' mypts=expand.grid(seq(-100.5,100.5,len=100), seq(-100.5,100.5,len=100))*1000
#' ff=okada_tsunami(0,0,30,0,15,1,1,0,1,
#'                 mypts[,1],mypts[,2])
#' plot(mypts[,1], mypts[,2],col=heat.colors(10)[cut(ff[[3]],10)],pch=19,asp=1)
#' 
#' ## More complex example
#' ff=okada_tsunami(c(0,0),c(0,50000.),c(30,25),c(0,0),c(15,12),c(1,1),c(1,1),c(0,0),c(1,1),
#'                  mypts[,1],mypts[,2])
#' 
#' plot(mypts[,1], mypts[,2],col=heat.colors(10)[cut(ff[[3]],10)],pch=19,asp=1)
okada_tsunami<-function(elon,elat,edep,strk,dip,lnth,wdt,
                        disl1,disl2,rlon,rlat,dstmx=9.0e+12,
                        verbose=FALSE){
    # Call Okada routines to compute surface displacement

    # Variables we already know
    alp=0.5
    n=length(elon)
    m=length(rlon)

    # Check input
    if(length(elat)!=n | length(edep)!=n | length(strk)!=n | length(dip)!=n |
       length(lnth)!=n | length(wdt)!=n  | length(disl1)!=n | length(disl2)!=n ){
       stop('sub-fault variables are not all of the same length ')
    }
    if(length(rlat)!=m) stop('length(rlon)!=length(rlat)')

    # Output variables
    #edsp=numeric(m)
    #ndsp=numeric(m)
    #zdsp=numeric(m)
    edsp=rep(1.0e-16,m)
    ndsp=rep(1.0e-16,m)
    zdsp=rep(1.0e-16,m)

    # Fortran call:
    #fault_disp(alp,elon,elat,edep,strk,dip,length,wdt,
    # &     disl1,disl2,rlon,rlat,dstmx,edsp,ndsp,zdsp,m,n)
    xout=.Fortran('fault_disp', as.double(alp),as.double(elon),as.double(elat),as.double(edep),
             as.double(strk),as.double(dip),as.double(lnth),as.double(wdt),
             as.double(disl1),as.double(disl2),as.double(rlon),as.double(rlat),as.double(dstmx),
             as.double(edsp),as.double(ndsp),as.double(zdsp),as.integer(m),as.integer(n), 
             DUP=TRUE, # DUP=FALSE is deprecated
             PACKAGE='EqSim')

    names(xout)[14:16]=c('edsp','ndsp','zdsp')
    if(verbose){
        print(paste('Max zdsp: ', max(xout$zdsp)))
        print(paste('Min zdsp: ', min(xout$zdsp)))
    }
    return(xout[14:16])
}

