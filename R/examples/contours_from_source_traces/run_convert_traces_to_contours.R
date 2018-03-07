# Code to make 'geometrically-defined' subduction interface contours
# These are used to assist manual editing of SLAB1.0 contours (where
# the latter do not extend as far as we would like), and also serve
# as contours in sites without any SLAB data
#
# Gareth Davies, Geoscience Australia 2018
#
# The algorithm takes source traces (defining the trench), with particular
# attributes
#
# BEWARE THE CONTOURS LIKELY NEED POST-HOC EDITING -- THIS CODE PROVIDES A
# STARTING POINT ONLY
#

source_traces = c('tanimbar/tanimbar.shp', 'timor/timor.shp', 'tolo_thrust/tolo_thrust.shp')

#
# Directory where output shapefile will be saved
#
output_dir = 'CONTOURS'

############################################################################
#
# Generate corresponding contours for the fault traces 
#
############################################################################

library(rptha)
# Key computational routines
source('collate_source_traces.R') 
#
# Numerical parameter -- contours are extended during calculation to avoid
# bottom-edge artefacts. You probably do not need to change this number
#
maxdepth_offset = 50 


dir.create(output_dir, showWarnings=FALSE)

store_geom_cont = list()
for(i in 1:length(source_traces)){

    # Read fault trace
    mylayer = strsplit(basename(source_traces[i]), ".shp")[[1]][1]
    mytrace = readOGR(dsn=dirname(source_traces[i]),
                    layer=mylayer, stringsAsFactors=FALSE)

    # Use the trace and the known min/max depth to extend it to a rupture surface 
    # The method is generates perpendicular lines down-dip from the contour,
    # with the correct distance to match the 'max-depth' given the dip. It
    # turns these into a polygon, and interpolates depth in this polygon
    # linearly (or parabolically) between the min/max based on the distance
    # down dip, and the provided dip

    #
    # PARSE SHAPEFILE INPUTS, based on format of files in these directories.
    #
    # We assume the input shapefile contains the max_depth (constant), the
    # concavity ('Up', 'Down', 'Linear'), the Dip (for 'Linear' concavity)
    # or values Dip_0 (giving the dip at the trench) and Dip_XX giving the dip
    # at depth XX (e.g. Dip_25). 

    # Parse max depth
    if(!('Max_depth' %in% names(mytrace))){
        stop(paste0('file ', source_traces[i], ' does not have a Max_depth attribute'))
    }
    maxdepth = as.numeric(mytrace$Max_depth)
    # max depth should be constant
    if(!all(maxdepth == maxdepth[1])){
        stop(paste0('file ', source_traces[i], 
            ' has varying Max_depth. It needs to be the same everywhere'))
    }
    # Ensure it is a constant
    maxdepth = maxdepth[1]

    # Parse concavity. Only 3 values are allowed. Depending on the concavity, the treatment
    # of Dip will vary
    if(!all(mytrace$Concavity %in% c('Up', 'Down', 'Linear'))){
        stop(paste0('file ', source_traces[i], 
            ' does not have a Concavity attribute %in% c("Up", "Down", "Linear") '))
    }

    # Check concavity is constant
    if(!all(mytrace$Concavity == mytrace$Concavity[1])){
        stop(paste0('file ', source_traces[i], 
            ' does not have the same Concavity everywhere. The code', 
            ' requires constant Concavity'))
    }

    # Extract dip, depending on how the concavity
    if(mytrace$Concavity[1] == 'Linear'){
        #
        # There can only be one dip
        #
        dip_cols = grep('Dip', substring(names(mytrace), 1, 3))
        if(length(dip_cols) != 1){
            stop(paste0('file ', source_traces[i], 
            ' has more than one attribute starting with "Dip", but also has ', 
            '"Linear" Concavity, which only allows one dip value. This is not allowed.')) 
        }else{
            # Put a 'Dip' column in the SpatialLinesDataFrame
            mytrace$Dip = as.numeric(mytrace@data[ , dip_cols])
        }

        # Other variables we need for this case
        downdip_profile = 'linear'
        dip_depth = NULL
        mindepth = 0
        mindepth_dip = NULL
            
    }else{
        #
        # Concavity = 'Up' or 'Down'
        #
        # Since concavity is nonlinear, the dip columns need to include a
        # '_depth', e.g. Dip_0, Dip_10, ...
        #
        dip_cols = grep('Dip_', substring(names(mytrace), 1, 4))
        if(length(dip_cols) != 2) stop(paste0('file ', source_traces[i], 
            ' has non-linear concavity, so must have exactly 2 attributes starting with',
            ' "Dip_" followed by the depth at which the dip applies. But it does not.'))

        # Extract depth from name
        dip_depths = as.numeric(sapply(names(mytrace)[dip_cols], 
            f<-function(x) as.numeric(strsplit(x, '_')[[1]][2])))

        # Dip_depths must include 0, and be less than maxdepth
        if(!(min(dip_depths) == 0)) stop('file ', source_traces[i], 
            ' does not have Dip_0, which is required for non-linear concavities.')
        if(!(all(dip_depths <= maxdepth)) ){
            stop(paste0('file ', source_traces[i], 
            ' has Dip specified at a depth that is greater than Max_depth.',
            ' This is not allowed'))
        }

        downdip_profile = 'parabolic' # 'linear'

        # Must have rupture to trench
        # Dip at trench
        mindepth = 0 # Constant
        mindepth_dip = as.numeric(mytrace$Dip_0) # Can vary
        # Dip at some depth
        dip_depth = as.numeric(mytrace$Dip_0) * 0 + max(dip_depths) # Fixed by the shapefile format
        mytrace$Dip = as.numeric(mytrace@data[, dip_cols[which.max(dip_depths)]]) # Can vary


        if(mytrace$Concavity[1] == 'Up'){
            if(any(mytrace$Dip > mindepth_dip)){
                stop('file ', source_traces[i], 
            ' is supposed to be concave-up, however, the deep dip is greater than the shallow dip in one or more locations') 
            }
        }

        if(mytrace$Concavity[1] == 'Down'){
            if(any(mytrace$Dip < mindepth_dip)){
                stop('file ', source_traces[i], 
            ' is supposed to be concave-down, however, the deep dip is less than the shallow dip in one or more locations') 
            }
        }
    }
    
    myconts = extend_trace_to_depth_contours(
        mytrace,
        maxdepth=maxdepth,
        mindepth=mindepth,
        maxdepth_offset=maxdepth_offset, 
        downdip_profile=downdip_profile,
        dip_depth=dip_depth,
        mindepth_dip = mindepth_dip)

    store_geom_cont[[i]] = c(myconts, mytrace)

    # Get rid of the contours > (maxdepth + 4.5)
    out_conts = myconts[[3]]
    keep = which( as.numeric(as.character(out_conts@data[,1])) <= 
        (maxdepth + 4.5))
    out_conts = out_conts[keep,]

    writeOGR(out_conts, dsn=output_dir, layer=paste0(mylayer, '_contours'), 
        driver='ESRI Shapefile', overwrite=TRUE)
}

save.image('session_end.Rdata')
