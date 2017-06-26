#
# Find the 'Bird' convergence rates along our unit-source top edges
#
suppressPackageStartupMessages(library(rptha))


#' Map Bird (2003) convergence data onto our unit-sources, and 
#
#' # Read Bird's data and setup key information required to make conditional
#' # probability functions
#' make_conditional_probability_function = event_conditional_probability_factory()
#' 
#' # Make the function for puysegur
#' source_name = 'puysegur' 
#' puysegur_conditional_prob_function = make_conditional_probability_function(source_name)
#'
#' # Suppose we have read the puysegur_uniform_slip_events table into a data.frame
#' Mw_8.0_events = puysegur_uniform_slip_events[which(puysegur_uniform_slip_events$Mw == 8.0),]
#' # Get the conditional probability of these "Mw = 8.0" events as:
#' Mw_8.0_conditional_prob = puysegur_conditional_prob_function(Mw_8.0_events)
#' 
event_conditional_probability_factory<-function(return_environment=FALSE){

    #
    # Parse bird's data (which is zipped for compression)
    # 

    bird_data = '../DATA/BIRD_PLATE_BOUNDARIES/PB2002_steps.dat.txt.zip'
    bd = read.table(
            unz(description=bird_data, 
                filename=basename(gsub('\\.zip', '', bird_data))) 
        )

    # See table 2 of the bird paper for definitions of the data columns
    names(bd) = c('id', 'plateboundary', 'lon1', 'lat1', 'lon2', 'lat2', 
        'length', 'azi', 'vel_L2R', 'vel_azi', 'vel_div', 'vel_rl', 
        'elev', 'age', 'class')
    bird_centroid = midPoint(
        as.matrix(bd[,c('lon1', 'lat1')]), 
        as.matrix(bd[,c('lon2', 'lat2')]), 
        f = 0)

    #
    # Make bird data as a SpatialLinesDataFrame, and save it for QC
    #
    lines_list = vector(mode='list', length=length(bd[,1]))
    for(i in 1:length(bd[,1])){

        line_mat = matrix(c(bd$lon1[i], bd$lat1[i], bd$lon2[i], bd$lat2[i]), ncol=2, byrow=T)

        line_mat[2,] = adjust_longitude_by_360_deg(line_mat[2,], line_mat[1,])

        lines_list[[i]] = Lines(list(Line(line_mat)), ID=as.character(i))
    }

    bd_sl = SpatialLines(lines_list, proj4string=CRS("+init=epsg:4326"))
    bd_sldf = SpatialLinesDataFrame(bd_sl, data=bd, match=FALSE)

    writeOGR(bd_sldf, dsn='bird_2003', layer='bird_2003', driver='ESRI Shapefile', overwrite=TRUE)


    #
    # Parse unit-source top edges
    #
    unit_source_files = Sys.glob(
        '../SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source_statistics*.nc')
    unit_source_tables = lapply(as.list(unit_source_files), read_table_from_netcdf)
    names(unit_source_tables) = unit_source_files

    #
    # Make 'top-edge-only' tables
    #
    top_edge_tables = unit_source_tables
    for(i in 1:length(unit_source_tables)){
        dd = top_edge_tables[[i]]$downdip_number
        top_edge_tables[[i]] = top_edge_tables[[i]][which(dd==1),]
    }

    #
    # Find Bird centroid nearest to top_edges
    #
    nearest_bird_point<-function(p){
        p_mat = bird_centroid*0
        p_mat[,1] = as.numeric(p[1])
        p_mat[,2] = as.numeric(p[2])

        # Bypass warnings about longitudes > 180 [which does not cause problems]
        suppressWarnings({ distances =  distHaversine(p_mat, bird_centroid) })
        k = which.min(distances)
        output = c(k, distances[k])
        return(output)
    }

    # For each table, loop over all 'trench' unit sources and find the nearest bird centroid to
    # the top edge. Make a SpatialLines object as we go for QC
    sldf_list = vector(mode='list', length=length(top_edge_tables))
    for(i in 1:length(top_edge_tables)){

        ti = top_edge_tables[[i]]
        di = ti[,1]*0 # Store distances to nearest bird point
        ki = ti[,1]*0 # Store index of nearest bird point

        lines_list = vector(mode='list', length=length(di))
        for(j in 1:nrow(ti)){

            top_point = as.numeric(ti[j,1:2])
            # Convert 'top_point' to the 'top-edge' point
            half_width_surface = as.numeric(ti$width[j])/2 * 1000 * cos(2*pi*ti$dip[j]/180)
            top_point_top_edge_approx = destPoint(top_point, ti$strike[j]-90, half_width_surface)

            output = nearest_bird_point(top_point_top_edge_approx)
            di[j] = output[2]
            ki[j] = output[1]

            # Prepare SpatialLines output, mapping the unit-source top edge to the Bird lines
            line_mat = matrix(c(bird_centroid[ki[j],1:2], top_point_top_edge_approx[1:2]), 
                nrow=2,byrow=TRUE)
            line_mat[2,] = adjust_longitude_by_360_deg(line_mat[2,], line_mat[1,])
            lines_list[[i]] = Lines(list(Line(line_mat)), ID=as.character(j))
        }

        top_edge_tables[[i]] = cbind(ti, 
            data.frame(
                'distance_bird' = di, 
                'bird_index' = ki, 
                'bird_vel_div' = bd$vel_div[ki], 
                'bird_vel_rl' = bd$vel_rl[ki])
            )
        sl = SpatialLines(lines_list, proj4string=CRS("+init=epsg:4326"))
    
        sldf_list[[i]] = SpatialLinesDataFrame(sl, data=top_edge_tables[[i]])
    }

    #
    # Use the top-edge data to populate the unit-source tables
    #
    for(i in 1:length(top_edge_tables)){

        ui = unit_source_tables[[i]]
        ti = top_edge_tables[[i]]

        kk = match(ui$alongstrike_number, ti$alongstrike_number)
        ui = cbind(ui, ti[kk, c('distance_bird', 'bird_index', 'bird_vel_div', 'bird_vel_rl')] )
        # The above cbind mangles the rownames, so fix that here.
        rownames(ui) = rownames(unit_source_tables[[i]])

        unit_source_tables[[i]] = ui
    }


    #'
    #' Create a function which gives the conditional probability of events with
    #' a fixed Mw on the chosen source-zone, which accounts for spatially variable
    #' slip and source-zone area
    #'
    make_conditional_probability_function_uniform_slip<-function(source_name){

        unit_source_match = grep( source_name, 
            basename(dirname(dirname(names(unit_source_tables)))) )

        if(length(unit_source_match) != 1){
            print(unit_source_match)
            stop( paste0('Could not uniquely match source_name ', source_name) )
        }

        # Get relevant part of 'unit_source_tables' from parent environment
        uss = unit_source_tables[[unit_source_match]]
        dim_uss = dim(uss)
       
        # Ensure table is correctly sorted 
        stopifnot(all(uss$subfault_number == 1:dim_uss[1]))

        #'
        #' Conditional probability function that can be passed to
        #' \code{get_event_probabilities_conditional_on_Mw} .
        #'
        conditional_probability_function<-function(events_with_Mw, debug_output=FALSE){

            if( 'event_slip_string' %in% names(events_with_Mw) ){
                msg = paste0('Tried to pass a stochastic slip events table to ',
                    'conditional_probability_function. Only uniform-slip tables are permitted')
                stop(msg)
            }

            if( !all(events_with_Mw$Mw == events_with_Mw$Mw[1]) ){
                stop('conditional probability function requires input events to have the same Mw')
            }

            # Get the unit sources in the event, and their area-weighted long-term slip [setting
            # non-convergent slip components to zero]
            long_term_slip_near_event = rep(NA, nrow(events_with_Mw))
            for(i in 1:nrow(events_with_Mw)){

                ui = get_unit_source_indices_in_event(events_with_Mw[i,])
                areas = uss$length[ui] * uss$width[ui]
                # Note with bird's data, negative 'vel_div' means convergence
                convergent_slip = pmax(0, -uss$bird_vel_div[ui])
                long_term_slip_near_event[i] = sum(areas * convergent_slip)/sum(areas)
            }
            # Set the conditional proability as proportional to [event area x long-term-slip]
            conditional_probability = (events_with_Mw$area * long_term_slip_near_event)
            conditional_probability = conditional_probability/sum(conditional_probability)

            if(!debug_output){
                # Default case
                return(conditional_probability)

            }else{
                # For debugging, return the information we used to make the table
                output = list(unit_sources = uss, 
                    long_term_slip_near_event = long_term_slip_near_event, 
                    conditional_probability=conditional_probability,
                    events_with_Mw = events_with_Mw)
                return(output)

            }

        }

        return(conditional_probability_function)
    }


    if(!return_environment){
        return(make_conditional_probability_function_uniform_slip)
    }else{
        # For debugging, it's useful to have the entire function environment
        return(environment())
    }

}


