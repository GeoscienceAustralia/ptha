#
# Copied from PTHA codes to get events similar to DART buoys
#

find_unit_sources_near_hypocentre<-function(
    event_hypocentre,
    unit_source_geometry,
    unit_source_statistics,
    event_magnitude,
    scaling_relation_type='Strasser'){

    # Find which events contain the hypocentre, by finding which unit source
    # contains it
    unit_source_containing_hypocentre = find_unit_source_index_containing_point(
        event_hypocentre, unit_source_geometry, unit_source_statistics) 

    # Allow events which 'touch' sites within uniform slip scaling law width
    # and half length
    local_AWL = Mw_2_rupture_size(event_magnitude, 
        relation=scaling_relation_type)
    expand_unit_source_alongstrike = ceiling(local_AWL[3] * 0.5/
        unit_source_statistics$length[unit_source_containing_hypocentre])
    expand_unit_source_downdip = ceiling(local_AWL[2] * 0.5/
        unit_source_statistics$width[unit_source_containing_hypocentre])

    # Find unit-source neighbours (within a few unit-sources in each direction)
    hypocentre_neighbours = c()
    for(j in (seq(-expand_unit_source_downdip,expand_unit_source_downdip))){
        for(i in (seq(-expand_unit_source_alongstrike,expand_unit_source_alongstrike))){
            usch = unit_source_containing_hypocentre # shorthand
            nbr = which(
                (unit_source_statistics$downdip_number == 
                    (unit_source_statistics$downdip_number[usch] + j)) &
                (unit_source_statistics$alongstrike_number == 
                    (unit_source_statistics$alongstrike_number[usch] + i))
                )
            if(length(nbr) > 1) stop('BUG! This should be impossible')
            if(length(nbr) == 1){
                hypocentre_neighbours = c(hypocentre_neighbours, nbr)
            }
        }
    }
    if(length(hypocentre_neighbours) == 0) stop('No unit source neighbours found')

    output = c(unit_source_containing_hypocentre, hypocentre_neighbours)

    return(output)

}

