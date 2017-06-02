./ELEV -- Pre-processed elevation datasets

./HAZARD_POINTS -- Hazard point files for modelling

./SOURCEZONE_CONTOURS -- Contours of source-zones, made by manually merging
    SLAB and planar contours. These are used to make the earthquake sources.

./SOURCEZONE_DOWNDIP_LINES -- Lines which traverse the source-zone contours in
    a down-dip direction, and are used to define the lateral boundaries of the
    unit-sources. They were generally created by using the rptha function
    'create_downdip_lines_on_source_contours_improved', and then converted
    to shapefile with 'downdip_lines_to_SpatialLinesDataFrame', and potentially
    then manually edited. (see make_initial_downdip_lines.R for a script which does 
    the non-manual parts).

