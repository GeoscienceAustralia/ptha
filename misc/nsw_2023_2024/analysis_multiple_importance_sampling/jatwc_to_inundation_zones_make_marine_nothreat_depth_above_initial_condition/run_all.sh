
ATWS_ZONES_ALL=(\
"Eden-Coast" \
"Batemans-Coast" \
"Illawarra-Coast" \
"Sydney-Coast" \
"Hunter-Coast" \
"Macquarie-Coast" \
"Coffs-Coast" \
"Byron-Coast" \
"Lord-Howe-Island" \
"Norfolk-Island" \
)

# Do the calculations for each zone.
for ATWS_ZONE in "${ATWS_ZONES_ALL[@]}"; do
    Rscript make_depth_above_initial_condition.R $ATWS_ZONE no_threat;
    Rscript make_depth_above_initial_condition.R $ATWS_ZONE marine_warning;
    done
