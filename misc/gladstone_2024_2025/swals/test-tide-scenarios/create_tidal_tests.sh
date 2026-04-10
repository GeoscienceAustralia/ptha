#!/bin/bash
# This script creates a set of tidal test scenarios for the SWALS tsunami model.
# The script is intended to be run from the directory containing the SWALS source code.
# It will create a set of PBS files that can be submitted to a PBS queueing system.
# It modifies a template PBS file to create the scenarios.


# iterate over locations (and their respective sea levels) and stage files
declare -A locations
locations=(\
["lady_elliot"]=1.41 \
["south_trees"]=2.42 \
["rosslyn_bay"]=2.85)

path_swals=.

declare -A stage_files
stage_files=(\
["southamerica"]=$path_swals/../sources/test/large_southamerica_148059/southamerica_row_0148059_Mw_96_HS.tif \
["large_kermadec"]=$path_swals/../sources/test/large_kermadectonga2_43783/kermadectonga2_row_0043783_Mw_95_HS.tif \
["solomon2"]=$path_swals/../sources/test/large_solomon2_16040/solomon2_row_0016040_Mw_94_HS.tif \
["kermadectonga2"]=$path_swals/../sources/test/medium_kermadectonga2_37695/kermadectonga2_row_0037695_Mw_87_HS.tif)


scenario_pbs=tidal_test_template.pbs
for location in "${!locations[@]}"; do
    sea_level=${locations[$location]}
    for stage_file in "${!stage_files[@]}"; do
        echo "Creating scenario for $location with sea level $sea_level and stage file $stage_file"
        sed -i "s/ambient_sea_level=.*/ambient_sea_level=$sea_level/g" $scenario_pbs
        sed -i "s/location=.*/location=$location/g" $scenario_pbs
        sed -i "s/event=.*/event=$stage_file/g" $scenario_pbs
        # sed -i "s/stage_file=.*/stage_file=${stage_files[$stage_file]}/g" $scenario_pbs
        # sed: -e expression #1, char 29: unknown option to `s'
        sed -i "s|stage_file=.*|stage_file=${stage_files[$stage_file]}|g" $scenario_pbs
        
        # create a new file with the scenario
        cp $scenario_pbs tides_${location}_${stage_file}.pbs
    done
done
