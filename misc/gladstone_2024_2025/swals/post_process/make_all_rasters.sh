#!/bin/bash
# source ../../modules_R_431.sh


# declare all variables to make
# declare -a vars=("elevation0" "max_stage" "max_speed")  # test case
# declare -a vars=("max_speed" "last_timestep_UH" "last_timestep_VH")  # small source
# declare -a vars=("max_stage" "min_stage" "max_speed" "last_timestep_UH" "last_timestep_VH")  # extreme source
# declare -a vars=("max_stage" "max_depth" "min_stage" "max_speed" "arrival_time" "max_flux" "time_of_max_stage")
declare -a vars=("max_stage" "max_depth" "min_stage" "max_speed" "arrival_time" "time_of_max_stage" "elevation0")  # validation

# iterate through all directories in ../OUTPUTS/tides
for dir in ../OUTPUTS/solomon2007_wei*; do
    # get the last RUN_* directory
    dir=$(ls -d $dir/RUN_* | tail -n 1)
    echo $dir

    # rm $dir/*.tif
    Rscript make_rasters.R $dir 48 ${vars[@]}
    
    # for each variable, make a vrt
    for var in ${vars[@]}; do
        if [ $var == "max_depth" ]; then
            gdalbuildvrt -resolution highest -overwrite $dir/depth.vrt $dir/depth*.tif
        elif [ $var == "last_timestep_UH" ]; then
            gdalbuildvrt -resolution highest -overwrite $dir/UH_last_timestep.vrt $dir/UH*.tif
        elif [ $var == "last_timestep_VH" ]; then
            gdalbuildvrt -resolution highest -overwrite $dir/VH_last_timestep.vrt $dir/VH*.tif
        else
            gdalbuildvrt -resolution highest -overwrite $dir/$var.vrt $dir/$var*.tif
        fi
    done
done

