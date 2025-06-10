#!/bin/bash
# Helper script for the makefile to set compiler flags

# Get the flag type and target
FLAG_TYPE=$1
TARGET=$2

if [ -z "$FLAG_TYPE" ]; then
  echo "Usage: $0 <flag_type> <target>"
  exit 1
fi

if [ -z "$TARGET" ]; then
  echo "Usage: $0 <flag_type> <target>"
  exit 1
fi

if [ "$FLAG_TYPE" == "preprocessor" ]; then
    SWALS_PREPROCESSOR_FLAGS=" -DTIMER -DSPHERICAL -DCOARRAY -DCOARRAY_PROVIDE_CO_ROUTINES -DCOARRAY_USE_MPI_FOR_INTENSIVE_COMMS -DLOCAL_TIMESTEP_PARTITIONED_DOMAINS"
    if [[ "$TARGET" == "debug" ]]; then
        SWALS_PREPROCESSOR_FLAGS+=" -DTRACK_MULTIDOMAIN_STABILITY -DEVOLVE_TIMER"
    elif [[ "$TARGET" == "old_nesting" ]]; then
        SWALS_PREPROCESSOR_FLAGS+=" -DOLD_PROCESS_DATA_TO_SEND_B4FEB22"
    fi
    echo $SWALS_PREPROCESSOR_FLAGS
elif [ "$FLAG_TYPE" == "arch" ]; then
    SWALS_FC_ARCH_FLAGS=" -xSAPPHIRERAPIDS"
    if [[ "$TARGET" == "skylake" ]]; then
        SWALS_FC_ARCH_FLAGS=" -xSAPPHIRERAPIDS -qopt-zmm-usage=high"
    elif [[ "$TARGET" == "cascadelake" ]]; then
        SWALS_FC_ARCH_FLAGS=" -xCASCADELAKE"
    elif [[ "$TARGET" == "broadwell" ]]; then
        SWALS_FC_ARCH_FLAGS=" -march=broadwell -axSKYLAKE,CASCADELAKE,SAPPHIRERAPIDS -qopt-zmm-usage=high"
    fi
    echo $SWALS_FC_ARCH_FLAGS
fi
