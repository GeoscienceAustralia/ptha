# Avoid errors when bundling as a library (not that we use the resulting library)
SWALS_AR :=/apps/intel-tools/intel-compiler-llvm/2025.2.0/bin/compiler/llvm-ar rcs

# Fortran compiler
SWALS_FC := mpif90
SWALS_FC_FLAGS := -fpp -O3 -qopenmp -ipo -g
#SWALS_FC_FLAGS := -cpp -Ofast -march=native -fopenmp -g -ffree-line-length-1024
#SWALS_FC_FLAGS := -cpp -O3 -march=native -fopenmp -g -ffree-line-length-1024

# Call a bash script to set the flags
MAKECMDGOALS_NO_TARGET := $(subst _target,,$(MAKECMDGOALS))
SWALS_PREPROCESSOR_FLAGS := $(shell ./set_flags.sh preprocessor $(MAKECMDGOALS_NO_TARGET))
SWALS_FC_ARCH_FLAGS := $(shell ./set_flags.sh arch sapphirerapids $(MAKECMDGOALS_NO_TARGET))

# Set the name for the main model file $(mymodel).f90
mymodel := model

# Include variables we typically need for compilation
SWALS_SRC := /g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/src
include $(SWALS_SRC)/src_standard_compiler_var

# wrapper to build and log the build 
build debug old_nesting both_debug_and_old_nesting cascadelake broadwell skylake:
	$(MAKE) $@_target > $(mymodel)_$@.log 2>&1

# Build the model of choice -- note the previous command appends _target to whatever
# the build type is.
build_target debug_target old_nesting_target both_debug_and_old_nesting_target cascadelake_target broadwell_target skylake_target: clean $(SWALS_LIBRARY) $(mymodel).o
	$(SWALS_FORTRAN) $(mymodel).o -o $(mymodel)_$(subst _target,,$@) $(SWALS_OBJECTS) $(SWALS_FC_LIBS)
	$(MAKE) clean

# Compile the driver script
$(mymodel).o: $(mymodel).f90
	$(SWALS_FORTRAN) -c $^

# Build the SWALS library
include $(SWALS_SRC)/src_make_commands

clean:
	rm -f *.o *.mod

clean_all: clean
	rm -f model_build model_debug model_old_nesting model_both_debug_and_old_nesting model_cascadelake model_broadwell model_skylake
