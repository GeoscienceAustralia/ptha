Here we run the conical island problem using FUNWAVE-TVD (nondispersive) to
compare the maximum runup with SWALS.

This is based on the FUNWAVE benchmark problem cited below, with the following changes for
a clearer comparison with SWALS:
- No dispersion
- Grid refined to 0.025m (vs original 0.05m -- I interpolated their input data)
- Wet-dry thresholds reduced to 1.0e-05 (previously 1.0e-03) -- except for caseC, where I used 1.0e-04 to avoid some numerical artefacts that caused substantially non-symmetric runup around the top/bottom island.
- Manning friction


See the original model setup and data here (accessed 10/11/2019) -- the full source-code is also in this repository:
https://github.com/fengyanshi/FUNWAVE-TVD/tree/master/benchmarks/car_conical_island

To run the code, you will have to:
1) unzip 'input.zip' so there is a folder 'input' in this directory
2) Compile the funwave executable (compiled to use manning friction, see their manual for details)
3) Change the 'run_model.sh' scripts in caseA/caseB/caseC to point to your funwave executable
4) Execute the run_model.sh scripts in caseA/caseB/caseC
5) After the models have finished, run 
    Rscript check_results.R
   from inside the caseA/caseB/caseC directories, which will make the max_island_runup.csv files.
    
