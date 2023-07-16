cd FujiSatake2007
Rscript Okada_vertical_component.R
cd ..
cd sumatra2005
Rscript Okada_exercise_solution.R
cd ..

# Make the kajiura versions
Rscript apply_kajiura_to_rasters.R 1 10 &
Rscript apply_kajiura_to_rasters.R 2 10 &
Rscript apply_kajiura_to_rasters.R 3 10 &
Rscript apply_kajiura_to_rasters.R 4 10 &
Rscript apply_kajiura_to_rasters.R 5 10 &
Rscript apply_kajiura_to_rasters.R 6 10 &
Rscript apply_kajiura_to_rasters.R 7 10 &
Rscript apply_kajiura_to_rasters.R 8 10 &
Rscript apply_kajiura_to_rasters.R 9 10 &
Rscript apply_kajiura_to_rasters.R 10 10 &
