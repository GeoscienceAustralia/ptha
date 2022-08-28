source ./shell_var.sh

qcc -Ofast -fopenmp -Wall -o simple_wave simple_wave.c -lm
./simple_wave > db.txt #> out.ppm 2>log
Rscript plot.R
