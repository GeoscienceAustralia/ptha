# Clean existing binary
rm ./dam_break
rm -r ./OUTPUTS
# Build the code
make -B -f make_dam_break > build_outfile.log

## Upstream 1.0, Downstream 0.1
# Run the job
./dam_break 1.0 0.1 > outfile.log
# Plot it and report tests
Rscript plot.R 1.0 0.1

## Upstream 1.0, Downstream 0.5
# Run the job
./dam_break 1.0 0.5 > outfile.log
# Plot it and report tests
Rscript plot.R 1.0 0.5

## Upstream 1.0, Downstream 0.01
# Run the job
./dam_break 1.0 0.01 > outfile.log
# Plot it and report tests
Rscript plot.R 1.0 0.01

## Upstream 1.0, Downstream 0.0001
# Run the job
./dam_break 1.0 0.0001 > outfile.log
# Plot it and report tests
Rscript plot.R 1.0 0.0001
