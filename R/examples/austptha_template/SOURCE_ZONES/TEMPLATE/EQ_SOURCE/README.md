Create a discretized source zone and compute unit-source tsunami initial conditions
-----------------------------------------------------------------------------------

Before running this, ensure you have correctly defined inputs in this csv file:
[../../../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv](../../../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv),

See the file
[../../../DATA/SOURCEZONE_PARAMETERS/README.md](../../../DATA/SOURCEZONE_PARAMETERS/README.md)
for more information.

In particular, the codes here rely on the `sourcename`, `rake`, and the
`approx_unit_source_length` and `approx_unit_source_width` being defined, and
not subsequently changed. The rake should either be 90 or -90 if you are
working with codes in other parts of these template scripts as well.
Furthermore, you should have created the SOURCEZONE_CONTOURS shapefile, and the
SOURCEZONE_DOWNDIP_LINES shapefile.

Beyond that, it is not necessary that other parameters in
sourcezone_parameters.csv are assigned 'final' values at this stage of the
analysis.

You also have to correctly set some parameters in [config.R](config.R). This is used in the main
run script, [produce_unit_sources.R](produce_unit_sources.R). 

On raijin.nci.org.au, the code can be run with
[run_produce_unit_sources.PBS](run_produce_unit_sources.PBS):

    qsub run_produce_unit_sources.PBS
