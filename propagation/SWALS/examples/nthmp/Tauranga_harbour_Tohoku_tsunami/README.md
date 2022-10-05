# NTHMP tsunami currents benchmark problem 3: Tohoku tsunami at Tauranga Harbour, NZ

We simulate currents and stage in Tauranga harbour (New Zealand) following the Tohoku tsunami. The data for this problem is available at the [NTHMP velocities benchmark suite](http://coastal.usc.edu/currents_workshop/problems.html).

The [SWALS model](tauranga.f90) uses two nested domains (Figure 1). The outer domain has a coarse cell size (167 m), with factor 3 refinement in the inner domain. The model is forced along the top boundary with a stage time-series observed at the "A Beacon" gauge and a radiation treatment of the fluxes. Alternative boundary treatments are implemented in the code; while not discussed here, the results are moderately sensitive to the boundary treatment because the boundary so close to the coast, making in nontrivial to prescribe the stage while properly radiating outgoing waves. 

We compare the modelled and observed water elevations at four water-level gauges (A Beacon, Moturiki, Tug Berth, Sulphur Point), as well as modelled and observed flow speeds at the ADCP. The ADCP location is subject to some uncertainty and so two points are included there (discussed below).

![Figure 1: Model elevation, domain structure and gauge locations](Model_elevation_and_gauge_locations.png)


gauges_plot_lowresolution_omp.png
