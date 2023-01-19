[combine_tsunami_sources.R](combine_tsunami_sources.R) -- Code to combine
tsunami unit-source initial conditions into event initial conditions. It also 
shows how to apply Kajiura filtering after combining the unit sources (although
activating that code requires setting kajiura=TRUE and pointing the code to
appropriate elevation data).

__NOTE__: In real applications you might not wish to use all random scenarios.
For example, in PTHA18 we rejected scenarios with overly high peak slip, and
imposed other constraints on the scenario magnitude (dependent on the source
geometry). See Section 3.2.3 of the PTHA18 report.
