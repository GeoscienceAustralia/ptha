[combine_tsunami_gauges.R](combine_tsunami_gauges.R) -- Code to combine tsunami
unit-source tide-gauges into events. Currently has hard-coded paths to input
data.

[combine_tsunami_sources.R](combine_tsunami_sources.R) -- Code to combine
tsunami unit-source initial conditions into event initial conditions. Kajiura
filter is applied after combining the unit sources, since numerically this
approach is less prone to artefacts caused by interpolation/discretization.
