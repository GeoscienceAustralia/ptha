Code template for the Australian Probabilistic Tsunami Hazard Assessment 2017
-----------------------------------------------------------------------------

[DATA](DATA) contains data that has been 'processed' or developed for modelling
  [e.g. sourcezone contours, edited DEMs, etc]. It should not contain RAW data,
  (the latter instead goes in '../../DATA/')

[SOURCE_ZONES](SOURCE_ZONES) contains unit source initial conditions and tsunami runs for
  every source-zone. This is the bulk of our modelling results.

[EVENT_RATES](EVENT_RATES) contains code to compute the rate of earthquakes on each source-zone,
  the tsunami wave height exceedance rates at each hazard point, etc.
