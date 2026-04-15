# Kalbarri to Coral Bay hazard calculations

This model includes a areas from south of Kalbarri to Coral Bay, with high resolution in most towns (with exceptions where preliminary runs showed the tsunami to be small).

See [doc](doc) for general guidance on using the code, along with further documentation and code in the subfolders.

Key subfolders are
* [breakwalls](breakwalls) - enforce elevation high points in the model
* [elevation](elevation) - specify location of elevation files used by the model
* [gauges](gauges) - locations where the model stores time-series.
* [initial_stage](initial_stage) - make the initial stage lower than the default value, typically to eliminate ponds on land.
* [inverts](inverts) - enforce elevation low points in the model grid
* [multidomain_design_B](multidomain_design_B) - Nesting geometry for **the model used for hazard calculations**.
  * [multidomain_design](multidomain_design) has the nesting geometry for **a preliminary model that was not ultimately used**. It covers more areas at high resolution, including north of Coral Bay, but was too expensive for probabilistic hazard calculations.
* [sources](sources) - Tsunami model initial conditions for testing and hazard calculations.
* [swals](swals) - Tsunami model run files.
* [tides](tides) - Basic plotting of tides using TPXO9
