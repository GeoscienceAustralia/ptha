# Kalbarri to Coral Bay hazard calculations

Key folders are
* [initial_stage](initial_stage) - Used to make the initial stage lower than the default value, typically to eliminate ponds on land.
* [inverts](inverts) - Used to enforce elevation low points in the model grid
* [multidomain_design_B](multidomain_design_B) - Nesting geometry for **the model used for hazard calculations**.
  * [multidomain_design](multidomain_design) has the nesting geometry for **a preliminary model that was not ultimately used**. It covers more areas at high resolution, including north of Coral Bay, but was too expensive for probabilistic hazard calculations.
* [sources](sources) - Tsunami model initial conditions for testing and hazard calculations.
* [swals](swals) - Tsunami model run files.
* [tides](tides) - Basic plotting of tides using TPXO9
