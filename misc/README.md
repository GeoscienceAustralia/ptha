# Miscellaneous project code

This folder contains code used for various tsunami projects, which are described in our papers and technical reports. 

The code structure usually reflects the computational environment in which our calculations were implemented. For instance, many scripts include hard-coded file paths and assume that particular software is already installed. We are not aiming to make the calculations trivial to rerun on other machines. 

* [SW_WA_2021_2024](SW_WA_2021_2024) - Code for some models in south-west Western Australia. Parts of this were used in [this ICCE 2022 conference paper](https://doi.org/10.9753/icce.v37.papers.18) and [this paper in the Australian Journal of Emergency Management](https://knowledge.aidr.org.au/resources/ajem-october-2024-science-informed-risk-reduction-for-earthquake-generated-tsunamis-in-western-australia/). The work is comprehensively described in [this technical report](https://dx.doi.org/10.26186/150015).

* [WA_2025_2027](WA_2025_2027) - Code for some models in the region around Kalbarri - Onslow in Western Australia.

* [./hunga_tonga_data_paper](./hunga_tonga_data_paper) - Code for a paper presenting [observations of the Hunga-Tonga tsunami and pressure wave in Australia](https://doi.org/10.1038/s41597-024-02949-2).

* [./monte_carlo_paper_2021](./monte_carlo_paper_2021) - Code and links to data for a [paper on Monte-Carlo techniques for PTHA](https://doi.org/10.1093/gji/ggac140). While the Monte-Carlo techniques used in this study are general, they are illustrated via a probabilistic inundation hazard analysis for Tongatapu. 

* [./nearshore_testing_2020](./nearshore_testing_2020) - Code and links to data for [a paper](https://www.frontiersin.org/articles/10.3389/feart.2020.598235/full) which simulated tsunamis using various finite-fault inversions, and compares with nearshore tsunami data in Australia.

* [./nearshore_testing_ptha_2025](./nearshore_testing_ptha_2025) - Code and links to data for [a paper](https://doi.org/10.1029/2025JB031949) which compares tsunamis observed at Australian tide gauges with stochastic tsunamis from the 2018 Australian Probabilistic Tsunami Hazard Assessment.

* [./nsw_2023_2024](./nsw_2023_2024) - Code for the NSW broadscale model, described [here](https://www.researchgate.net/publication/396141903_NSW-wide_probabilistic_tsunami_inundation_hazards_from_subduction_earthquakes) and [here](http://dx.doi.org/10.26186/150281)

* [./probabilistic_inundation_tonga2020](./probabilistic_inundation_tonga2020) - A probabilistic inundation hazard analysis for Tongatapu. I suggest to use an improved variant of this analysis which is available in [./monte_carlo_paper_2021](./monte_carlo_paper_2021)

* [./gladstone_2024_2025](./gladstone_2024_2025) - Code for the Gladstone, Queensland model, described [here](https://doi.org/10.26186/150488) and [here](https://www.researchgate.net/publication/396142067_Accommodating_Spatially_Varying_Tidal_Planes_in_Tsunami_Hazard_Assessments). The SWALS model showcases variable friction, makefiles to compile inputs, elevation adjustments for an uneven but quiescent initial water level, a domains nested with a tree structure, breakwalls, inverts, setting initial stages dry for inland areas below sea level, and running at various sea levels.
