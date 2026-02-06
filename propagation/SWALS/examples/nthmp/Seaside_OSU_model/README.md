# NTHMP tsunami currents benchmark problem 4: Seaside OSU experiment of tsunami inundation through a city

We model the Seaside OSU experiment from Park et al. (2013), which simulates tsunami inundation through a city. A movie of this problem is [available here](https://www.youtube.com/watch?v=nj98sHcTGOo).

Data for this problem was obtained from the [NTHMP velocities benchmark suite](http://coastal.usc.edu/currents_workshop/problems.html). That page also links to presentations from other modelling groups, and can help give an idea of the performance we might expect. See also the NTHMP currents workshop report (NTHMP, 2017).

Figures 1 and 2 below show the [SWALS model](model.f90) structure, elevation data and locations of gauges. The model uses two nested domains with 4x grid refinement near the "city". A wave time-series is applied from the west. Gauges WG3, WG4 are used to check that the boundary forcing is reasonable. The city includes 4 rows of gauges (A1-9, B1-9, C1-9, D1-4), where the model is compared with data on depths, cross-shore velocities, and cross-shore momentum fluxes $(u^2h)$. Here "cross-shore" is taken to be the x-direction as specified in the NTHMP problem description.

![Figure 1: Multidomain structure (red lines), elevation and gauge locations. The wave arrives from the west.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/Model_elevation_and_gauges.png)

![Figure 2: Inner domain elevation and gauge locations.](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/Model_elevation_and_gauges_zoom.png)

The model is forced using a provided synthetic stage time-series at x=5m, which bypasses the need for a wavemaker forcing around x=0m. The boundary speeds are specified with the plane wave relation $(uh = \sqrt{g h} stage)$. Both are provided to a Flather boundary condition which aims to enforce the wave while minimising spurious boundary reflections. Figure 3 compares the modelled and prescribed boundary stage, along with a gauge observation at x=2m (outside our domain, but inside the experimental domain). The model does a reasonable job of representing the boundary forcing, and increases the stage at later times (consistent with the observed gauge) to reduce boundary reflections.

![Figure 3: Boundary time-series](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/boundary_check.png)

Figure 4 compares the model and observations at gauges WG3 and WG4. They are in general agreement once the model is time-shifted to the left by 0.7s, which previous studies also found to be needed (Macais et al., 2020) and attributed to an error in the forcing. We do not expect perfect agreement because the flow is physically affected by non-hydrostatic physics, which is not represented with the nonlinear shallow water equations. This can be seen most clearly in the latter part of Figure 4, where the shallow water model predicts a shock wave that arrives earlier than the observed dispersive (short period) waves.

![Figure 4: Time-series at gauges WG3 and WG4](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/gauges_wg3_wg4.png)

# Model performance at the urban gauges

Below we compare the modeled and observed depths, cross-shore speeds, and cross-shore momentum fluxes. The model is in general agreement with the data. There are a few cases where the model does less well (e.g. the momentum fluxes for gauges B1 and C1), similar to reports in other studies. 

![Figure 5: Modelled and observed depths at gauges A1-A9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_A_depth.png)
![Figure 6: Modelled and observed cross-shore speed at gauges A1-A9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_A_speed.png)
![Figure 7: Modelled and observed cross-shore momentum flux at gauges A1-A9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_A_convective_flux_hv2.png)

![Figure 8: Modelled and observed depths at gauges B1-B9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_B_depth.png)
![Figure 9: Modelled and observed cross-shore speed at gauges B1-B9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_B_speed.png)
![Figure 10: Modelled and observed cross-shore momentum flux at gauges B1-B9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_B_convective_flux_hv2.png)

![Figure 11: Modelled and observed depths at gauges C1-C9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_C_depth.png)
![Figure 12: Modelled and observed cross-shore speed at gauges C1-C9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_C_speed.png)
![Figure 13: Modelled and observed cross-shore momentum flux at gauges C1-C9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_C_convective_flux_hv2.png)

![Figure 14: Modelled and observed depths at gauges D1-D9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_D_depth.png)
![Figure 15: Modelled and observed cross-shore speed at gauges D1-D9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_D_speed.png)
![Figure 16: Modelled and observed cross-shore momentum flux at gauges D1-D9](https://github.com/GeoscienceAustralia/ptha/blob/figures/propagation/SWALS/examples/nthmp/Seaside_OSU_model/urban_gauge_group_D_convective_flux_hv2.png)

## Issues

- The proper treatment of the landward boundary in this problem is unclear. Park et al. (2013) develop a numerical model of their experiment, and note that its landward boundary differs from the physical model, such that reflected waves arrive earlier in the model than in the experiment. We also see this in the model of Macias et al (2020) (see gauge B9 in their paper). Herein we prevent wave reflection from the landward boundary by adding an artificial trough that catches the water, but suggest the correct treatment is unclear.
- Many modellers report under-prediction of depth maxima at gauge 'B1' (Macias et al. 2020; NTHMP 2017), and also 'A1 and C1' (Gao et al., 2020), relatively good results at gauges B4, B6, and varying results at B9 (e.g. over-estimation in Park et al; underestimation in Macais et al (2020)). In our plotting script I have added approximate "PASS/FAIL" critera based on the observed performance, and while useful for detecting future regressions in the model performance, the critera do not have strong justication.
- This benchmark problem can be forced using either a wave-maker forcing, or a provided model time-series. We do the latter, although it leads to some phase-offset at gauges WG3 and WG4, as noted previously by Macias et al (2020). The latter authors suggest this is probably due to a mistake in the provided initial condition. In future it would be good to try a wavemaker forcing.
- The velocity and flux data are in the "cross-shore" direction. The problem indicates this is the x direction, which we use here, although it might instead be the along-street direction which is slightly different?

## References

National Tsunami Hazard Mitigation Program. 2017. Proceedings and Results of the National Tsunami Hazard Mitigation Program 2015 Tsunami Current Modeling Workshop, February 9-10, 2015, Portland, Oregon: compiled by Patrick Lynett and Rick Wilson. 194 p., downloaded at https://nws.weather.gov/nthmp/documents/NTHMP_Currents_Workshop_Report.pdf

Gao, S.; Collecutt, G.; Syme, W. J. & Ryan, P. High resolution numerical modelling of tsunami inundation using quadtree method and GPU acceleration Proceedings of the 22nd IAHR-APD Congress 2020, Sapporo, Japan, 2020

Macías, J.; Castro, M. J. & Escalante, C. Performance assessment of the Tsunami-HySEA model for NTHMP tsunami currents benchmarking. Laboratory data Coastal Engineering, Elsevier BV, 2020, 158, 103667

Park, H.; Cox, D. T.; Lynett, P. J.; Wiebe, D. M. & Shin, S. Tsunami inundation modeling in constructed environments: A physical and numerical comparison of free-surface elevation, velocity, and momentum flux Coastal Engineering, Elsevier BV, 2013, 79, 9-21

