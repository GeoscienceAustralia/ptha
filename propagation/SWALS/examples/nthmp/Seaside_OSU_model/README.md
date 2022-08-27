# Seaside OSU experiment of tsunami inundation through a city

Model of the Seaside OSU experiment from Park et al. (2013), which simulates tsunami inundation through a city.

The data for this problem was obtained from the NTHMP velocities benchmark suite. http://coastal.usc.edu/currents_workshop/problems.html. From that page there is also a link to presentations which discuss solutions from other modelling groups, and is good to look at to get an idea of the performance we might expect. In addition there is discussion of the results from different modelling groups in the workshop report (NTHMP, 2017).

## Issues

- The proper treatment of the landward boundary in this problem is unclear. Park et al. (2013) develop a numerical model of their experiment, and note that its landward boundary differs from the physical model, such that reflected waves arrive earlier in the model than in the experiment. We also see this in the model of Macias et al (2020) which is setup based on the description here (see gauge B9 in their paper). Herein we prevent wave reflection from the landward boundary by adding an artificial trough that catches the water, but suggest the correct treatment is unclear.
- Many modellers report under-prediction of depth maxima at gauge 'B1' (Macias et al. 2020; NTHMP 2017), and also 'A1 and C1' (Gao et al., 2020), relatively good results at gauges B4, B6, and varying results at B9 (e.g. over-estimation in Park et al; underestimation in Macais et al (2020)). In our plotting script I have added approximate "PASS/FAIL" critera based on the observed performance, and while useful for detecting future regressions in the model performance, the critera do not have strong justication.
- This benchmark problem can be forced using either a wave-maker forcing, or a provided model time-series. We do the latter, although it leads to some phase-offset at gauges WG3 and WG4, as noted previously by Macias et al (2020). The latter authors suggest this is probably due to a mistake in the provided initial condition. In future it would be good to try a wavemaker forcing.


## References

National Tsunami Hazard Mitigation Program. 2017. Proceedings and Results of the National Tsunami Hazard Mitigation Program 2015 Tsunami Current Modeling Workshop, February 9-10, 2015, Portland, Oregon: compiled by Patrick Lynett and Rick Wilson. 194 p., downloaded at https://nws.weather.gov/nthmp/documents/NTHMP_Currents_Workshop_Report.pdf

Gao, S.; Collecutt, G.; Syme, W. J. & Ryan, P. High resolution numerical modelling of tsunami inundation using quadtree method and GPU acceleration Proceedings of the 22nd IAHR-APD Congress 2020, Sapporo, Japan, 2020

Macías, J.; Castro, M. J. & Escalante, C. Performance assessment of the Tsunami-HySEA model for NTHMP tsunami currents benchmarking. Laboratory data Coastal Engineering, Elsevier BV, 2020, 158, 103667

Park, H.; Cox, D. T.; Lynett, P. J.; Wiebe, D. M. & Shin, S. Tsunami inundation modeling in constructed environments: A physical and numerical comparison of free-surface elevation, velocity, and momentum flux Coastal Engineering, Elsevier BV, 2013, 79, 9-21

