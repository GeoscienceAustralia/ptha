This folder contains several de-tided tide-gauge records from Nuku'alofa, as
well as code to create them from the original data and a  function for use in
our other plotting scripts, which converts the data to a consistent format.

*The original data is not included to keep file-sizes down -- instead we provide
the post-processed data - but the original data is available from the National
Tidal Unit of the Australian Bureau of Meterorology, tides@bom.gov.au.*

Key codes are:
* [./get_gauge_data_for_event.R](./get_gauge_data_for_event.R) provides a consistent interface to the data for our plotting codes in [../../swals/plots/](../../swals/plots).
* [./spectral_highpass_filter.R](./spectral_highpass_filter.R) contains a simple fft-based highpass filter function, which by default separates a time-series into a long-period component with wave periods longer than 3 hours, and a high-frequency component with shorter periods.

Key directories are:
* [./Tonga2006_BOM](./Tonga2006_BOM) includes a record of the 2006 Tonga tsunami
* [./BOM_2009_2010](./BOM_2009_2010) includes a record of the 2010 Chile Mw 8.8 tsunami
* [./Tohoku2011_IOC_Sealevel](./Tohoku2011_IOC_Sealevel) includes a record of the 2011 Tohoku tsunami
* [./BOM_2014_2020](./BOM_2014_2020) includes a record from late 2014 to mid 2020, which includes the 2015 Mw 8.3 Chile tsunami

