Tide-gauge data
---------------

The script [gauge_data_links.R](gauge_data_links.R) is used in this study to provide a consistent interface to our tide-gauge data, and
 is repeatedly called in other scripts in [../swals](../swals) and [../analysis](../analysis). It makes sure all the data is converted to a consistent time-zone (GMT), and has consistent column names. It also stores the gauge coordinates, and gives us a convenient data-structure to look-up the data.

The script [gauge_data_links.R](gauge_data_links.R) relies on processed tide-gauge data that can be [downloaded here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/DATA.zip). The resulting `DATA.zip` file needs to be unzipped into the current directory. At that point there should be a folder ./DATA inside this directory.

In general the processed data files contain the time, the observed tidal-stage, and the inferred tsunami after removal of tides and long-period components (the latter computed by Gareth Davies as discussed in the paper). The time-zone varies but is always converted to GMT by [gauge_data_links.R](gauge_data_links.R), which also contains the gauge coordinates. The tidal-stage values are in meters with respect to some unspecified vertical datum.

The processed tide-gauge data was derived from long-term tide-gauge records provided by various organisations, in particular:
* The Western Australia Department of Transport: https://catalogue.data.wa.gov.au/dataset/tide-stations
* The Bureau of Meteorology: http://www.bom.gov.au/oceanography/projects/ntc/ntc.shtml
* Manly Hydraulics Laboratory and the Department of Planning, Industry and Environment: https://mhl.nsw.gov.au/
* The NSW Port Authority: http://wavewindtide.portauthoritynsw.com.au/

In general the original data can be obtained from these institutions, and is available under the Creative Commons Attribution 4.0 licence (except for the commercial Port Authority Data). See the README files in ./DATA for additional information.

The processed data available at the download link can be used under the Creative Commons Attribution 4.0 Licence, consistent with the underlying data sources.

## Why are some records missing the raw tidal-stage values ? 
In the case of data obtained from the NSW Port Authority, we cannot freely distribute the observed tidal-stage values because of commerical constraints on the data. That data can be purchased from the NSW Port Authority via their web-page.

We have an agreement with the NSW Port Authority that allows us to freely distribute the detided tsunami signal (extracted from the observed tidal-stage by detiding and removing long-period waves), but not the observed tidal-stage. Geoscience Australia is responsible for the detiding, and the files we derived from the Port Authority gauge data have NA in the stage column. 
