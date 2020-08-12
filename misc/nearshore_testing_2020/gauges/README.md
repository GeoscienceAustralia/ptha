Tide-gauge data
---------------

The script [gauge_data_links.R](gauge_data_links.R) is used in this study to provide a consistent interface to our tide-gauge data, and
 is repeatedly called in other scripts in [../swals](../swals) and [../analysis](../analysis). It makes sure all the data is converted to a consistent time-zone (GMT), and has consistent column names. It also stores the gauge coordinates, and gives us a convenient data-structure to look-up the data.

The script [gauge_data_links.R](gauge_data_links.R) relies on processed tide-gauge data that can be downloaded here:
http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/DATA.zip
and needs to be unzipped into the current directory (i.e. so there is a folder ./DATA inside this directory). In general, this processed data contains the time, the observed tidal-stage, and inferred tsunami after removal of tides and long-period components. The time-zone varies but is always converted to GMT by [gauge_data_links.R](gauge_data_links.R), which also contains the gauge coordinates. 

The processed tide-gauge data was derived from tide-gauge records provided by various organisations (see details in the README files in the ./DATA folder), in particular:
* The Western Australia Department of Transport: 
* The Bureau of Meteorology: 
* Manly Hydraulics Laboratory and the Department of Planning, Industry and Environment
* The NSW Port Authority

## Why are some records missing the raw tidal-stage values ? 
In the case of data obtained from the NSW Port Authority, we cannot provide the observed tidal-stage values because of commerical constraints on the data. 

We have agreement from the Port Authority to freely distribute the detided tsunami signal (extracted from the observed tidal-stage by detiding and removing long-period waves), but not the observed stage. Geoscience Australia is responsible for the detiding, and files derived from the Port Authority gauge data have NA in the stage column. 
