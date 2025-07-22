# Scripts for comparing specific source-zones with deep-ocean tsunami measurements

The scripts in this directory were used to compare source-zone results against
tsunami observations of DART buoys. It is included to illustrate how we did
this. However, please beware of the issues below.

**Note 1:** For the PTHA18 analysis, the scripts were not run in this
directory. To see how they were used, recall that for each source-zone we make
a copy of the [TEMPLATE](../TEMPLATE/) directory (as discussed
[here](../README.md) ). For example, to make the puysegur2 source-zone codes,
we rename ../TEMPLATE to ../puysegur2, and then run the codes inside that
directory. The check_dart_puysgur2.R code would have been placed into the
directory ../puysgur2/TSUNAMI_EVENTS/, alongside the other template codes, and
scenario files. It would have been run from that location after all of the
puysgur2 tsunami scenarios had been created.

**Note 2:** The scripts refer to data in locations starting with
"../../../../../DATA/TIDES/DART/dart_extract/". This was the location in our
local file system when running the PTHA18 analysis. A copy of this data has
been posted [here](https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/DATA/dart_extract.zip).
The codes used to process the DART data are in [this folder](./dart_process),
although at the time of the PTHA18 analysis they were located in
"../../../../../DATA/TIDES/DART/".


