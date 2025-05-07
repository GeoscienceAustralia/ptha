This folder contains initial conditions for historical tsunamis derived from published source inversions.

As the contents are large, in most cases the resulting source tifs have been removed from this repository, and can be downloaded separately here: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/sources.tar.bz2

This file should be extracted (e.g. `tar -jxf sources.tar.bz2`) and merged with the current folder.

## Notes
Most source models use the Kajiura filter code in the current directory. Exceptions are:
* The source models with a rise time used the Kajiura filter code in `./sources_with_rise_time`
    * The sumatra2005 model also used the latter code
* The Kermdec2021 source model used code provided by Fabrizio Romano (INGV) which already included a Kajiura filter, so no extra filtering was needed.
* The Chile 1960 source model directly inverts the water surface so does not need a Kajiura filter
