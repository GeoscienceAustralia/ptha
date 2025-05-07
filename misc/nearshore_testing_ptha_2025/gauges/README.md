# Uniform interface to a range of detided tide gauge data
---------------------------------------------------------

The script `gauge_data_links.R` provides a uniform interface to detided tide-gauge data that we have obtained from many sources over the years. It contains a list of tide gauges with:
* Location information
* A list of files containing tide gauge observations for different tsunami events. These were created by the author using datasets that were obtained as described in the `./DATA/` folder near each dataset (download as explained below).
* A function which reads this data and converts it to a uniform format. This ensures a uniform time zone and consistent column names, which the underlying data does not have.

To use the script you have to download and extract the associated data files here: https://thredds.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2025/DATA.tar.bz2
* They can be extracted with (e.g.) `tar -jxf DATA.tar.bz2`
* You should end up with a folder `./DATA/` containing all the files that `gauge_data_links.R` points to.

Note that code in other folders uses a different file path to refer to the `gauge_data_links.R` script (since our tide gauge data is stored in a central location). 
