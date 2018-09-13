This folder contains various scripts + data that were used during PTHA18 for
working with DART buoy records.

[scrape_data.R](scrape_data.R) downloads the historical records from the web
- For the 2006/05/03 Tonga earthquake, I separately downloaded data at
  https://www.ngdc.noaa.gov/thredds/catalog/dart_bpr/rawdata/51407/catalog.html
  and processed it into a format similar to the other data before analysis. It
  was put into file "dart_extract/51407Bt2006.txt". This was done because the
  'regular' data source for that event had low temporal resolution (for some
  unknown reason).

[read_dart_metadata.R](read_dart_metadata.R) parses the buoy metadata.

[read_dart_files.R](read_dart_files.R) parses the files and has code for plotting their results.

[extract_tsunami_from_dart.R](extract_tsunami_from_dart.R) runs our manual de-tiding
of the data. For each event it includes choices of start times, smoothing parameters, etc,
which were derived by trial and error.
