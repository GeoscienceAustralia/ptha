The following programs might be useful, although there can be many ways to do things.

* R (or Rscript for non-interactive usage)
  * Many components within the language
  * Useful scripts in subdirectories, but might need modification.
  * Pretty good for GIS processing, scanning log files, plotting, etc
  * [plot.R](../../../propagation/SWALS/plot.R) (in the SWALS github head directory) has useful routines for working with SWALS outputs, creating load balance files, scanning logs, etc.
* Some kind of GIS for interactive work (I use QGIS)
* `scp` to move files to/from NCI
* `git` for version control:
  * `git init` to allow git to begin.
  * `git add file1.txt file2.txt` to add two files to the project
  * `git commit -a -m "my commit message here"` to make a commit, which will snapshot the state of the files that have been added.
  * `git ls-files` to see which files are under version control
  * `git log` to see past commit messages
  * `git diff` to see how the tracked files differ from the last commit.
  * For working with online repositories
    * `git remote-add`
    * `git pull`
* GDAL for many things
  * Make a VRT file (collection of multiple rasters). 
    * `gdalbuildvrt -resolution highest myvrt.vrt inputfiles*.tif`
    * Any raster overlaps should be missing data for all but one raster.
  * `gdal_translate` to efficiently convert file formats 
  * `gdalwarp` for reprojection of rasters
  * Don't forget to add the option `-co COMPRESS=DEFLATE` when creating tifs with the above programs.
* General terminal commands `cd, cp, mv, ls, ls -halt, top` and their combination with `|`
* `grep` to scan files for a pattern
  * `grep pattern file.txt` to print lines matching `pattern` in `file.txt`
  * `grep -C4 pattern file.txt` for 4 lines around match (vs `-B` before and `-A` after).
  * `for i in *.log; do echo $i; grep -A2 pattern $i; done` to apply grep to search for `pattern` in all the log files.
* Quick viewing of netcdf files
    * `ncview name_of_netcdf_file.nc`
* Print a summary of information in netcdf files
    * `ncdump -h name_of_netcdf_file.nc`
* Get a file hash with `md5sum myfile`. 
  * This is useful if you've got 2 files in different locations and want to see if they differ.
  * For example: to check that you're using the most up-to-date files from another machine.
