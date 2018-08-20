PB2002_steps_dat.zip
--------------------

    This data is part of the supplementary information for the paper:

    Bird, P. (2003) An updated digital model of plate boundaries. Geochemistry, Geophysics, Geosystems 4(3), doi:10.1029/2001GC000252

    It was downloaded from Peter Bird's website in 06/2017. See the ftp link on:
    http://peterbird.name/publications/2003_PB2002/2003_PB2002.htm



reformat_edited_Bird_convergence.R
----------------------------------

    Script to process subduction zone convergence information in a consistent way.
    
    This script uses an edited version of Bird's geometries (where many sources were deleted, but no other changes were made), to
    produce output in a format that is consistent with a post-processed version of Jonathan Griffin's source-zone traces. It
    also adds info where missing (for Arakan + Seram-south). See the script for further information & justifications.

combine_traces.R
-----------------
    This merges the intermediate datasets made above.   
 
    Note not all intermediate datasets are included in the repository (to
    reduce file size). However, sourcezone_traces_table_merged_uncorrected.csv is the
    'final' output of that script

fix_outer_rise_rates.R
-----------------------
    This script updates the 'sourcezone_traces_table_merged_uncorrected.csv' file so that
    the outer-rise rates are 'effective divergence rates' (see the PTHA18 report for more info).
    
    It produces the final 'sourcezone_traces_table_merged.csv'.

