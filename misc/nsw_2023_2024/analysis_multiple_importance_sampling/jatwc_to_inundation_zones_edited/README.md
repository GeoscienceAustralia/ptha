# Edited JATWC-Style land-warning polygons

The JATWC-style inundation polygons made for the NSW study included all areas
within the 1/2500 84% inundation zone. Some of these areas have negligably
small modelled tsunami at this exceedance rate. The report recommends manual
editing of the evacutation zones, and part of this process will involve
removing areas with an insignificantly small tsunami size (while also extending
zones elsewhere for ease of communication, to treat model shortcomings due to
the 1/54 arcminute grid size, to give site specific consideration of risk
tradeoffs for sites such as hospitals and nursing homes, etc).

Here we edited the land-warning zones to remove sites where the 1/2500 84%
max-stage does not exceed 1 cm above the 1.1m background sea level  (i.e.
max-stage less than 1.11 m). This was done with the included R script.

The approach used here is not aiming to replace manual editing of the
evacuation zones. Rather, the aim is to provide a quick-and-dirty adjustement
that removes many areas with low hazard which would likely be addressed with
manual edits, and thereby somewhat reduce the impression that extensive "more
inland" areas are subject to significant hazard. 
