# Maximum depth above initial condition for marine-warning and no-threat scenarios.

For each ATWS coastal zone, this folder contains:
- The maximum depth (above initial condition) attained by any marine warning scenario
- The maximum depth (above initial condition) attained by any no threat scenario

The results are limited to sites with elevation > 0. They are also clipped to
polygons defining the marine warning / no threat zones. Herein we use versions
of the latter polygons that were limited to sites where the 1/2500 84% tsunami
exceeded 1 cm above the initial sea level, to match what NSW SES have been
using.

Run the calculations with
```
source run_all.sh
```

Our initial JATWC-style calculations included the maximum waterlevels attained
by marine warning and no threat scenarios in each coastal zone. The
calculations here first converted them to the maximum depth (for marine warning
/ no threat respectively) by subtracting the elevation, and then convert them
to the "depth above initial condition", which is the minimum of the depth and
the "maximum waterlevel above the model's background sea level (1.1 m AHD)".


## Background on "depth above initial condition"

The "depth above initial condition" is only different to the depth at sites
where the elevation is below the model's background sea level (1.1 m AHD). But
at these sites it gives a better indication of the tsunami size. 

To see why, consider a site with an elevation of 0.5 m and a very small
tsunami, reaching just 0.01 m above the model's background sea level. The
maximum waterlevel is thus 1.11 m, while the "depth above initial condition" is
0.01 m (reflecting the small tsunami size), but the depth is 0.61 m. The latter
might be accidently interpreted as a significant tsunami. Thus we use "depth
above initial condition" to make it easier to identify sites that are flooded
at 1.1 m AHD and yet have a small tsunami. 

In sites where the model's background sea level is know to be too large, it is
reasonable to suppose that a model using the local background sea level would
have smaller modelled depths. For instance, suppose the 1.1m AHD sealevel used
by the model exceeds the local high tide (0.6 m AHD) by 50 cm. Say the modelled
depth above initial condition is only 30 cm for a marine warning (i.e.
substantially less than the overestimate of the tide). Then if we made a new
model which correctly represented the 0.6 m high tide (with no other changes),
it probably would not predict flooding during marine warnings at this site.
More precise results would require the modelling to account for spatial
variations in the tide.

