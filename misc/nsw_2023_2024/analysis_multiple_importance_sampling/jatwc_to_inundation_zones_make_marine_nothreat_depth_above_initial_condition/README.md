# Maximum depth above initial condition for marine-warning and no-threat scenarios.

This folder contains:
- The maximum depth (above initial condition) attained by any marine warning scenario
- The maximum depth (above initial condition) attained by any no threat scenario

Our initial JATWC-style calculations included the maximum waterlevels attained
by marine warning and no threat scenarios. These can be converted to the
maximum depth (for marine warning / no threat respectively) by subtracting the
elevation. Similarly they can be converted to "depth above initial condition",
which is the minimum of the depth and the "maximum waterlevel above the model's
background sea level (1.1 m AHD)". 

The "depth above initial condition" is only different to the depth at sites
where the elevation is below the model's background sea level (1.1 m AHD). But
at these sites it gives a better indication of the tsunami size. To see why,
consider a site with an elevation of 0.5 m and a very small tsunami, reaching
just 0.01 m above the model's background sea level. The maximum waterlevel is
thus 1.11 m, while the "depth above initial condition" is 0.01 m (reflecting
the small tsunami size), but the depth is 0.61 m. The latter might be
accidently interpreted as a significant tsunami. Thus we use "depth above
initial condition" to make it easier to identify sites that are flooded at 1.1
m AHD and yet have a small tsunami. 

In sites where the model's background sea level is know to be too large, it is
reasonable to suppose that the actual tsunami depth will be less than the
modelled depths here (since real tsunamis will occur with a lower background
sea level). For instance, if the modelled 1.1m AHD sealevel exceeds the local
high tide (say 0.6 m AHD) by 50 cm, and the modelled depth above initial
condition is only 30 cm for a marine warning (i.e. substantially less than the
overestimate of the tide), then it seems likely that the site won't be flooded
during any marine warning. More precise results would require the modelling to
account for spatial variations in the tide.

