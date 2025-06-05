# Nearshore points. PTHA18 can't perform here.
# Rscript extract_max_stage_at_a_point.R 150.789448242126 -23.1609999987272 # Rossyln Bay. 3 m deep. 3rd level nesting. TPXO9 2.56 m.
# Rscript extract_max_stage_at_a_point.R 151.33117821582 -23.8475383691669 # South Trees. 3 m deep. 3rd level nesting. TPXO9 2.15 m.

Rscript extract_max_stage_at_a_point.R 150.789448242126 -23.1609999987272 # Rossyln Bay. 3 m deep. 3rd level nesting. TPXO9 2.56 m.
Rscript extract_max_stage_at_a_point.R 151.33117821582 -23.8475383691669 # South Trees. 3 m deep. 3rd level nesting. TPXO9 2.15 m.
Rscript extract_max_stage_at_a_point.R 152.713302612305 -24.0557460784912  # 6 km N of Lady Elliot. 55 m deep. 2nd leve nesting. TPXO9 1.46 m.


# Matt, this point fails with an error (before GD made any code changes), see the logs.
# The problem is:
#  Error: Finding NA max-stage values. This can happen if target_point is right on
#  the boundary of two domains, given how the raster lookup works here. The
#  simplest workaround is to use a nearby point that isn't right on the boundary.
#  A more complex fix would be to make a vrt file with multiple neighbouring tifs,
#  and extract from that.
Rscript extract_max_stage_at_a_point.R 150.8625145 -23.5841507  # Port Alma Storm Surge. TPXO9 2.67 m.

# Offshore but near Gladstone/Yeppoon
Rscript extract_max_stage_at_a_point.R 151.180465698242 -23.1141948699951  # 20 km E of Keppel Islands. 28 m deep. 2nd level nesting. TPXO9 2.47 m.
Rscript extract_max_stage_at_a_point.R 151.516738891602 -23.8272743225098 # 15 km East of Facing Island (Gladstone). 25 m deep. 

# Deaggregation target point.
Rscript extract_max_stage_at_a_point.R 152.66667175293 -23.333333961162  # 60 km E of Heron Island. 360 m deep.  TPXO9 1.45 m.
# Somewhat close to the deaggregation point.
Rscript extract_max_stage_at_a_point.R 153.6666 -22.666666 # Depth about 350m, NE of the deaggrregation point 
Rscript extract_max_stage_at_a_point.R 155.0 -24.666666 # Depth about 4314m, SE of the deaggrregation point 

# These points are a fairly long way from Gladstone and don't work especially well.
# Rscript extract_max_stage_at_a_point.R 153.3333333 -18.0 # Depth around 2367 m, well offshore and further north. 
# Rscript extract_max_stage_at_a_point.R 153.5372222 -14.71472222 # DART buoy 55023, depth about 4629m
# Rscript extract_max_stage_at_a_point.R 151.0 -18.666666 # Depth about 1480m, but inside some shallow topography beyond the outer reef, well north.
