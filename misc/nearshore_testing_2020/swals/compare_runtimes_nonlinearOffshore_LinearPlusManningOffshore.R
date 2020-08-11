source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')

# Get runtimes for linear-plus-manning model
all_wallclock_times_LM = lapply(Sys.glob('OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-linear_with_manning-0.035-highres_NSW/RUN_20200610_073112227/*.log'), get_domain_wallclock_times_in_log)
all_wallclock_times_LM = do.call(rbind, all_wallclock_times_LM)
all_wallclock_times_LM = aggregate(all_wallclock_times_LM$times, list(all_wallclock_times_LM$index), sum)
all_wallclock_times_LM
## Domains 1-4 are linear+Manning using a leapfrog solver
#     Group.1         x
#  1        1 20283.177
#  2        2 20163.383
#  3        3 20208.387
#  4        4 19902.936
#  5        5 31996.415
#  6        6 16082.635
#  7        7 16714.720
#  8        8  5745.297
#  9        9  9079.012
#  10      10  2215.998
#  11      11  3612.888
#  12      12  1512.726
#  13      13 36103.048

# Get runtimes for nonlinear model
all_wallclock_times_nl = lapply(Sys.glob('OUTPUTS/Tohoku2011_YamakaziEtAl2018-risetime_0-full-leapfrog_nonlinear-0.035-highres_NSW/RUN_20200615_123551400/*.log'), get_domain_wallclock_times_in_log)
all_wallclock_times_nl = do.call(rbind, all_wallclock_times_nl)
all_wallclock_times_nl = aggregate(all_wallclock_times_nl$times, list(all_wallclock_times_nl$index), sum)
all_wallclock_times_nl
## Domains 1-4 are full nonlinear using a leapfrog solver
#   Group.1          x
# 1        1 134361.901
# 2        2 134557.985
# 3        3 133401.738
# 4        4 134219.769
# 5        5  31900.313
# 6        6  16046.189
# 7        7  16516.690
# 8        8   5825.144
# 9        9   8944.345
# 10      10   2200.755
# 11      11   3510.673
# 12      12   1510.135
# 13      13  36045.603

# Here we see the ratio of ~ 0.15 for the first 4 domains -- that means a speedup of 1/0.15 = 6.666 due to the reduced-physics solver.
all_wallclock_times_LM$x / all_wallclock_times_nl$x
# [1] 0.1509593 0.1498490 0.1514852 0.1482862 1.0030126 1.0022713 1.0119897
# [8] 0.9862928 1.0150561 1.0069262 1.0291156 1.0017156 1.0015937
