#
# In practice this code was run interactively in a terminal.
#

source('stats_under_null_hypothesis.R')
# Compute the main statistics used in the paper
stats_for_paper()

# Convert to tables in Latex format that are actually used in the paper.
convert_stats_to_latex_table(tsunami_size_type='tsunami_maxima', downsample=FALSE)
convert_stats_to_latex_table(tsunami_size_type='tsunami_maxima', downsample=TRUE)
convert_stats_to_latex_table(tsunami_size_type='stage_range', downsample=FALSE)
convert_stats_to_latex_table(tsunami_size_type='stage_range', downsample=TRUE)

