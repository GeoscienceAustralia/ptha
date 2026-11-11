# Make animations of the model running.

The three main scripts are [make_frames.R](make_frames.R), [make_intro.R](make_intro.R) to generate the png frame files then [merge_frames.R](merge_frames.R) to combine them into a video. Inputs such as [input/scenes.in.R](input/scenes.in.R), and [input/sumatra_2004.in.R](input/sumatra_2004.in.R) are sourced from within these scripts. They set up file paths to the SWALS simulation (with temporal stage saved) and scenes to make for each scenario. Use [run_frames.pbs](run_frames.pbs) to run the make_frames.R for all shots. It needs to be run once with an internet connection to download imagery and will then halt if on a login node on Gadi. If having issues, try running each `Rscript` line of `run_frames.pbs` on a login node. Then submit as `qsub run_frames.pbs` to make all the frames in parallel.

Setting the aspect ratio in the scenes data frame can be a little iterative. The plotting aspect ratio is around 0.675, whereas the requested scene's (i.e. downloaded imagery's) lower left and upper right might produce a different ratio. This can be resolved by adding a buffer `dbox` or iteratively adjusting the shot's bounding lower left and upper right. 

Subfunctions are all put inside the [R/](R/) directory.

## Useful files for previous iterations:
- [for Gladstone animations](../../../gladstone_2024_2025/animation/)
