# Make some pretty animations of the model

## Inputs
Inputs are stored in `.R` files which are sourced as needed. You'll need two input files:
1. a [input/scenes.in.R](input/scenes.in.R) to define the scenes (lat-lon extents) you want. If only a subset of all the domains cover your area, it's efficient to specify which ones in the `scenes$domains` (by changing `NULL` to a vector).
2. a scenario specific file e.g. [input/kermadec_94.in.R](input/kermadec_94.in.R). This sets up file paths to the SWALS simulation (which must temporal stage saved) and extracts the times it saved. It also defines when to cut to each scene using a `shots` tibble/dataframe.

The colourmaps are also defined here in the [input/custom_colours.R](inputs/custom_colours.R) file.

## The four main scripts
1. [make_frames.R](make_frames.R) is used to generate the png frames for one scene e.g.
```
Rscript make_frames.R sw_pacific
```
where the first argument `sw_pacific` is the `dir` (directory) of the scene you want.
  - Run once with an internet connection (i.e. on an NCI login node) to download the background imagery. If running on NCI this will automatically stop if after downloading the imagery it detects you're on a login node (by checking the hostname).
  - Then run on a compute node. Making frames is "embarrassingly parallel" so benefits from having more cores (up to one per frame).
  
2. [run_frames.pbs](run_frames.pbs) is a PBS script that automates running [make_frames.R](make_frames.R) for many scenes. You could use it as a bash script to download all the imagery on a login node.
  ```
  . run_frames.pbs
  ```
And then actually make the frames by submitting it as a PBS job
  ```
  qsub run_frames.R
  ```

3. [make_intro.R](make_intro.R) to make intro slides.
4. [merge_frames.R](merge_frames.R) to put all the frames into a single directory with frame numbers and combine them into a video with ffmpeg.

## R Subfunctions
Subfunctions are all put inside the [R/](R/) directory. These files define functions used in the main scripts: `osm_backdrop.R` to download and cache the OpenStreetMap, `make_image_2d.R` to generate the png image, `combine_frame.R` to combine elements (side panels and headers) in each frame and `make_ic.R` to plot scenario initial conditions. 

## Choose scenario
The [choose_scenario/plot_scenario.R](choose_scenario/plot_scenario.R) script plots the max stage for tarred runs, which is useful for selecting a scenario before rerunning with temporal saving (which can take significant disk space).
