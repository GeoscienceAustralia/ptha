# Make animations of the model's stage

The main scripts are [make_frames.R](make_frames.R) then [merge_frames.R](merge_frames.R) to combine them into a video.

First set up the inputs, which are sourced from within these scripts:
- [input/scenes.in.R](input/scenes.in.R) defines information for each location, e.g. the bounding box coordinates,
- [input/pusegur_84.in.R](input/pusegur_84.in.R) containing information about the scenario (like file paths) and the sequence of shots (how long each scene lasts).

Once set-up, run `make_frames.R` once with an internet connection to download imagery. Pass one shot's directory name as the first argument:
```{bash}
Rscript make_frames.R se_australia
```
This will then halt if it detects it's on a login node on Gadi.

Making the frames is embarrisingly parallel, so get an interactive node (or a computer with multiple cores) and run it again `Rscript make_frames.R se_australia`, making sure to set the mc.cores variable appropriately.

Setting the aspect ratio in the scenes data frame can be a little iterative. The plotting aspect ratio is around 0.675, whereas the requested scene's (i.e. downloaded imagery's) lower left and upper right might produce a different ratio. This can be resolved by adding a buffer `dbox` or iteratively adjusting the shot's bounding lower left and upper right. Although imagery isn't downloaded correctly if it crosses the international date line.

Once the frames are made, `module load ffmpeg` and combine the frames into a mp4 video file with:
```{bash}
Rscript merge_frames.R
```
