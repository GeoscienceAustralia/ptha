# Create domains
----------------

The command
    `Rscript create_boxes.R`
will run code to create the domain bounding boxes, using various user-created inputs defined therein.

The command
    `Rscript edit_boxes.R`
will create variants of the files created with the `create_boxes.R` command. The 
new files will contain information on which domains can be coarsened, which are needed
for the WA model.



## Example of the logic to seting-up complicated nesting
--------------------------------------------------------

It requires care to setup nesting while respecting the SWALS constraint that,
when 2 domains overlap, the coarser domain cells must be completely covered by
finer domain cells. 

Here is an example that shows the reasoning (although might not correspond to our actual model).

* We have 1-arcmin outer grids with boundaries aligned to integer degrees or half degrees.
    * These have 60 cells/degree
* Next finer level is chosen to have 9x refinement and 0.5x0.5deg patches
    * This implies each patch would enclose 30x30 parent-grid cells. These numbers are integers, so it will work -- non-integers would cause problems. 
    * There would be 270x270 refined cells in these patches
* Next finer level is chosen to have 9x refinement (coarsened to only 3x in deeper cells) and (0.5/3)x(0.5/3) degree patches
    * This implies each patch would enclose 90x90 parent-grid cells. These numbers must be integers.
    * There would be 810x810 refined cells there (or 270x270 for coarsen-cells)
    * Notice that here we could NOT have set the patch size to be (0.1/4)x(0.1/4) -- because it would not contain an integer number of parent grid cells. 
* Next finer level is chosen to have 3x refinement in (0.1/3)x(0.1/3) degree patches
    * This implies each patch would enclose 162x162 parent-grid cells, or 54x54 coarsend-cells. Both are integers, so we are good --non-integers would cause problems.
    * Notice that here we could NOT have set the patch size to be (0.1/4)x(0.1/4) -- because it would not contain an integer number of parent grid cells.


