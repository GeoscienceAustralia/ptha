# Define the multidomain nesting geometry
----------------------------------------

If everything is correctly setup, then run with:
```
Rscript create_boxes.R
Rscript edit_boxes.R
```

The command `Rscript create_boxes.R` will run code to create the domain
bounding boxes. Within that script you need to define some variables that
describe the multidomain structure.

The current domains have been setup to have 3 levels of nesting, reflected in
the directory name `domains_1_3_15_15`. The outer level is 1-arcmin in each
direction, the second level is 3x finer, and the next level is 5x finer than
that. The fourth level is 5 times finer than its parent (the second level),
like the third level. Note that the cell size in each domain is defined in the
`../swals/multidomain_design` files, the corresponding refinement factors are 9
(from the global grid tp the first level), 6, 3. Which are defined relative to
the global grid of 1-arcmin cells.

The command `Rscript edit_boxes.R` will create variants of the files created
with the `create_boxes.R` script. The new files will contain information on
which domains can be coarsened. In some situations this is useful for speeding
up models (to prevent overly fine domains in deep areas).



## Example of the logic to setting up complicated nesting:
----------------------------------------------------------

It requires care to setup nesting while respecting the SWALS constraint that,
when 2 domains overlap, the coarser domain cells must be completely covered by
finer domain cells. 

Here is an example that shows reasoning (might not correspond to our actual model).

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
    * This implies each patch would enclose 162x162 parent-grid cells, or 54x54 coarsened-cells. Both are integers, so we are good --non-integers would cause problems.
    * Notice that here we could NOT have set the patch size to be (0.1/4)x(0.1/4) -- because it would not contain an integer number of parent grid cells.


