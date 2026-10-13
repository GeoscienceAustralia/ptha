# Single grid tsunami model with animation

This code was used to make an animation of a single-grid tsunami scenario at the Puysegur trench. 

It uses the single grid SWALS code here [../../propagation/SWALS/examples/generic_example/](../../propagation/SWALS/examples/generic_example/), and thus the animation code we often use needs modification, because no multidomain is used in the SWALS model so outputs are somewhat different. 

To save space the tsunami source is not included. But can be extracted from PTHA18 (`puysegur2_row_0003046_Mw_84_HS.tif` -- i.e. row 3046 of the `puysegur2` HS scenarios) [in the usual way](../../ptha_access).
