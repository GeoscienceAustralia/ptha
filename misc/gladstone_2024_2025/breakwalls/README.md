Breakwalls are terrain features that are enforced as high points in the elevation grid. These are to ensure no flow through thin terrain features that might not otherwise be well represented by the elevation in the model.
They are burned into the grid, just as for inverts in the [../swals/model_initial_conditions_mod.f90](../swals/model_initial_conditions_mod.f90).

We employ two techniques to generate them, either by specifying a [specific elevation](fixed_height) or by taking the [maximum elevation within a buffer zone](buffer_height) for each point on the line.
