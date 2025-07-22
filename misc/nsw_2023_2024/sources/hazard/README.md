# Multiple importance sampling of random scenarios

Directories named like `scenarios_IDXXXXX` contain code to sample a batch of scenarios from PTHA18, using importance sampling.

The folder `multiple_importance_sampling_weights` contains code to compute the weights of each batch of scenarios, as used for multiple importance sampling. This involves computing the optimal solution for PTHA18 points, and interpolating/extrapolating these to get weights for the full model.
