Generation of tsunami events
-----------------------------

To make tsunami events, run the code in the following folders in order (full
explanation is provided in the README in each folder).

# Step 1

[EQ_SOURCE](EQ_SOURCE) contains code to make tsunami unit source initial
conditions

# Step 2

[TSUNAMI_UNIT_SOURCE](TSUNAMI_UNIT_SOURCE) contains code to run tsunami
propagation models for each unit source

# Step 3

[TSUNAMI_EVENTS](TSUNAMI_EVENTS) contains code to create uniform and stochastic
slip tsunami events by linearly summing the unit-source tsunami propagation models.
