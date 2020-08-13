The script here compares a scenario using the leapfrog_nonlinear solver offshore, vs the same scenario using linear_with_manning. 

* The model with leapfrog_nonlinear  offshore took 21578.58134782 seconds [FIRST VERSION] and 21269.0753 [REVISED NESTING]. The answers were not visually different, but there were minor differences higher decimal places.
* The model with linear_with_manning offshore took  6940.42678905 seconds

These were both done using 8 nodes of Gadi with the same compiler options.


