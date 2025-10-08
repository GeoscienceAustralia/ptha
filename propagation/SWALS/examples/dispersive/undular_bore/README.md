# Development of an undular bore in a constant depth channel.

This problem is from an example in [Madsen et al. (2008, their Figure 10)](https://doi.org/10.1029/2008JC004932).

A long channel of constant depth (-20 m) is forced with a sinusoidal wave on the left hand side (amplitude 2 m, period 13 minutes).
As the wave propagates into the channel, nonlinearity causes it to steepen, and in models with dispersion this triggers the formation 
of an undular bore. 

The test problem compares results at a range of grid resolutions with results digitized from Madsen et al. (2008), focussing on the
leading wave about 21 km upstream of the boundary. Madsen et al. (2008) use a more complex dispersive model than currently
in SWALS, but at high resolutions the results are quite similar, although the undular bore wavelenght is slightly shorter in SWALS.

# Midpoint 

The midpoint solver gives similar results at 2m and 5m resoltion, while at 10m resolution the undular bores have reduced amplitude, and at
20 m resolution there are no bores at all.
![Figure 1: Mid point solver](FIXME)

# rk2

The rk2 solver is expected to be lower-order-accuracy than midpoint when dispersion matters, because it employs separate shallow water and dispersive steps (whereas midpoint includes dispersion in the predictor as well as the corrector step). 
![Figure 2: rk2 solver](FIXME)

# Leapfrog nonlinear

![Figure 3: Leapfrog Nonlinear solver](FIXME)

# Cliffs

![Figure 4: Cliffs solver](FIXME)
