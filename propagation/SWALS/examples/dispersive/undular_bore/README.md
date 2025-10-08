# Development of an undular bore in a constant depth channel.

This problem is from an example in [Madsen et al. (2008, their Figure 10)](https://doi.org/10.1029/2008JC004932).

A long channel of constant depth (-20 m) is forced with a sinusoidal wave on the left hand side (amplitude 2 m, period 13 minutes).
As the wave propagates into the channel, nonlinearity causes it to steepen, and in models with dispersion this triggers the formation 
of an undular bore. 

The test problem compares results at a range of grid resolutions with results digitized from Madsen et al. (2008), focussing on the
leading wave about 21 km upstream of the boundary. Madsen et al. (2008) use a more complex dispersive model than currently
in SWALS, but at high resolutions the results are quite similar, although the undular bore wavelength is slightly shorter in SWALS.

An interesting feature of this problem is that for the finite volume schemes,
the resolution needs to be fine enough for noticable undular bore formation to
occur. A similar phenomina (where undular bore formation only occurred at high
model resolution) was noted for the FUNWAVE model in a test by [Yuan et al.
(2020)](https://doi.org/10.1029/2020JB020344).

# Midpoint 

The midpoint solver gives similar results at 2m and 5m resolution, while at 10m
resolution the undular bores have reduced amplitude, and at 20 m resolution
they are barely present.
![Figure 1: Mid point solver](FIXME)

# rk2

The rk2 solver is expected to be lower-order-accuracy than midpoint when
dispersion matters, because it employs separate shallow water and dispersive
steps (whereas midpoint includes dispersion in the predictor as well as the
corrector step). 
![Figure 2: rk2 solver](FIXME)

# Leapfrog nonlinear

This solver includes first-order-accurate momentum advection, and since
momentum advection matters for this problem, at high resolution it is expected
to be less accurate than the finite volume schemes. Interestingly at low
resolution this solver does a better job of reflecting the presence of bores,
likely because the leapfrog nonlinear solver does not include fancy treatments
of dissipation near discontinuities (i.e. limiters).

![Figure 3: Leapfrog Nonlinear solver](FIXME)

# Cliffs

This seems slower to converge at high resolution, likely reflecting limitations of
how dispersion is integrated into the CLIFFS timestepping.

![Figure 4: Cliffs solver](FIXME)
