# Catalogue duration (years)
D = 100000 

# Exceedance-rate of interest (events/year)
lambda = 1/c(1000, 2500, 10000)

# The number of samples is like a Poisson random variable with rate=lambda*D, 
# and the Monte-Carlo rate is equal to the latter divided by D.
#
# Below we compute 95% confidence intervals on (1/monte_carlo_exceedance_rate)
# for various cases
#
# lambda = 1/1000
1/(qpois(c(0.025, 0.975), lambda[1]*D)/D)
# [1] 1234.567901234568  833.333333333333
#
# lambda = 1/2500
1/(qpois(c(0.025, 0.975), lambda[2]*D)/D)
# [1] 3571.42857142857 1886.79245283019
#
# lambda = 1/10000
1/(qpois(c(0.025, 0.975), lambda[3]*D)/D)
# [1] 25000.00000000000  5882.35294117647

# Generate many random scenarios for 1/2500, to check the result with lambda[2]=1/2500
random_scenarios = rpois(1e+06, lambda=lambda[2]*D)
1/(quantile(random_scenarios, c(0.025, 0.5, 0.975))/D)
#            2.5%              50%            97.5% 
# 3571.42857142857 2500.00000000000 1886.79245283019 
# Notice the above result is pretty-well 'exact' because random_scenarios is full of integers,
# so the quantile gets the answer exactly right.

