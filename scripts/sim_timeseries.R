library(SimDesign)

corr = 0.8
autocorr = c(0.23, 0.23)
cv = 0.6
mn = 0.5*cv^2
nyrs = 100
npops = 2


simulated_eps = sim_correlated_ar_ts(corr, autocorr, cv, mn, nyrs, npops, ind_pops = NULL)
