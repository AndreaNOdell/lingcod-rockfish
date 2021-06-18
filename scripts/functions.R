# This function simulates an arbitrary number of positive time series that are correlated with one another and are also autocorrelated.
# corr: correlation between the two time series, on the log-scale 
# ar: autoregressive parameter for each time series, on the log-scale
# cv: SD of each time series on the log-scale (approximate CV on the actual scale)
# mn: mean of the time series, on the log-scale. No bias correction is done, so bias correction should be done prior
#     to running this function
# ind_pops: Index of any populations that are independent of the others
sim_correlated_ar_ts <- function(corr, autocorr, cv, mn, nyrs, npops, ind_pops = NULL) {
  sigma.mat <- matrix(0, nrow = npops, ncol = npops)
  diag(sigma.mat) <- cv^2 * (1-autocorr^2)
  for(pop1 in 1:(npops-1)) {
    for(pop2 in (pop1+1):npops) {
      sigma.mat[pop1, pop2] <- sigma.mat[pop2, pop1] <- corr * (1-autocorr[pop1]*autocorr[pop2]) * cv^2 * !(pop1 %in% ind_pops | pop2 %in% ind_pops)
    }
  }
  
  
  burn.in <- 300
  eps <- rmvnorm(nyrs + burn.in, rep(0,npops), sigma.mat) %>%
    as.data.frame()
  
  out.ts <- map2(eps, autocorr, function(.x, .y)
    as.vector(arima.sim(list(ar = .y), nyrs, innov = .x[burn.in+1:nyrs], start.innov = .x[1:burn.in]))) %>% 
    map2(mn, ~ exp(.x + .y)) %>%
    bind_cols() %>%
    as.matrix()
  colnames(out.ts) <- names(mn)
  return(out.ts)
}


# Function to calculate length using age
calc_lengthxage = function(Linf, k, age) {
  length_a = Linf*(1 - exp(-k*age))
  print(length_a)
}

# Function to calculate weight using lengths
# w = aL^b
calc_weightxlength = function(L, a, b) {
  w = a*(L^b)
  print(w)
}

BevHolt = function(phi, h, r0, SBL) ((4*h / (phi*(1-h)))*SBL) / (1 + ((5*h-1) / (phi*r0*(1-h)))*SBL)

# Phi calculation - calculate survival into each age class (L(a) * exp(-M))
phi_calc = function(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE) {
  if(lingcod) {
    phi_n0 = 1 # Initial recruit
    phi_struct = matrix(NA, nrow = 2, ncol = length(age)) # Empty matrix to fill
    phi_struct[,1] = phi_n0 # Fill in initial recruit
    rownames(phi_struct) = c("female", "male") # name rows
    for(ling.sex in c("female", "male")) { # Calculate survival across ages
      for(j in 1:(length(age)-1)) {
        phi_struct[ling.sex,j+1] = phi_struct[ling.sex,j] * exp(-nat.mort[ling.sex])
      }
    }
    # Then multiply each age class by weights and maturity and inverse
    sb_r = phi_struct * weight.at.age * mat.at.age # Spawners per recruit
    phi = 1/rowSums(sb_r) # Recruits per spawner
  } else {
    phi_n0 = 1 # Initial recruit
    phi_struct = numeric(length(age)) # Empty matrix to fill
    phi_struct[1] = phi_n0 # Fill in initial recruit
    # Calculate survival across ages
    for(j in 1:(length(age)-1)) {
      phi_struct[j+1] = phi_struct[j] * exp(-nat.mort)
    }
    # Then multiply each age class by weights and maturity and inverse
    sb_r = phi_struct * weight.at.age * mat.at.age # Spawners per recruit
    phi = 1/sum(sb_r)# Recruits per spawner
    #print(phi)
  }
}