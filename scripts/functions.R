library(SimDesign)
# This function simulates an arbitrary number of positive time series that are correlated with one another and are also autocorrelated.
# corr: correlation between the two time series, on the log-scale 
# ar: autoregressive parameter for each time series, on the log-scale
# cv: SD of each time series on the log-scale (approximate CV on the actual scale)
# mn: mean of the time series, on the log-scale. No bias correction is done, so bias correction should be done prior
#     to running this function
# ind_pops: Index of any populations that are independent of the others
sim_correlated_ar_ts <- function(corr, autocorr, cv, mn, tf, npops, ind_pops = NULL) {
  sigma.mat <- matrix(0, nrow = npops, ncol = npops)
  diag(sigma.mat) <- cv^2 * (1-autocorr^2)
  for(pop1 in 1:(npops-1)) {
    for(pop2 in (pop1+1):npops) {
      sigma.mat[pop1, pop2] <- sigma.mat[pop2, pop1] <- corr * (1-autocorr[pop1]*autocorr[pop2]) * cv^2 * !(pop1 %in% ind_pops | pop2 %in% ind_pops)
    }
  }
  
  
  burn.in <- 300
  eps <- rmvnorm(tf + burn.in, rep(0,npops), sigma.mat) %>%
    as.data.frame()
  
  out.ts <- map2(eps, autocorr, function(.x, .y)
    as.vector(arima.sim(list(ar = .y), tf, innov = .x[burn.in+1:tf], start.innov = .x[1:burn.in]))) %>% 
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



# Deterministic burn in 
burn.in = function(lingcod = lingcod, rockfish = rockfish, phi.l = phi.l, phi.r = phi.r, weight.at.age = weight.at.age, 
                   mat.at.age = mat.at.age, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling,
                   a_ij = a_ij.29) {
  
  burn.in.nmat = matrix(NA, nrow = 2*lingcod$nage+rockfish$nage, ncol = 150) # create empty matrix to fill in
  n0 =c(n = c(rep(200, 2*lingcod$nage), rep(70, rockfish$nage))) # Add initial abundances
  burn.in.nmat[,1] = n0
  
  # Parameters for lsoda
  bycatch = as.matrix(c(rep(1,2*lingcod$nage), rep(0.5, rockfish$nage)))
  fish.mort = 0.3
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, lingcod_harvest = lingcod_harvest, 
               selectivity = selectivity, handling = handling, a_ij = a_ij, weight = weight.at.age)
  times = 1:2
  
  for(burn.in in 2:150) {
    # Calculate recruitment based on previous time step
    SB = burn.in.nmat[,burn.in-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
    burn.in.nmat[1,burn.in] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum(SB[1:nage])))/2) # Lingcod Female Recruitment
    burn.in.nmat[with(lingcod, (nage+1)), burn.in] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Male Recruitment
    burn.in.nmat[with(lingcod, (2*nage+1)), burn.in] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(burn.in.nmat)]))) # Rockfish Recruitment
    
    # Numerically integrate abundance after one time step
    parms$time = 1 # assign time step ID to select time-varying parameter for lsoda
    out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
    
    # now move age classes up one step
    burn.in.nmat[2:with(lingcod, nage-1), burn.in] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
    burn.in.nmat[with(lingcod, nage), burn.in] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
    burn.in.nmat[with(lingcod, (nage+2):(2*nage-1)), burn.in] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
    burn.in.nmat[with(lingcod, (2*nage)), burn.in] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
    burn.in.nmat[with(lingcod, (2*nage+2)):(nrow(burn.in.nmat)-1), burn.in] = out[2, with(lingcod, (2*nage+2)):(nrow(burn.in.nmat)-1)] # Rockfish age increase
    burn.in.nmat[nrow(burn.in.nmat), burn.in] = out[2, (nrow(burn.in.nmat))] + out[2,(nrow(burn.in.nmat)+1)] #rockfish plus group
    
    # Now set the vector of abundances as the new n0 for next lsoda run
    n0 = burn.in.nmat[,burn.in]
  }
  output = burn.in.nmat[,150]
  return(output)
}


