# equilibrium age-structure (N_a/N_total i.e. the number of age a individuals divided by total individuals) for each fishing x bycatch scenario for lingcod and rockfish
# proportion of older fish (ages 12-20 - I should think more about this range) N[12-20]/sum(N) across fishing x bycatch scenarios for lingcod and rockfish
# total yield in number and biomass as a function across fishing x bycatch scenarios for lingcod only
# Comparison of biomass to unfished biomass (B_equil/B_unfished) to rockfish and lingcod

# Set up -------------------------------------------------------------------------
library(deSolve)
library(tidyverse)
source("scripts/functions.R")
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")

# numerical estimation of continuous model between time step t and t+1 ---------------------
# The output will be the number of individuals just before the next time step
predprey_int = function(t, n, parms) {
  parms$fish.mort <- parms$fish.mort[parms$time]
  parms$bycatch <- parms$bycatch[,parms$time]
  with(as.list(parms), { # extract parameters from parms vector
    dn = rep(0, length(n)) # initialize dn/dt vector
    dn = -(M + bycatch*selectivity*c(rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40), rep(fish.mort,rockfish$nage))) * n - c(rep(0, 40), rowSums(t(t((-log(1 - (consump_yelloweye_n / n[41:105])) * n[41:105]) %*% diag(n[1:40])) / (1 + handling * colSums(-log(1 - (consump_yelloweye_n / n[41:105])) * n[41:105]) + handling*otherprey_n)))) # continuous dynamics
    dc = rowSums(t(t((-log(1 - (consump_yelloweye_n / n[41:105])) * n[41:105]) %*% diag(n[1:40])) / (1 + handling * colSums(-log(1 - (consump_yelloweye_n / n[41:105])) * n[41:105]) + handling*otherprey_n))) * rockfish$weight.at.age  # calculate predation in biomass
    df = sum(bycatch[1:40]*selectivity[1:40]*rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40) * n[1:40]) # calculate fishery yield in number
    dfb = sum((bycatch[1:40]*selectivity[1:40]*rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40)) * (n[1:40]*weight[1:40])) # calculate fishery yield in biomass
    return(list(dn, dc, df, dfb)) # return dn/dt and df/dt as a list
  })
}


# Deterministic model --------------------------------------------------------------------------

get_pop_det = function(rockfish, lingcod, init.l = 200, init.r = 70, tf = (150+20+300), mpa.yr = 20, hist.f = 0.3, hist.by = 0.5, 
                       f, b, a_ij = binned.size.spec.29, handling = 0.7, times = 1:2, min.selectivity = TRUE,
                       min.age.consumed = 4, rockfish.prop = 0.03, yelloweye.prop = 0.05) {
  
  # Create empty matrix for data
  nmat = matrix(NA, nrow = 2*lingcod$nage+rockfish$nage, ncol = tf)
  yield_N = numeric(tf)
  yield_B = numeric(tf)
  
  # Create weight and maturity at age vectors
  weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
  mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
  
  # Selecting harvest scenario
  if(min.selectivity) {
    lingcod_harvest = lingcod$min.selectivity
  }  else {
    lingcod_harvest = lingcod$harvest.slot
  }
  
  # Parameters for lsoda
  M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                               rep(nat.mort["male"],nage))), # Lingcod male natural mortality
               with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
  post.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                         rep(1,nage))), # Lingcod male bycatch
                                         with(rockfish, rep(b,nage)))), (tf-150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by, nage)))), (150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod_harvest, # Lingcod fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr))))
  # diet 
  min.age.consumed = min.age.consumed
  rockfish.prop = rockfish.prop
  yelloweye.prop = yelloweye.prop
  diet.frac.yelloweye <- rep(c(seq(0, rockfish.prop, length = min.age.consumed)*yelloweye.prop, rep(rockfish.prop*yelloweye.prop, (lingcod$nage-4))), 2)
  # Lingcod total annual consumption per capita at age (kg) distributed across rockfish ages using size-spectrum 
  # of preferences. Resulting in age-specific annual consumption for each lingcod age and sex.
  consump_tot_kg = a_ij %*% diag(lingcod$consump.at.age)
  # age/sex-specific annual lingcod consumption that is specifically yelloweye rockfish (kg) 
  consump_yelloweye_kg = (consump_tot_kg %*% diag(diet.frac.yelloweye))
  # Convert kg (biomass) to number (abundance)
  consump_yelloweye_n = sweep(consump_yelloweye_kg, MARGIN = 1, FUN = "/", STATS = rockfish$weight.at.age)
  
  otherprey_n = sum(consump_yelloweye_n) / (yelloweye.prop*rockfish.prop)
  
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling, consump_yelloweye_n = consump_yelloweye_n, weight = weight.at.age, otherprey_n = otherprey_n)
  id = 1:tf
  
  
# filling nmat in
    n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
    nmat[,1] = n0
    
    for(t in 2:tf) {
      # Calculate recruitment based on previous time step
      SB = nmat[,t-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
      nmat[1,t] = nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l, h, r0, sum((SB[1:nage]))))/2) # Lingcod Recruitment
      nmat[with(lingcod, (2*nage+1)), t] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat)]))) # Rockfish Recruitment
      
      # Numerically integrate abundance after one time step
      parms$time = id[t-1] # assign time step ID to select time-varying parameter for lsoda
      out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
      parms$fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr)))) # re-enter the vector of fish.mort
      parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch) # re-enter the vector of bycatch
      
      # now move age classes up one step
      nmat[2:with(lingcod, nage-1), t] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
      nmat[with(lingcod, nage), t] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
      nmat[with(lingcod, (nage+2):(2*nage-1)), t] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
      nmat[with(lingcod, (2*nage)), t] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
      nmat[with(lingcod, (2*nage+2)):(nrow(nmat)-1), t] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat)-1)] # Rockfish age increase
      nmat[nrow(nmat), t] = out[2, (nrow(nmat))] + out[2,(nrow(nmat)+1)] #rockfish plus group
      
      # Add fishery yield into fishery yield vector
      yield_N[(t-1)] = out[1, (ncol(out)-1)]
      yield_B[(t-1)] = out[1, ncol(out)]
      
      # Now set the vector of abundances as the new n0 for next lsoda run
      n0 = nmat[,t]
    }
    
    # Collecting and saving the necessary information
      #lingcod_SAD = (nmat[1:lingcod$nage,tf] + nmat[(lingcod$nage+1):(2*lingcod$nage), tf])/sum(nmat[1:(2*lingcod$nage),tf]) # Stable age distribution
      #rockfish_SAD = nmat[(2*lingcod$nage+1):nrow(nmat),tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
    lingcod_SBeq = sum(nmat[1:(lingcod$nage),tf] * lingcod$weight.at.age[1:(lingcod$nage)] * lingcod$mat.at.age[1:(lingcod$nage)]) # total Spawning biomass at equilibrium
    rockfish_SBeq = sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf] * rockfish$weight.at.age * rockfish$mat.at.age)
    lingcod_oldprop = (nmat[lingcod$nage,tf] + nmat[2*lingcod$nage,tf])/sum(nmat[1:(2*lingcod$nage), tf]) # proportion of pop that are plus age 
    rockfish_oldprop = nmat[nrow(nmat), tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
    yield = rbind(yield_N, yield_B) # total yield 
    biomass_rockfish = colSums(nmat[((2*lingcod$nage)+1):nrow(nmat), ] * rockfish$weight.at.age)

    # Send output to global environment
    output = list(lingcod_SBeq, rockfish_SBeq, lingcod_oldprop, rockfish_oldprop, yield, biomass_rockfish)
    names(output) = c("lingcod_SBeq", "rockfish_SBeq", "lingcod_oldprop", "rockfish_oldprop", "yield", "biomass_rockfish")
    return(output) # sends as list
}


# Quantify balance between fishing and bycatch in generating unfished biomass and 40% fished biomass --------------

get_fxb_balance = function(rockfish, lingcod, init.l = 200, init.r = 70, tf = (150+20+300), mpa.yr = 20, 
                           hist.f = 0.3, hist.by = 0.5, f, b, a_ij = a_ij.29, handling = 0.7, 
                           times = 1:2, min.selectivity = TRUE, min.age.consumed = 4, rockfish.prop = 0.03, 
                           yelloweye.prop = 0.05) {
  
  # Create empty matrix for data
  nmat = matrix(NA, nrow = 2*lingcod$nage+rockfish$nage, ncol = tf)
  
  # Calculate phi and create vector
  phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))[1]
  phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))
  # Create weight and maturity at age vectors
  weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
  mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
  
  # Selecting harvest scenario
  if(min.selectivity) {
    lingcod_harvest = lingcod$min.selectivity
  }  else {
    lingcod_harvest = lingcod$harvest.slot
  }
  
  # Parameters for lsoda
  M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                               rep(nat.mort["male"],nage))), # Lingcod male natural mortality
               with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
  post.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                         rep(1,nage))), # Lingcod male bycatch
                                         with(rockfish, rep(b,nage)))), (tf-150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by, nage)))), (150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod_harvest, # Lingcod fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr))))
  # diet 
  min.age.consumed = min.age.consumed
  rockfish.prop = rockfish.prop
  yelloweye.prop = yelloweye.prop
  diet.frac.yelloweye <- rep(c(seq(0, rockfish.prop, length = min.age.consumed)*yelloweye.prop, rep(rockfish.prop*yelloweye.prop, (lingcod$nage-4))), 2)
  # Lingcod total annual consumption per capita at age (kg) distributed across rockfish ages using size-spectrum 
  # of preferences. Resulting in age-specific annual consumption for each lingcod age and sex.
  consump_tot_kg = a_ij %*% diag(lingcod$consump.at.age)
  # age/sex-specific annual lingcod consumption that is specifically yelloweye rockfish (kg) 
  consump_yelloweye_kg = (consump_tot_kg %*% diag(diet.frac.yelloweye))
  # Convert kg (biomass) to number (abundance)
  consump_yelloweye_n = sweep(consump_yelloweye_kg, MARGIN = 1, FUN = "/", STATS = rockfish$weight.at.age)
  
  otherprey_n = sum(consump_yelloweye_n) / (yelloweye.prop*rockfish.prop)
  
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling, consump_yelloweye_n = consump_yelloweye_n, weight = weight.at.age, otherprey_n = otherprey_n)
  id = 1:tf
  
  # filling nmat in
  n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
  nmat[,1] = n0
  
  for(t in 2:tf) {
    # Calculate recruitment based on previous time step
    SB = nmat[,t-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
    nmat[1,t] = nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Recruitment
    nmat[with(lingcod, (2*nage+1)), t] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat)]))) # Rockfish Recruitment
    
    # Numerically integrate abundance after one time step
    parms$time = id[t-1] # assign time step ID to select time-varying parameter for lsoda
    out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
    parms$fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr)))) # re-enter the vector of fish.mort
    parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch) # re-enter the vector of bycatch
    
    # now move age classes up one step
    nmat[2:with(lingcod, nage-1), t] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
    nmat[with(lingcod, nage), t] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
    nmat[with(lingcod, (nage+2):(2*nage-1)), t] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
    nmat[with(lingcod, (2*nage)), t] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
    nmat[with(lingcod, (2*nage+2)):(nrow(nmat)-1), t] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat)-1)] # Rockfish age increase
    nmat[nrow(nmat), t] = out[2, (nrow(nmat))] + out[2,(nrow(nmat)+1)] #rockfish plus group
    
    # Now set the vector of abundances as the new n0 for next lsoda run
    n0 = nmat[,t]
  }
  
  # Collecting and saving the necessary information
  rock_equil_bio = sum(nmat[(2*lingcod$nage+1):nrow(nmat),(tf-1)] * rockfish$weight.at.age * rockfish$mat.at.age)
  # Send output to global environment
  return(rock_equil_bio)
}




# Stochastic --------------------------------------------------------------------------------------
# calculates the average and coefficient of variation  
get_pop_stoch = function(rockfish, lingcod, nsim, corr, autocorr = c(0.23,0.23), cv = 0.6,
                         tf = (30+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                         f, b, a_ij = binned.size.spec.29, handling = 2, times = 1:2, min.selectivity = TRUE,
                         min.age.consumed = 4, rockfish.prop = 0.03, yelloweye.prop = 0.05) {
  
  # Create empty matrix for data
  nmat_sims = array(NA, dim = c(2*lingcod$nage+rockfish$nage, tf, nsim), 
                    dimnames=list(spc.age=NULL, 
                                  year=NULL, 
                                  simulation=NULL))
  lingcod_biomass_ts = rockfish_biomass_ts = matrix(NA, nrow = nsim, ncol = tf)
  lingcod_consumption_ts = matrix(NA, nrow = nsim, ncol = tf)
  
  # Create weight and maturity at age vectors
  weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
  mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
  
  # Selecting harvest scenario
  if(min.selectivity) {
    lingcod_harvest = lingcod$min.selectivity
  }  else {
    lingcod_harvest = lingcod$harvest.slot
  }
  
  
  # Parameters for lsoda
  M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                               rep(nat.mort["male"],nage))), # Lingcod male natural mortality
               with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
  post.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                         rep(1,nage))), # Lingcod male bycatch
                                         with(rockfish, rep(b,nage)))), (tf-30)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by, nage)))), (30)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod_harvest, # Lingcod fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 30), rep(0, mpa.yr), rep(f, (tf-(30+mpa.yr))))
  # diet 
  min.age.consumed = min.age.consumed
  rockfish.prop = rockfish.prop
  yelloweye.prop = yelloweye.prop
  diet.frac.yelloweye <- rep(c(seq(0, rockfish.prop, length = min.age.consumed)*yelloweye.prop, rep(rockfish.prop*yelloweye.prop, (lingcod$nage-4))), 2)
  # Lingcod total annual consumption per capita at age (kg) distributed across rockfish ages using size-spectrum 
  # of preferences. Resulting in age-specific annual consumption for each lingcod age and sex.
  consump_tot_kg = a_ij %*% diag(lingcod$consump.at.age)
  # age/sex-specific annual lingcod consumption that is specifically yelloweye rockfish (kg) 
  consump_yelloweye_kg = (consump_tot_kg %*% diag(diet.frac.yelloweye))
  # Convert kg (biomass) to number (abundance)
  consump_yelloweye_n = sweep(consump_yelloweye_kg, MARGIN = 1, FUN = "/", STATS = rockfish$weight.at.age)
  
  otherprey_n = sum(consump_yelloweye_n) / (yelloweye.prop*rockfish.prop)
  
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling, consump_yelloweye_n = consump_yelloweye_n, weight = weight.at.age, otherprey_n = otherprey_n)
  id = 1:tf
  
  
  det.burn.in = burn.in(lingcod = lingcod, rockfish = rockfish, phi.l = phi.l, phi.r = phi.r, weight.at.age = weight.at.age, M = M,
                        mat.at.age = mat.at.age, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling,
                        a_ij = a_ij, consump_yelloweye_n = consump_yelloweye_n, otherprey_n = otherprey_n)
  
  for(i in 1:nsim) {
    
    # filling nmat_sims in
    nmat_sims[,,i] = NA
    n0 = n = det.burn.in # vector of initial abundances  
    nmat_sims[,1,i] = n0
    mn = 1-(0.5*cv^2)
    list.of.seeds = 1:nsim
    set.seed(list.of.seeds[i])
    corr_eps = sim_correlated_ar_ts(corr, autocorr, cv = log(cv), mn = log(mn), tf, npops = 2, ind_pops = NULL) # each simulation needs a new stochastic recrtuitment time series
    
    for(t in 2:tf) {
      # Calculate recruitment based on previous time step
      SB = nmat_sims[,t-1, i] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
      nmat_sims[1,t, i] = nmat_sims[with(lingcod, (nage+1)), t, i] = with(lingcod, (BevHolt(phi.l, h, r0, sum(SB[1:nage]))*corr_eps[t,1])/2) # Lingcod Recruitment
      nmat_sims[with(lingcod, (2*nage+1)), t, i] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat_sims)]))*corr_eps[t,2]) # Rockfish Recruitment
      
      # Numerically integrate abundance after one time step
      parms$time = id[t-1] # assign time step ID to select time-varying parameter for lsoda
      out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
      parms$fish.mort = c(rep(hist.f, 30), rep(0, mpa.yr), rep(f, (tf-(30+mpa.yr)))) # re-enter the vector of fish.mort
      parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch) # re-enter the vector of bycatch
      
      # now move age classes up one step
      nmat_sims[2:with(lingcod, nage-1), t, i] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
      nmat_sims[with(lingcod, nage), t, i] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
      nmat_sims[with(lingcod, (nage+2):(2*nage-1)), t, i] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
      nmat_sims[with(lingcod, (2*nage)), t, i] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
      nmat_sims[with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1), t, i] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1)] # Rockfish age increase
      nmat_sims[nrow(nmat_sims), t, i] = out[2, (nrow(nmat_sims))] + out[2,(nrow(nmat_sims)+1)] #rockfish plus group
      
      # Save lingcod consumption information from ODE
      lingcod_consumption_ts[i, (t-1)] = out[1, (ncol(out)-2)]
      
      # Now set the vector of abundances as the new n0 for next lsoda run
      n0 = nmat_sims[,t,i]
    }
    
    # Collecting and saving the necessary information
    lingcod_biomass_ts[i,] = colSums(nmat_sims[1:lingcod$nage,,i]*lingcod$weight.at.age[1:lingcod$nage] * lingcod$mat.at.age[1:lingcod$nage]) # time series of spawningbiomass (columns) across simulations (rows)
    rockfish_biomass_ts[i,] = colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i]*rockfish$weight.at.age * rockfish$mat.at.age)
 }
  
  # Calculate cv of biomass timeseries for each simulation and take the mean of cv (i.e average cv across simulation)
    #lingcod_cv = mean(apply(lingcod_biomass_ts[,-(1:100)], 1, function(x) sd(x) / mean(x) * 100))
    #lingcod_avg = mean(apply(lingcod_biomass_ts[,-(1:100)], 1, mean))
  rockfish_cv = mean(apply(rockfish_biomass_ts[,-(1:100)], 1, function(x) sd(x) / mean(x) * 100))
  rockfish_avg = mean(apply(rockfish_biomass_ts[,-(1:100)], 1, mean))
  lingcod_consumption_ts = lingcod_consumption_ts[, -(ncol(lingcod_consumption_ts))]
  consumption_cv = mean(apply(lingcod_consumption_ts[,-(1:100)], 1, function(x) sd(x) / mean(x) * 100))
  consumption_avg = mean(apply(lingcod_consumption_ts[,-(1:100)], 1, mean))
  
  # Send output to global environment
  output = list(rockfish_cv, rockfish_avg, consumption_cv, consumption_avg)
  names(output) = c("rockfish_cv", "rockfish_avg", "consumption_cv", "consumption_avg")
  return(output) # sends as list
}



# Extract random simulations of one fishing x bycatch scenario but contrasting -----------
# recruitment responses between populations.
get_sample_ts = function(rockfish, lingcod, nsim, corr, autocorr = c(0.23,0.23), cv = 0.6,
                         tf = (30+20+400), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                         f = 0.1, b = 0.1, a_ij = binned.size.spec.29, handling = 0.7, times = 1:2, 
                         min.selectivity = TRUE, min.age.consumed = 4, rockfish.prop = 0.03, yelloweye.prop = 0.05) {
  
  # Create empty matrix for data
  nmat_sims = array(NA, dim = c(2*lingcod$nage+rockfish$nage, tf, nsim), 
                    dimnames=list(spc.age=NULL, 
                                  year=NULL, 
                                  simulation=NULL))
  lingcod_biomass_ts = rockfish_biomass_ts = matrix(NA, nrow = nsim, ncol = tf)
  lingcod_consumption_ts = matrix(NA, nrow = nsim, ncol = tf)
  lingcod_SB_ts = rockfish_SB_ts = matrix(NA, nrow = nsim, ncol = tf)
  
  # Create weight and maturity at age vectors
  weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
  mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
  
  # Selecting harvest scenario
  if(min.selectivity) {
    lingcod_harvest = lingcod$min.selectivity
  }  else {
    lingcod_harvest = lingcod$harvest.slot
  }
  
  # Parameters for lsoda
  M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                               rep(nat.mort["male"],nage))), # Lingcod male natural mortality
               with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
  post.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                         rep(1,nage))), # Lingcod male bycatch
                                         with(rockfish, rep(b,nage)))), (tf-30)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by, nage)))), (30)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod_harvest, # Lingcod fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 30), rep(0, mpa.yr), rep(f, (tf-(30+mpa.yr))))
  # diet 
  min.age.consumed = min.age.consumed
  rockfish.prop = rockfish.prop
  yelloweye.prop = yelloweye.prop
  diet.frac.yelloweye <- rep(c(seq(0, rockfish.prop, length = min.age.consumed)*yelloweye.prop, rep(rockfish.prop*yelloweye.prop, (lingcod$nage-4))), 2)
  # Lingcod total annual consumption per capita at age (kg) distributed across rockfish ages using size-spectrum 
  # of preferences. Resulting in age-specific annual consumption for each lingcod age and sex.
  consump_tot_kg = a_ij %*% diag(lingcod$consump.at.age)
  # age/sex-specific annual lingcod consumption that is specifically yelloweye rockfish (kg) 
  consump_yelloweye_kg = (consump_tot_kg %*% diag(diet.frac.yelloweye))
  # Convert kg (biomass) to number (abundance)
  consump_yelloweye_n = sweep(consump_yelloweye_kg, MARGIN = 1, FUN = "/", STATS = rockfish$weight.at.age)
  
  otherprey_n = sum(consump_yelloweye_n) / (yelloweye.prop*rockfish.prop)
  
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling, consump_yelloweye_n = consump_yelloweye_n, weight = weight.at.age, otherprey_n = otherprey_n)
  id = 1:tf
  
  det.burn.in = burn.in(lingcod = lingcod, rockfish = rockfish, phi.l = phi.l, phi.r = phi.r, weight.at.age = weight.at.age, M = M,
                        mat.at.age = mat.at.age, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling,
                        a_ij = a_ij, consump_yelloweye_n = consump_yelloweye_n, otherprey_n, otherprey_n)
  
  for(i in 1:nsim) {
    
    # filling nmat_sims in
    nmat_sims[,,i] = NA
    n0 = n = det.burn.in # vector of initial abundances  
    nmat_sims[,1,i] = n0
    mn = 1-(0.5*cv^2)
    corr_eps = sim_correlated_ar_ts(corr, autocorr, cv = log(cv), mn = log(mn), tf, npops = 2, ind_pops = NULL) # each simulation needs a new stochastic recrtuitment time series
    
    for(t in 2:tf) {
      # Calculate recruitment based on previous time step
      SB = nmat_sims[,t-1, i] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
      nmat_sims[1,t, i] = nmat_sims[with(lingcod, (nage+1)), t, i] = with(lingcod, (BevHolt(phi.l, h, r0, sum(SB[1:nage]))*corr_eps[t,1])/2) # Lingcod Female Recruitment
      nmat_sims[with(lingcod, (2*nage+1)), t, i] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat_sims)]))*corr_eps[t,2]) # Rockfish Recruitment
      
      # Numerically integrate abundance after one time step
      parms$time = id[t-1] # assign time step ID to select time-varying parameter for lsoda
      out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
      parms$fish.mort = c(rep(hist.f, 30), rep(0, mpa.yr), rep(f, (tf-(30+mpa.yr)))) # re-enter the vector of fish.mort
      parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch) # re-enter the vector of bycatch
      
      # now move age classes up one step
      nmat_sims[2:with(lingcod, nage-1), t, i] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
      nmat_sims[with(lingcod, nage), t, i] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
      nmat_sims[with(lingcod, (nage+2):(2*nage-1)), t, i] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
      nmat_sims[with(lingcod, (2*nage)), t, i] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
      nmat_sims[with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1), t, i] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1)] # Rockfish age increase
      nmat_sims[nrow(nmat_sims), t, i] = out[2, (nrow(nmat_sims))] + out[2,(nrow(nmat_sims)+1)] #rockfish plus group
      
      # Save lingcod consumption information from ODE
      lingcod_consumption_ts[i, (t-1)] = out[1, (ncol(out)-2)]
      
      # Now set the vector of abundances as the new n0 for next lsoda run
      n0 = nmat_sims[,t,i]
    }
    
    # Collecting and saving the necessary information
    lingcod_SB_ts[i,] = colSums(nmat_sims[1:lingcod$nage,,i]*lingcod$weight.at.age[1:lingcod$nage] * lingcod$mat.at.age[1:lingcod$nage]) # time series of spawningbiomass (columns) across simulations (rows)
    rockfish_SB_ts[i,] = colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i]*rockfish$weight.at.age * rockfish$mat.at.age)
    
    lingcod_biomass_ts[i,] = colSums(nmat_sims[1:(2*lingcod$nage),,i]*lingcod$weight.at.age)
    rockfish_biomass_ts[i,] = colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i]*rockfish$weight.at.age)
  }
  
  # Extract random simulations
  sim_sample = sample(1:nsim, 3)
  SB_ts_1 = rbind(lingcod_SB_ts[sim_sample[1],], rockfish_SB_ts[sim_sample[1],])
  SB_ts_2 = rbind(lingcod_SB_ts[sim_sample[2],], rockfish_SB_ts[sim_sample[2],])
  SB_ts_3 = rbind(lingcod_SB_ts[sim_sample[3],], rockfish_SB_ts[sim_sample[3],])
  biomass_ts_1 = rbind(lingcod_biomass_ts[sim_sample[1],], rockfish_biomass_ts[sim_sample[1],])
  biomass_ts_2 = rbind(lingcod_biomass_ts[sim_sample[2],], rockfish_biomass_ts[sim_sample[2],])
  biomass_ts_3 = rbind(lingcod_biomass_ts[sim_sample[3],], rockfish_biomass_ts[sim_sample[3],])
  lingcod_consumption = lingcod_consumption_ts[sim_sample,]
  age_ts_1 = nmat_sims[,,sim_sample[1]]
  age_ts_2 = nmat_sims[,,sim_sample[2]]
  age_ts_3 = nmat_sims[,,sim_sample[3]]
    
  
  # Send output to global environment
  output = list(SB_ts_1, SB_ts_2, SB_ts_3, biomass_ts_1, biomass_ts_2, biomass_ts_3, 
                lingcod_consumption, age_ts_1, age_ts_2, age_ts_3)
  names(output) = c("SB_ts_1", "SB_ts_2", "SB_ts_3", "biomass_ts_1", "biomass_ts_2", "biomass_ts_3", 
                    "lingcod_consumption", "age_ts_1", "age_ts_2", "age_ts_3")
  return(output) # sends as list
}



# trial_pos = get_sample_ts(rockfish, lingcod, nsim = 10, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6,
#                           tf = (30+20+400), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
#                           f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE)
# 
# trial_neg = get_sample_ts(rockfish, lingcod, nsim = 10, corr = -0.8, autocorr = c(0.23,0.23), cv = 0.6,
#                           tf = (30+20+400), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
#                           f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE)
# 
# plot(1:450, trial_neg$ts_3[2,], type = "l", ylim = c(0,3500), xlab = "time", ylab = "rockfish SB", main = "negative corr timeseries")
# lines(1:450, trial_neg$ts_1[2,], col = "blue")
# lines(1:450, trial_neg$ts_2[2,], col = "red")
# lines(1:450, trial_neg$ts_4[2,], col = "green")
# lines(1:450, trial_neg$ts_5[2,], col = "orange")

# plot(1:450, trial_neg$ts_3[1,], type = "l", xlab = "time", ylim = c(0,43000), ylab = "lingcod SB", main = "negative corr timeseries")
# lines(1:450, trial_neg$ts_1[1,], col = "blue")
# lines(1:450, trial_neg$ts_2[1,], col = "red")
# lines(1:450, trial_neg$ts_4[1,], col = "green")
# lines(1:450, trial_neg$ts_5[1,], col = "orange")


  