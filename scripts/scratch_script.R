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

# numerical estimation of continuous model between time step t and t+1
# The output will be the number of individuals just before the next time step
predprey_int = function(t, n, parms) {
  parms$fish.mort <- parms$fish.mort[parms$time]
  parms$bycatch <- parms$bycatch[,parms$time]
  with(as.list(parms), { # extract parameters from parms vector
    dn = rep(0, length(n)) # initialize dn/dt vector
    dn = ( -(M + bycatch*selectivity*fish.mort) * n ) - c(rep(0, 40), rowSums( (t(t(a_ij) %*% diag(n[41:105])) %*% diag(n[1:40])) / ((1 + handling * colSums(t(t(a_ij) %*% diag(n[41:105]))))[col(a_ij)]))) # continuous dynamics
    df = sum((bycatch[1:40]*selectivity[1:40]*fish.mort) * n[1:40]) # calculate fishery yield in number
    dfb = sum((bycatch[1:40]*selectivity[1:40]*fish.mort) * (n[1:40]*weight[1:40])) # calculate fishery yield in biomass
    return(list(dn, df, dfb)) # return dn/dt and df/dt as a list
  })
}

# model --------------------------------------------------------------------------

get_pop_det = function(rockfish, lingcod, init.l = 200, init.r = 70, tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                       f, b, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE) {
  
  # Create empty matrix for data
  nmat = matrix(NA, nrow = 2*lingcod$nage+rockfish$nage, ncol = tf)
  yield_N = numeric(tf)
  yield_B = numeric(tf)
  
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
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, selectivity = selectivity, handling = handling, a_ij = a_ij, weight = weight.at.age)
  id = 1:tf
  
  
# filling nmat in
    n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
    nmat[,1] = n0
    
    for(t in 2:tf) {
      # Calculate recruitment based on previous time step
      SB = nmat[,t-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
      nmat[1,t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Female Recruitment
      nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Male Recruitment
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
    lingcod_SAD = (nmat[1:lingcod$nage,tf] + nmat[(lingcod$nage+1):(2*lingcod$nage), tf])/sum(nmat[1:(2*lingcod$nage),tf]) # Stable age distribution
    rockfish_SAD = nmat[(2*lingcod$nage+1):nrow(nmat),tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
    lingcod_Beq = sum(nmat[1:(2*lingcod$nage),tf] * lingcod$weight.at.age) # total biomass at equilibrium
    rockfish_Beq = sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf] * rockfish$weight.at.age)
    lingcod_oldprop = (nmat[lingcod$nage,tf] + nmat[2*lingcod$nage,tf])/sum(nmat[1:(2*lingcod$nage), tf]) # proportion of pop that are plus age 
    rockfish_oldprop = nmat[nrow(nmat), tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
    yield = rbind(yield_N, yield_B) # total yield 

    # Send output to global environment
    output = list(lingcod_SAD, rockfish_SAD, lingcod_Beq, rockfish_Beq, lingcod_oldprop, rockfish_oldprop, yield)
    names(output) = c("lingcod_SAD", "rockfish_SAD", "lingcod_Beq", "rockfish_Beq", "lingcod_oldprop", "rockfish_oldprop", "yield")
    return(output) # sends as list
}


# Quantify balance between fishing and bycatch in generating unfished biomass and 40% fished biomass --------------

get_fxb_balance = function(rockfish, lingcod, init.l = 200, init.r = 70, tf = (150+20+300), mpa.yr = 20, 
                           hist.f = 0.5, hist.by = 0.5, f, b, a_ij = a_ij.29, handling = 0.01, 
                           times = 1:2, min.selectivity = TRUE) {
  
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
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, selectivity = selectivity, handling = handling, a_ij = a_ij, weight = weight.at.age)
  id = 1:tf
  
  
  # filling nmat in
  n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
  nmat[,1] = n0
  
  for(t in 2:tf) {
    # Calculate recruitment based on previous time step
    SB = nmat[,t-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
    nmat[1,t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Female Recruitment
    nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Male Recruitment
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
  rock_equil_bio = sum(nmat[(2*lingcod$nage+1):nrow(nmat),(tf-1)] * rockfish$weight.at.age)
  
  # Send output to global environment
  return(rock_equil_bio)
}




trial = get_fxb_balance(rockfish, lingcod, init.l = 200, init.r = 70, tf = (150+20+300), mpa.yr = 20, 
  hist.f = 0, hist.by = 0, f = 0, b = 0, a_ij = a_ij.29, handling = 0.01, 
  times = 1:2, min.selectivity = TRUE)



#-------------------------------

library(miceadds)
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish) # with interaction


harvest.slot.trial = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
            f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = FALSE)

# Quantify balance between fishing and bycatch rate ---------------------------------------------------

f = seq(0.01,0.3, by = 0.03)
b = seq(0,0.3, by = 0.01)
harvest.scenarios = as.data.frame(expand_grid(f,b))

fxb_balance = mapply(get_fxb_balance, f = harvest.scenarios[,1], b = harvest.scenarios[,2], 
                     MoreArgs = list(rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE)

fxb_balance_df = cbind(harvest.scenarios, round(unlist(fxb_balance))/round(rockfish_unfished_equil))
names(fxb_balance_df) = c("f", "b", "rockfish_recovery")

fxbplot = fxb_balance_df %>% 
  filter(rockfish_recovery > 0.98) %>% 
  filter(rockfish_recovery < 1.02) %>% 
  ggplot(aes(x = f, y = b)) +
  geom_point() +
  theme_classic()

#lapply -- mclapply
# purrr - perhaps easier but slightly slower map and map2 furrr
# try with a simple example to start
# struct analagous to list

##########################

# Unfished equilibrium
unfished_vals = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
            f = 0, b = 0, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough
lingcod_unfished_equil = unfished_vals$lingcod_Beq
rockfish_unfished_equil = unfished_vals$rockfish_Beq
lingcod_unfished_oldprop = lingcod_oldprop
rockfish_unfished_oldprop = rockfish_oldprop
lingcod_unfished_SAD = lingcod_SAD
rockfish_unfished_SAD = rockfish_SAD

##############################
# struggling with mapply - so here's the forloop version
f = seq(0,0.5, by = 0.1)
b = seq(0,0.5, by = 0.1)
scenarios = unname(as.matrix(expand_grid(f,b)))

# minimum size limit
lingcod_Beq_scenarios = rep(NA, nrow(expand_grid(f,b)))
rockfish_Beq_scenarios = rep(NA, nrow(expand_grid(f,b)))
lingcod_oldprop_scenarios = rep(NA, nrow(expand_grid(f,b)))
rockfish_oldprop_scenarios = rep(NA, nrow(expand_grid(f,b)))

for(s in 1:nrow(scenarios)) {
  get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
              f = scenarios[s,1], b = scenarios[s,2], a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE)
  
  lingcod_Beq_scenarios[s] = lingcod_Beq/lingcod_unfished_equil
  rockfish_Beq_scenarios[s] = rockfish_Beq/rockfish_unfished_equil
  lingcod_oldprop_scenarios[s] = lingcod_oldprop/lingcod_unfished_oldprop
  rockfish_oldprop_scenarios[s] = rockfish_oldprop/rockfish_unfished_oldprop
}

# harvest slot
lingcod_Beq_harvestslot = rep(NA, nrow(expand_grid(f,b)))
rockfish_Beq_harvestslot = rep(NA, nrow(expand_grid(f,b)))
lingcod_oldprop_harvestslot = rep(NA, nrow(expand_grid(f,b)))
rockfish_oldprop_harvestslot = rep(NA, nrow(expand_grid(f,b)))

for(s in 1:nrow(scenarios)) {
  get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
              f = scenarios[s,1], b = scenarios[s,2], a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = FALSE)
  
  lingcod_Beq_harvestslot[s] = lingcod_Beq/lingcod_unfished_equil
  rockfish_Beq_harvestslot[s] = rockfish_Beq/rockfish_unfished_equil
  lingcod_oldprop_harvestslot[s] = lingcod_oldprop/lingcod_unfished_oldprop
  rockfish_oldprop_harvestslot[s] = rockfish_oldprop/rockfish_unfished_oldprop
}

# now combine data and visualize
Beq_scenarios = cbind(scenarios,lingcod_Beq_scenarios, rockfish_Beq_scenarios, lingcod_Beq_harvestslot, rockfish_Beq_harvestslot)
colnames(Beq_scenarios) = c("f", "b", "lingcodMSL", "rockfishMSL", "lingcodHS", "rockfishHS")
save(Beq_scenarios, file = "results/Beq_scenarios_det.Rdata")
oldprop_scenarios = cbind(scenarios,lingcod_oldprop_scenarios, rockfish_oldprop_scenarios, lingcod_oldprop_harvestslot, rockfish_oldprop_harvestslot)
colnames(oldprop_scenarios) = c("f", "b", "lingcodMSL", "rockfishMSL", "lingcodHS", "rockfishHS")
save(oldprop_scenarios, file = "results/oldprop_scenarios_det.Rdata")

library(scales)
library(gridExtra)
# rockfish biomass
grid.arrange(
  ggplot(as.data.frame(Beq_scenarios), aes(x = f, y = b, fill = rockfishMSL)) +
               geom_tile() +
               labs(title = "Rockfish equilibrium biomass", subtitle = "Minimum size limit") +
               theme_classic() + 
               scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 1, limits = c(0,2.5)), 
  ggplot(as.data.frame(Beq_scenarios), aes(x = f, y = b, fill = rockfishHS)) +
    geom_tile() +
    labs(title = "Rockfish equilibrium biomass", subtitle = "Harvest Slot") +
    theme_classic() + 
    scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 1, limits = c(0,2.5)), 
  ncol=2)
# Lingcod biomass
as.data.frame(Beq_scenarios) %>% 
  select("f", "b", "lingcodMSL", "lingcodHS") %>% 
  pivot_longer(cols = starts_with("lingcod"), names_to = "strategy", values_to = "RelativeBiomass") %>% 
  ggplot(aes(x = f, y = RelativeBiomass, color = strategy)) +
    geom_point() +
    theme_classic() +
    labs(title = "Lingcod equilibrium biomass", subtitle = "across fishing and harvest scenarios", y = "B_eq/B_unfished")
# Lingcod old fish
as.data.frame(oldprop_scenarios) %>% 
  select("f", "b", "lingcodMSL", "lingcodHS") %>% 
  pivot_longer(cols = starts_with("lingcod"), names_to = "strategy", values_to = "RelativeOldProp") %>% 
  ggplot(aes(x = f, y = RelativeOldProp, color = strategy)) +
  geom_point() +
  theme_classic() +
  labs(title = "Lingcod long term relative % in plus group", subtitle = "across fishing and harvest scenarios", y = "relative % in plus group")

# rockfish old fish 
grid.arrange(
  ggplot(as.data.frame(oldprop_scenarios), aes(x = f, y = b, fill = rockfishMSL)) + # minimum size limit
    geom_tile() +
    labs(title = "Rockfish long term relative % in plus group", subtitle = "Minimum size limit") +
    theme_classic() + 
    scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 1, limits = c(0,1.4)), 
  ggplot(as.data.frame(oldprop_scenarios), aes(x = f, y = b, fill = rockfishHS)) + # harvest slot
    geom_tile() +
    labs(title = "Rockfish long term relative % in plus group", subtitle = "Harvest slot") +
    theme_classic() + 
    scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 1, limits = c(0,1.4)),
  ncol=2)



# Stochastic --------------------------------------------------------------------------------------

get_pop_stoch = function(rockfish, lingcod, init.l = 200, init.r = 70, nsim = 1, 
                         corr, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1-0.5*cv^2),
                         tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                         f, b, a_ij = a_ij, handling = 0.01, times = 1:2, min.selectivity = TRUE) {
  
  # Create empty matrix for data
  nmat_sims = array(NA, dim = c(2*lingcod$nage+rockfish$nage, tf, nsim), 
                    dimnames=list(spc.age=NULL, 
                                  year=NULL, 
                                  simulation=NULL))
  yield_N = numeric(tf)
  yield_B = numeric(tf)
  
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
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, selectivity = selectivity, handling = handling, a_ij = a_ij, weight = weight.at.age)
  id = 1:tf
  
  
  # filling nmat in
  n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
  nmat[,1] = n0
  
  for(t in 2:tf) {
    # Calculate recruitment based on previous time step
    SB = nmat[,t-1] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
    nmat[1,t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Female Recruitment
    nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Male Recruitment
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
  lingcod_SAD = (nmat[1:lingcod$nage,tf] + nmat[(lingcod$nage+1):(2*lingcod$nage), tf])/sum(nmat[1:(2*lingcod$nage),tf]) # Stable age distribution
  rockfish_SAD = nmat[(2*lingcod$nage+1):nrow(nmat),tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
  lingcod_Beq = sum(nmat[1:(2*lingcod$nage),tf] * lingcod$weight.at.age) # total biomass at equilibrium
  rockfish_Beq = sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf] * rockfish$weight.at.age)
  lingcod_oldprop = (nmat[lingcod$nage,tf] + nmat[2*lingcod$nage,tf])/sum(nmat[1:(2*lingcod$nage), tf]) # proportion of pop that are plus age 
  rockfish_oldprop = nmat[nrow(nmat), tf]/sum(nmat[(2*lingcod$nage+1):nrow(nmat),tf])
  yield = rbind(yield_N, yield_B) # total yield 
  
  # Send output to global environment
  output = list(lingcod_SAD, rockfish_SAD, lingcod_Beq, rockfish_Beq, lingcod_oldprop, rockfish_oldprop, yield)
  names(output) = c("lingcod_SAD", "rockfish_SAD", "lingcod_Beq", "rockfish_Beq", "lingcod_oldprop", "rockfish_oldprop", "yield")
  return(output) # sends as list
}



cv = 0.6










  