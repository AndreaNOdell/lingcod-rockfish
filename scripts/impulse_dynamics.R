
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
    df = sum((bycatch[1:40]*selectivity[1:40]*fish.mort) * n[1:40]) # calculate fishery yield
    dfb = sum((bycatch[1:40]*selectivity[1:40]*fish.mort) * (n[1:40]*weight[1:40]))
    return(list(dn, df, dfb)) # return dn/dt and df/dt as a list
  })
}


# model ------------------------------------------------------------------------
cv = 0.6

# Create empty matrix

get_pop_n = function(rockfish, lingcod, nsim, init.l, init.r, corr, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1-0.5*cv^2), 
                     tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = 0.1, a_ij = a_ij, handling = 0.01, times = 1:2, 
                     harvest.strat = min.selectivity, stochastic = TRUE, biomass = FALSE) {
  # Create empty matrix
  nmat_sims = array(NA, dim = c(2*lingcod$nage+rockfish$nage, tf, nsim), 
                    dimnames=list(spc.age=NULL, 
                                  year=NULL, 
                                  simulation=NULL))
  lingcod_tot = rockfish_tot = matrix(NA, nrow = nsim, ncol = tf)
  lingcod_tot_B = rockfish_tot_B = matrix(NA, nrow = nsim, ncol = tf)
  fishery_yield_sims_N = matrix(NA, nrow = nsim, ncol = tf)
  fishery_yield_sims_B = matrix(NA, nrow = nsim, ncol = tf)
  
  # Calculate phi and create vector
  phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))[1]
  phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))
  # Create weight and maturity at age vectors
  weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
  mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
  # Parameters for lsoda
  M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                               rep(nat.mort["male"],nage))), # Lingcod male natural mortality
               with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
  post.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                     rep(1,nage))), # Lingcod male bycatch
                     with(rockfish, rep(b,nage)))), (tf-150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by ,nage)))), (150)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod$harvest.strat, # Lingcod female fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr))))
  parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, selectivity = selectivity, handling = handling, a_ij = a_ij, weight = weight.at.age)
  id = 1:tf
  
  for(i in 1:nsim) {
   if(stochastic) { # for stochastic recruitment. multiply by autoregressive lognormally distributed deviations
   corr_eps = sim_correlated_ar_ts(corr, autocorr, cv, mn, tf, npops = 2, ind_pops = NULL) # each simulation needs a new stochastic recrtuitment time series
   } else { # for deterministic recruitment. multiply by 1
     corr_eps = matrix(1, nrow = tf, ncol = 2) 
   }
      nmat_sims[,,i] = NA
      n0 = c(n = c(rep(init.l, 2*lingcod$nage), rep(init.r, rockfish$nage))) # vector of initial abundances  
      nmat_sims[,1,i] = n0
    
   for(t in 2:tf) {
      # Calculate recruitment based on previous time step
        # Caclulate spawning biomass
          SB = nmat_sims[,t-1, i] * weight.at.age * mat.at.age
          nmat_sims[1,t, i] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage])))*corr_eps[t,1])/2) # Lingcod Female Recruitment
          nmat_sims[with(lingcod, (nage+1)), t, i] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage])))*corr_eps[t,1])/2) # Lingcod Male Recruitment
          nmat_sims[with(lingcod, (2*nage+1)), t, i] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat_sims)]))*corr_eps[t,2]) # Rockfish Recruitment
      
      # Numerically integrate abundance after one time step
        parms$time = id[t-1]
        out = as.matrix(lsoda(n0, times, predprey_int, parms))
        parms$fish.mort = c(rep(hist.f, 150), rep(0, mpa.yr), rep(f, (tf-(150+mpa.yr))))
        parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
      
      # now move age classes up one step
        nmat_sims[2:with(lingcod, nage-1), t, i] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
        nmat_sims[with(lingcod, nage), t, i] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
        nmat_sims[with(lingcod, (nage+2):(2*nage-1)), t, i] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
        nmat_sims[with(lingcod, (2*nage)), t, i] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
        nmat_sims[with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1), t, i] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1)] # Rockfish age increase
        nmat_sims[nrow(nmat_sims), t, i] = out[2, (nrow(nmat_sims))] + out[2,(nrow(nmat_sims)+1)] #rockfish plus group
      
      # Add fishery yield into fishery yield matrix
        fishery_yield_sims_N[i, (t-1)] = out[1, (ncol(out)-1)]
        fishery_yield_sims_B[i, (t-1)] = out[1, ncol(out)]
        
      # Now set the vector of abundances as the new n0 for next lsoda run
        n0 = nmat_sims[,t,i]
   }
    # Now take the total lingcod and rockfish abundances each year (collapse structure)
      lingcod_tot[i,] = colSums(nmat_sims[1:with(lingcod, 2*nage),,i])
      rockfish_tot[i,] = colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i])
      
      if(biomass) {
        lingcod_tot_B[i,] = colSums(nmat_sims[1:(2*lingcod$nage),,i] * weight.at.age[1:(2*lingcod$nage)])
        rockfish_tot_B[i,] = colSums(nmat_sims[(2*lingcod$nage+1):nrow(nmat_sims),,i] * weight.at.age[(2*lingcod$nage+1):length(weight.at.age)])
      }
  }  
  
  if(biomass) {
    output = list(lingcod_tot, rockfish_tot, fishery_yield_sims_N, fishery_yield_sims_B, 
                  lingcod_tot_B, rockfish_tot_B)
    names(output) = c("lingcod_pop", "rockfish_pop", "fishery_yields_N", "fishery_yields_B","lingcod_pop_B", "rockfish_pop_B")
    list2env(output, envir = .GlobalEnv) 
    } else {
      output = list(lingcod_tot, rockfish_tot, fishery_yield_sims_N, fishery_yield_sims_B)
      names(output) = c("lingcod_pop", "rockfish_pop", "fishery_yields_N", "fishery_yields_B")
      list2env(output, envir = .GlobalEnv) }
}




