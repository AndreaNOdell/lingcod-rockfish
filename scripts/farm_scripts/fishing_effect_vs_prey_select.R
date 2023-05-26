library("SimDesign", lib.loc = "../R_Packages/R4.2.2")
library("purrr", lib.loc = "../R_Packages/R4.2.2")
library("deSolve", lib.loc = "../R_Packages/R4.2.2")
library(parallel, lib.loc = "../R_Packages/R4.2.2")
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/functions.R")

# model ----------

# numerical estimation of continuous model between time step t and t+1
# The output will be the number of individuals just before the next time step
predprey_int = function(t, n, parms) {
  parms$fish.mort <- parms$fish.mort[parms$time]
  parms$bycatch <- parms$bycatch[,parms$time]
  with(as.list(parms), { # extract parameters from parms vector
    dn = rep(0, length(n)) # initialize dn/dt vector
    dn = -(M + bycatch*selectivity*c(rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40), rep(fish.mort,rockfish$nage))) * n - c(rep(0, 40), rowSums(t(t((-log(1 - pmin((consump_yelloweye_n / n[41:105]), 0.99999)) * n[41:105]) %*% diag(n[1:40])) / (1 + handling * colSums(-log(1 - pmin((consump_yelloweye_n / n[41:105]), 0.99999)) * n[41:105]) + handling*otherprey_n)))) # continuous dynamics
    dc = rowSums(t(t((-log(1 - pmin((consump_yelloweye_n / n[41:105]), 0.99999)) * n[41:105]) %*% diag(n[1:40])) / (1 + handling * colSums(-log(1 - pmin((consump_yelloweye_n / n[41:105]), 0.99999)) * n[41:105]) + handling*otherprey_n))) * rockfish$weight.at.age  # calculate predation in biomass
    df = sum(bycatch[1:40]*selectivity[1:40]*rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40) * n[1:40]) # calculate fishery yield in number
    dfb = sum((bycatch[1:40]*selectivity[1:40]*rep((fish.mort*sum(n[c(5:20,25:40)]))/sum(lingcod_harvest*n[1:40]), 40)) * (n[1:40]*weight[1:40])) # calculate fishery yield in biomass
    return(list(dn, dc, df, dfb)) # return dn/dt and df/dt as a list
  })
}

get_pop_ts = function(rockfish, lingcod, nsim = 100, corr = 0, autocorr_L = 0.23, autocorr_R = 0.23, cv = 0.5,
                      tf = (5+20+300), rec.yr = 20, mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                      f, b = 0.1, quant95 = 0.29, handling, min.selectivity = TRUE,
                      min.age.consumed = 4, rockfish.prop = 0.03, yelloweye.prop = 0.05) {
  
  # Create empty matrix for data
  nmat_sims = array(NA, dim = c(2*lingcod$nage+rockfish$nage, tf, nsim), 
                    dimnames=list(spc.age=NULL, 
                                  year=NULL, 
                                  simulation=NULL))
  rockfish_biomass_ts = matrix(NA, nrow = nsim, ncol = tf)
  rockfish_age_str_ts = matrix(NA, nrow = nsim, ncol = tf)
  lingcod_biomass_ts = matrix(NA, nrow = nsim, ncol = tf)
  early_recruit = numeric(nsim) # empty vector
  
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
                                         with(rockfish, rep(b,nage)))), (tf-5)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch 
  pre.mpa.bycatch = matrix(rep(unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                                        rep(1,nage))), # Lingcod male bycatch
                                        with(rockfish, rep(hist.by, nage)))), (5)), nrow = 2*lingcod$nage+rockfish$nage) # rockfish bycatch
  bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch)
  selectivity = c(lingcod_harvest, # Lingcod fishing selectivity
                  rockfish$selectivity) # rockfish fishing selectivity
  fish.mort = c(rep(hist.f, 5), rep(0, mpa.yr), rep(f, (tf-(5+mpa.yr))))
  # diet 
  min.age.consumed = min.age.consumed
  rockfish.prop = rockfish.prop
  yelloweye.prop = yelloweye.prop
  diet.frac.yelloweye <- rep(c(seq(0, rockfish.prop, length = min.age.consumed)*yelloweye.prop, rep(rockfish.prop*yelloweye.prop, (lingcod$nage-4))), 2)
  a_ij = get_binnedsizespec(quant95)
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
  times = 1:2
  
  
  det.burn.in = burn.in(lingcod = lingcod, rockfish = rockfish, phi.l = phi.l, phi.r = phi.r, weight.at.age = weight.at.age, M = M,
                        mat.at.age = mat.at.age, lingcod_harvest = lingcod_harvest, selectivity = selectivity, handling = handling,
                        a_ij = a_ij, consump_yelloweye_n = consump_yelloweye_n, otherprey_n = otherprey_n, hist.f = hist.f)
  
  for(i in 1:nsim) {
    
    # filling nmat_sims in
    nmat_sims[,,i] = NA
    n0 = n = det.burn.in # vector of initial abundances  
    nmat_sims[,1,i] = n0
    mn = 1
    cv_rec = sqrt(log(cv^2+1))
    list.of.seeds = 1:nsim
    set.seed(list.of.seeds[i])
    autocorr = c(autocorr_L, autocorr_R)
    corr_eps = sim_correlated_ar_ts(corr, autocorr, cv = cv_rec, mn = log(mn) - cv_rec^2/2, tf, npops = 2, ind_pops = NULL) # each simulation needs a new stochastic recrtuitment time series
    
    # calculate average recruitment in earlier years
    early_recruit[i] = mean(corr_eps[1:rec.yr,2]) # take the average of the early subset of recruitment for rockfish for each simulation
    
    for(t in 2:tf) {
      # Calculate recruitment based on previous time step
      SB = nmat_sims[,t-1, i] * weight.at.age * mat.at.age # calculate spawning biomass from prev. year
      nmat_sims[1,t, i] = nmat_sims[with(lingcod, (nage+1)), t, i] = with(lingcod, (BevHolt(phi.l, h, r0, sum(SB[1:nage]))*corr_eps[t,1])/2) # Lingcod Recruitment
      nmat_sims[with(lingcod, (2*nage+1)), t, i] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat_sims)]))*corr_eps[t,2]) # Rockfish Recruitment
      
      # Numerically integrate abundance after one time step
      parms$time = id[t-1] # assign time step ID to select time-varying parameter for lsoda
      out = as.matrix(lsoda(n0, times, predprey_int, parms)) # run lsoda
      parms$fish.mort = c(rep(hist.f, 5), rep(0, mpa.yr), rep(f, (tf-(5+mpa.yr)))) # re-enter the vector of fish.mort
      parms$bycatch = cbind(pre.mpa.bycatch, post.mpa.bycatch) # re-enter the vector of bycatch
      
      # now move age classes up one step
      nmat_sims[2:with(lingcod, nage-1), t, i] = out[2, 2:with(lingcod, nage-1)] # Female lingcod
      nmat_sims[with(lingcod, nage), t, i] = out[2, with(lingcod, nage)] + out[2, with(lingcod, (nage+1))] # Female plus group
      nmat_sims[with(lingcod, (nage+2):(2*nage-1)), t, i] = out[2, with(lingcod, (nage+2):(2*nage-1))] # Male lingcod age increase
      nmat_sims[with(lingcod, (2*nage)), t, i] = out[2, with(lingcod, (2*nage))] + out[2, with(lingcod, (2*nage+1))] # Male plus group
      nmat_sims[with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1), t, i] = out[2, with(lingcod, (2*nage+2)):(nrow(nmat_sims)-1)] # Rockfish age increase
      nmat_sims[nrow(nmat_sims), t, i] = out[2, (nrow(nmat_sims))] + out[2,(nrow(nmat_sims)+1)] #rockfish plus group
      
      # Save lingcod consumption information from ODE
      #lingcod_consumption_ts[i, (t-1)] = out[1, (ncol(out)-2)]
      
      # Now set the vector of abundances as the new n0 for next lsoda run
      n0 = nmat_sims[,t,i]
    }
    
    # Collecting and saving the necessary information
    rockfish_biomass_ts[i,] = colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i]*rockfish$weight.at.age * rockfish$mat.at.age) # spawning biomass
    rockfish_age_str_ts[i,] = nmat_sims[nrow(nmat_sims),,i] / colSums(nmat_sims[with(lingcod, (2*nage+1)):nrow(nmat_sims),,i])
    lingcod_biomass_ts[i,] = colSums(nmat_sims[1:with(lingcod, (2*nage)),,i]*lingcod$weight.at.age) # spawning biomass
  }
  
  # Replace NAs with 0s in the rockfish biomass timeseries
  rockfish_biomass_ts[is.na(rockfish_biomass_ts)] <- 0
  rockfish_age_str_ts[is.na(rockfish_age_str_ts)] <- 0
  
  # Calculate cv of spawning biomass timeseries for each simulation and take the mean of cv (i.e average cv across simulation)
  rockfish_avg = mean(apply(rockfish_biomass_ts[,-(1:150)], 1, mean))
  rockfish_cv = mean(apply(rockfish_biomass_ts[,-(1:150)], 1, function(x) sd(x) / mean(x) * 100) %>% 
                       replace(is.na(.), 0))
  rockfish_age = mean(apply(rockfish_age_str_ts[,-(1:150)], 1, mean))
  
  # Send output to global environment
  output = list(rockfish_biomass_ts[,6:tf], lingcod_biomass_ts[,6:tf], rockfish_cv, rockfish_avg, rockfish_age, early_recruit)
  names(output) = c("rockfish_biomass_ts", "lingcod_biomass_ts", "rockfish_cv", "rockfish_avg", "rockfish_age", "early_recruit")
  return(output) # sends as list
}


##  Select scenarios --------------------------------------------
SBt = 0.4 # Specified SBt/SB0eq value of interest. In most cases, 0.4
nsims = 10 # number of simulations
cv = 0.5
autocorr_R = 0.23
mpa.year = 0
rec.year = 25
handl = 0.3
bycatch = 0.05
rock.prop = 0.05 # base prey selectivity
yellow.prop = 0.02 #base prey selectivity


fishing = seq(0, 0.4, length.out = 3)

fishing_base_selec <- mcmapply(get_pop_ts, f = fishing,
                                    MoreArgs = list(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                                    cv = cv, tf = (5+mpa.year+150), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, b = bycatch, 
                                                    quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                                                    rockfish.prop = rock.prop, yelloweye.prop = yellow.prop), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)

save(fishing_base_selec, file = "fishing_data/fishing_base_selec.Rdata")
save(fishing, file = "fishing_data/fishing.Rdata")



rock.prop = 0.1 # extreme prey selectivity
yellow.prop = 1 # extreme prey selectivity

fishing_extreme_selec <- mcmapply(get_pop_ts, f = fishing,
                               MoreArgs = list(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                               cv = cv, tf = (5+mpa.year+350), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, b = bycatch, 
                                               quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                                               rockfish.prop = rock.prop, yelloweye.prop = yellow.prop), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)

save(fishing_base_selec, file = "fishing_data/fishing_extreme_selec.Rdata")

