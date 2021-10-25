
# Set up -------------------------------------------------------------------------
library(deSolve)
library(tidyverse)
source("scripts/functions.R")
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")

# numerical estimation of continuous model between time step t and t+1
# The output will be the number of individuals just before the next time step
predprey_int = function(t, n, parms) {
  with(as.list(parms), { # extract parameters from parms vector
    dn = rep(0, length(n)) # initialize dn/dt vector
    dn = ( -(M + bycatch*selectivity*fish.mort) * n ) - c(rep(0, 40), rowSums( (t(t(a_ij) %*% diag(n[41:105])) %*% diag(n[1:40])) / ((1 + handling * colSums(t(t(a_ij) %*% diag(n[41:105]))))[col(a_ij)]))) # continuous dynamics
    return(list(dn)) # return dn/dt as a list
  })
}


# model ------------------------------------------------------------------------

# Parameters for recruitment 
  # Calculate phi and create vector
    phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))
    phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))
  # Create weight and maturity at age vectors
    weight.at.age = c(lingcod$weight.at.age, rockfish$weight.at.age)
    mat.at.age = c(lingcod$mat.at.age, rockfish$mat.at.age)
    
# Consumption
load("cleaned_data/binned.size.spec.Rdata") # Load in rockfish size-spectra in size-specific lingod diet
# binned.size.spec[65,] <- 0
diet.frac.rockfish <- c(rep(0, 4), rep(0.03, (lingcod$nage-4)), rep(0, 4), rep(0.03, (lingcod$nage-4)))
handling <- 0.01
a_ij = binned.size.spec %*% diag(diet.frac.rockfish)
# rowSums((a_ij * n0[41:105] * n0[1:40])/(1 + handling * rowSums(a_ij) * n0[41:105]))

consumption <- c(rep(0, 40), rowSums((t(t(a_ij) %*% diag(n[41:105])) %*% diag(n[1:40])) / ((1 + handling * colSums(t(t(a_ij) %*% diag(n[41:105]))))[col(a_ij)])))


# Parameters for lsoda
fish.mort = 0.1
M = unname(c(with(lingcod, c(rep(nat.mort["female"],nage), # Lingcod female natural mortality
                             rep(nat.mort["male"],nage))), # Lingcod male natural mortality
             with(rockfish, rep(nat.mort,nage))))         # Rockfish natural mortality
bycatch = unname(c(with(lingcod, c(rep(1,nage), # Lingcod female bycatch
                                   rep(1,nage))), # Lingcod male bycatch
                   with(rockfish, rep(0.1,nage)))) # rockfish bycatch
selectivity = c(with(lingcod, rep(selectivity, length(nat.mort))), # Lingcod female fishing selectivity
                with(rockfish, selectivity)) # rockfish fishing selectivity
parms = list(fish.mort = fish.mort, M = M, bycatch = bycatch, selectivity = selectivity, handling = handling, a_ij = a_ij)



# Now for model-building
times = 1:2 # we only want to run the lsoda for one time step because we're going to pulse recruitment and age
n_name <- c(with(lingcod,c(paste0("LF_", age), paste0("LM_", age))),
           with(rockfish, paste0("R_", age))) # add names for the matrix to signify what is lingcod and rockfish
n0 = c(n = rep(10, length(M))) # vector of initial abundances 
names(n0) = n_name # add names to the initial population sizes

# Create empty matrix to fill in individuals per age (row) through time (column)
# The rows of our matrix will have LF_1-LF20, LM1-LM20, R_1-R_65. columns will be each time step
tf = 200 # run for 100 time steps
nmat = matrix(NA, nrow = length(n0), ncol = tf) # generate empty matrix to fill in with projections
rownames(nmat) = n_name # add names
nmat[,1] = n0 # add initial abundance to time step 1


for(t in 2:tf) {
  
  # Calculate recruitment based on previous time step
    # Caclulate spawning biomass
      SB = nmat[,t-1] * weight.at.age * mat.at.age
      nmat[1,t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Female Recruitment
      nmat[with(lingcod, (nage+1)), t] = with(lingcod, (BevHolt(phi.l[1], h, r0, sum((SB[1:nage]))))/2) # Lingcod Male Recruitment
      nmat[with(lingcod, (2*nage+1)), t] = with(rockfish, BevHolt(phi.r, h, r0, sum(SB[with(lingcod, (2*nage+1)):nrow(nmat)]))) # Rockfish Recruitment
  
  # Numerically integrate abundance after one time step
    out = as.matrix(lsoda(n0, times, predprey_int, parms))
  
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


lingcod_tot = colSums(nmat[1:with(lingcod, 2*nage),])
rockfish_tot = colSums(nmat[with(lingcod, (2*nage+1)):nrow(nmat),])

plot(1:tf, lingcod_tot, type = "l")
plot(1:tf, rockfish_tot, type = "l")








