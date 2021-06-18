# Set up ----------------------------------------------------------------------

library(tidyverse)
library(PNWColors)
theme_set(theme_classic()) # set all ggplot graphs to classic theme
source("scripts/functions.R") #Load in functions


# Lingcod Parameters ----------------------------------------------------------
# From 2017 stock assessment
age = 1:20 # age classes
# Weight calculation
Linf = c(100.9, 86.3) # Linf values for lingcod (south) from stock assessment
k = c(0.191, 0.214) # growth rate coefficient for lingcod (south) from stock assessment
a = c(3.308*10^-6, 2.179*10^-6) # allometric scaling parameter
b = c(3.248, 3.36) # allometric scaling parameter
names(Linf) = names(k) = names(a) = names(b) = c("female", "male")

# Calculate length at age into matrix
length.at.age = matrix(NA, nrow = length(Linf), ncol = length(age))
rownames(length.at.age) = c("female", "male")
for(i in c("female","male")) {
  length.at.age[i,] = calc_lengthxage(Linf[i], k[i], age)
}
# Calculate length at age into matrix
weight.at.age = calc_weightxlength(length.at.age, a, b) # weight in kg
# My own estimated vector of maturity for each age class and sex
mat.at.age = rbind(c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)),
           c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)))
rownames(mat.at.age) = c("female", "male")

lingcod = list(length.at.age = length.at.age, # vector of length at age
                weight.at.age = weight.at.age, # vector of weight at age
                mat.at.age = mat.at.age, # vector of maturity at age
                nat.mort = c(0.18, 0.32), # Natural mortality estimate
                h = 0.8,# Steepness
                age = age, #age classes
                nage = length(age), # of age classes
                r0 = 4848,  # recruitment at unfished biomass
                susceptibility = rbind(c(rep(0, 4), rep(1, 16)),c(rep(0, 4), rep(1, 16))))
rownames(lingcod$susceptibility) = names(lingcod$nat.mort) = c("female", "male") 


# Rockfish Parameters ----------------------------------------------------------
# From 2017 stock assessment
age = 0:65 # age classes

# Weight calculation
Linf = 63.9 # Linf value for rockfish from Kiva's paper
k = 0.049 # growth rate coefficient for rockfish from stock assessment
a = 7.312807*10^-6 # allometric scaling parameter
b = 3.242482 # allometric scaling parameter

# Calculate length at age into matrix
length.at.age = calc_lengthxage(Linf, k, age)
# Calculate weight at age using length.at.age vector
weight.at.age = calc_weightxlength(length.at.age, a, b) # weight in kg
# My own estimated vector of maturity for each age class
mat.at.age = c(rep(0, 15), 0.1, 0.2, 0.27, 0.33, 0.4, 0.47, 0.53, 0.6,  0.67, 0.73, 0.8, 0.87, 
       0.9, rep(1, 38))


rockfish = list(length.at.age = length.at.age, # vector of length at age
                weight.at.age = weight.at.age, # vector of weight at age
                mat.at.age = mat.at.age, # vector of maturity at age
                nat.mort = 0.044, # Natural mortality estimate
                h = 0.718, # Steepness
                age = age, # age classes
                nage = length(age), # of age classes
                r0 = 220,
                selectivity = c(rep(0, 4), rep(1, 62))) # recruitment at unfished biomass

rm(mat.at.age, length.at.age, weight.at.age, age, a, b, k, Linf) #clean up environment


# model ------------------------------------------------------------------------

get_popn = function(rockfish, lingcod, weight.at.age, mat.at.age, age, selectivity, nat.mort, fish.mort = 0.1, bycatch = 0, r_sd = 0.6, h, r0, nage, tf = 100, n.init.l = 100, n.init.r = 50, nsims = 10) {
  # Calculate phi for bev holt curve
  phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))
  phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))
  # empty matrix for age-structured population across years
  n0.l = rep(n.init.l, length(lingcod$age)) # Initial starting size
  nmat.l = with(lingcod, array(NA, dim = c(length(age), tf, 2), 
                 dimnames=list(age=NULL, 
                               year=NULL, 
                               ling.sex=c('female', 'male'))))
  nmat.l[,1,] = n0.l 
  
  n0.r = rep(n.init.r, length(rockfish$age)) # Initial starting size
  nmat.r = with(rockfish, matrix(NA, nrow = length(age), ncol = tf))
  nmat.r[,1] = n0.r

# empty array/matrix for total population size across years 
 ntot_l = array(NA, dim = c(nsims, tf, 2), # empty array for total lingcod population size across years
                            dimnames=list(simulations=1:nsims, 
                                          year=1:tf, 
                                          ling.sex=c('female', 'male')))
 ntot_r = matrix(NA, nrow = nsims, ncol = tf) # empty matrix for total rockfish population size across years

# Simulations
  for(n in 1:nsims) {
    
# lingcod 
  # run for loop
    for(t in 2:tf) { # loop through time
      eps_r = rlnorm(1, meanlog = 0.5*r_sd^2, sdlog = r_sd)# lognormal distribution for varying r - same value for both popn
      for(ling.sex in c("female", "male")){
        M = with(lingcod, nat.mort[ling.sex] + fish.mort*susceptibility[ling.sex,])
        # At new time step, first calculate recruitment via spawning biomass and input into first row
        SBLs = numeric(2) # create empty SBL vector
        names(SBLs) = c("female", "male") # name the columns
        for(lingcod.sex in c("female", "male")){ # Calculate Spawning Biomass for each sex
          SBLs[lingcod.sex] = with(lingcod, sum((0.5*nmat.l[,t-1,lingcod.sex]) * weight.at.age[lingcod.sex,] * mat.at.age[lingcod.sex,])) 
        }
        SBL = sum(SBLs) # sum of male and female spawning biomass for total spawning biomass
        nmat.l[1,t,ling.sex] = with(lingcod, 0.5*(BevHolt(phi.l, h, r0, SBL)[ling.sex])*eps_r) # input Bev Holt recruitment into first row of time t
        # Then calculate number of individuals in subsequent ages
        nmat.l[2:lingcod$nage, t, ling.sex] = nmat.l[1:(lingcod$nage-1), t-1, ling.sex] * exp(-M[1:(length(M)-1)])
        nmat.l[lingcod$nage, t, ling.sex] = (nmat.l[lingcod$nage-1, t-1, ling.sex] * exp(-M[length(M)])) + (nmat.l[lingcod$nage, t-1, ling.sex] * exp(-M[length(M)]))
      }
#rockfish
      M = with(rockfish, nat.mort + bycatch*fish.mort*selectivity)
      # eps_r = rlnorm(1, meanlog = -0.5*r_sd^2, sdlog = r_sd) # lognormal distribution for varying r
      SBL = with(rockfish, sum((nmat.r[,t-1]) * weight.at.age * mat.at.age))
      nmat.r[1,t] = with(rockfish, (BevHolt(phi.r, h, r0, SBL))*eps_r) # input Bev Holt recruitment into first row of time t
        
      # Then calculate number of individuals in subsequent ages
      nmat.r[2:rockfish$nage, t] = nmat.r[1:(rockfish$nage-1), t-1] * exp(-M[1:(length(M)-1)])
      nmat.r[rockfish$nage, t] = (nmat.r[rockfish$nage-1, t-1] * exp(-M[length(M)]) + (nmat.r[rockfish$nage, t-1] * exp(-M[length(M)])))
  } # end of population modeling loop
  
#lingcod total population (output)
  # first let's set up our data into a long format
  ntot_l[n,,] = colSums(nmat.l)
  #ntot_wide = as.data.frame(ntot_l)
  #ntot_long <- gather(ntot_wide, sex, abundance, female:male, factor_key=TRUE) # change to long format
  
# rockfish total population (output)
  ntot_r[n,] = colSums(nmat.r)
  #ntot_r = as.data.frame(ntot_r) %>% 
    #rename(abundance = ntot_r)
  
  } # end of simulation loop
  
# send output to the global environment
  abundance_output = list(ntot_l, ntot_r)
  names(abundance_output) = c("lingcod_pop", "rockfish_pop")
  list2env(abundance_output, envir = .GlobalEnv)
}

# Get rockfish equilibrium
get_popn(rockfish, lingcod, weight.at.age, mat.at.age, age, selectivity, nat.mort, fish.mort = 0, 
         bycatch = 0.1, r_sd = 0, h, r0, nage, tf = 200, n.init.l = 100, n.init.r = 10, nsims = 1) 
rockfish_equilibrium = rockfish_pop[,150]
rockfish_40 = 0.4*rockfish_equilibrium

# now run with 1000 simulations
get_popn(rockfish, lingcod, weight.at.age, mat.at.age, age, selectivity, nat.mort, fish.mort = 0.1, 
         bycatch = 0.1, r_sd = 0.6, h, r0, nage, tf = 100, n.init.l = 100, n.init.r = 10, nsims = 1000)


# now let's graph!
pdf("plots/rockfish_recovery.pdf")
matplot(t(rockfish_pop), type = "l", xlab = "Years", ylab = "Abundance")
abline(h = rockfish_40, lty = 2, lwd = 2)
abline(h = rockfish_equilibrium, lty = 2, lwd = 2)
dev.off()

#ggplot(ntot_r, aes(x = time, y = abundance)) +
#  geom_line() +
#  geom_line(data = ntot_long, aes(x = time, y = abundance, col = sex)) +
#  scale_color_manual(values=pnw_palette("Sunset2", 3,"discrete"))

# What fraction of simulated populations made full recovery?
fully.recovered <- logical(1000)
for(n in 1:1000) {
  fully.recovered[n] <- max(rockfish_pop[n,]) >= rockfish_equilibrium 
  }
mean(fully.recovered) # 0.317

# now increase bycatch
get_popn(rockfish, lingcod, weight.at.age, mat.at.age, age, selectivity, nat.mort, fish.mort = 0.1, 
         bycatch = 0.3, r_sd = 0.6, h, r0, nage, tf = 100, n.init.l = 100, n.init.r = 10, nsims = 1000)


# now let's graph!
pdf("plots/recovery_high_bycatch.pdf")
matplot(t(rockfish_pop), type = "l", xlab = "Years", ylab = "Abundance")
abline(h = rockfish_40, lty = 2, lwd = 2)
abline(h = rockfish_equilibrium, lty = 2, lwd = 2)
dev.off()

fully.recovered <- logical(1000)
for(n in 1:1000) {
  fully.recovered[n] <- max(rockfish_pop[n,]) >= rockfish_equilibrium 
}
mean(fully.recovered) # 0.028







