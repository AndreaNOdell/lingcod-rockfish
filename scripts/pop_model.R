# Set up ----------------------------------------------------------------------

library(tidyverse)
library(PNWColors)
theme_set(theme_classic()) # set all ggplot graphs to classic theme


# Functions --------------------------------------------------------------------

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
                r0 = 4848) # recruitment at unfished biomass
names(lingcod$nat.mort) = c("female", "male") 


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
                r0 = 220) # recruitment at unfished biomass

rm(mat.at.age, length.at.age, weight.at.age, age, a, b, k, Linf)


# model ------------------------------------------------------------------------


get_popn = function(rockfish, lingcod, weight.at.age, mat.at.age, age, nat.mort, fish.mort = 0.1, bycatch = 0, r_sd = 0.6, h, r0, nage, tf = 100, n.init.l = 100, n.init.r = 50) {
  # Calculate phi for bev holt curve
  phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))
  phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))
  # empty matrix
  n0.l = rep(n.init.l, length(lingcod$age)) # Initial starting size
  nmat.l = with(lingcod, array(NA, dim = c(length(age), tf, 2), 
                 dimnames=list(age=NULL, 
                               year=NULL, 
                               ling.sex=c('female', 'male'))))
  nmat.l[,1,] = n0.l 
  
  n0.r = rep(n.init.r, length(rockfish$age)) # Initial starting size
  nmat.r = with(rockfish, matrix(NA, nrow = length(age), ncol = tf))
  nmat.r[,1] = n0.r
  
# lingcod 
  # run for loop
  for(t in 2:tf) { # loop through time
    eps_r = rlnorm(1, meanlog = -0.5*r_sd^2, sdlog = r_sd)# lognormal distribution for varying r - same value for both popn
    for(ling.sex in c("female", "male")){
      M = with(lingcod, nat.mort[ling.sex] + fish.mort)
      # At new time step, first calculate recruitment via spawning biomass and input into first row
      SBLs = numeric(2) # create empty SBL vector
      names(SBLs) = c("female", "male") # name the columns
      for(lingcod.sex in c("female", "male")){ # Calculate Spawning Biomass for each sex
        SBLs[lingcod.sex] = with(lingcod, sum((0.5*nmat.l[,t-1,lingcod.sex]) * weight.at.age[lingcod.sex,] * mat.at.age[lingcod.sex,])) 
      }
      SBL = sum(SBLs) # sum of male and female spawning biomass for total spawning biomass
      nmat.l[1,t,ling.sex] = with(lingcod, 0.5*(BevHolt(phi.l, h, r0, SBL)[ling.sex])*eps_r) # input Bev Holt recruitment into first row of time t
      # Then calculate number of individuals in subsequent ages
      nmat.l[2:lingcod$nage, t, ling.sex] = nmat.l[1:(lingcod$nage-1), t-1, ling.sex] * exp(-M)
      nmat.l[lingcod$nage, t, ling.sex] = (nmat.l[lingcod$nage-1, t-1, ling.sex] * exp(-M)) + (nmat.l[lingcod$nage, t-1, ling.sex] * exp(-M))
    }
#rockfish
    M = with(rockfish, nat.mort + bycatch*fish.mort)
    # eps_r = rlnorm(1, meanlog = -0.5*r_sd^2, sdlog = r_sd) # lognormal distribution for varying r
    SBL = with(rockfish, sum((nmat.r[,t-1]) * weight.at.age * mat.at.age))
    nmat.r[1,t] = with(rockfish, (BevHolt(phi.r, h, r0, SBL))*eps_r) # input Bev Holt recruitment into first row of time t
      
    # Then calculate number of individuals in subsequent ages
    nmat.r[2:rockfish$nage, t] = nmat.r[1:(rockfish$nage-1), t-1] * exp(-M)
    nmat.r[rockfish$nage, t] = (nmat.r[rockfish$nage-1, t-1] * exp(-M) + (nmat.r[rockfish$nage, t-1] * exp(-M)))
    }
  
  # send output to the global environment
  abundance_output = list(nmat.l, nmat.r)
  names(abundance_output) = c("lingcod_pop", "rockfish_pop")
  list2env(abundance_output, envir = .GlobalEnv)
}


get_popn(rockfish, lingcod, weight.at.age, mat.at.age, age, nat.mort, fish.mort = 0.1, 
         bycatch = 0.1, r_sd = 0.6, h, r0, nage, tf = 100, n.init.l = 100, n.init.r = 10)

#lingcod 
  # first let's set up our data into a long format
tf = 100
ntot_l = colSums(lingcod_pop)
ntot_wide = as.data.frame(ntot_l) %>% 
  mutate(time = 1:tf) # change from matrix to dataframe and add a column for time
ntot_long <- gather(ntot_wide, sex, abundance, female:male, factor_key=TRUE) # change to long format

# rockfish
ntot_r = colSums(rockfish_pop)
ntot_r = as.data.frame(ntot_r) %>% 
  mutate(time = 1:tf) %>% 
  rename(abundance = ntot_r)

# now let's graph!
ggplot(ntot_r, aes(x = time, y = abundance)) +
  geom_line() +
  geom_line(data = ntot_long, aes(x = time, y = abundance, col = sex)) +
  scale_color_manual(values=pnw_palette("Sunset2", 3,"discrete"))




