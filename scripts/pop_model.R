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
phi_calc = function(age, nat.mort, weight.at.age, mat.at.age) {
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
  print(phi)
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
                h = 0.8, # Steepness
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
                nage = length(age), # of age classes
                r0 = 220) # recruitment at unfished biomass












