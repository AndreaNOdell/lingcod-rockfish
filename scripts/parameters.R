# Parameters
library(tidyverse)
source("scripts/functions.R") #Load in functions

# Rockfish ----------------------------------------------------------------------
# From 2017 stock assessment
age = 1:65 # age classes

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
mat.at.age = c(rep(0, 14), 0.1, 0.2, 0.27, 0.33, 0.4, 0.47, 0.53, 0.6,  0.67, 0.73, 0.8, 0.87, 
               0.9, rep(1, 38))


rockfish = list(length.at.age = length.at.age, # vector of length at age
                weight.at.age = weight.at.age, # vector of weight at age
                mat.at.age = mat.at.age, # vector of maturity at age
                nat.mort = 0.044, # Natural mortality estimate
                h = 0.718, # Steepness
                age = age, # age classes
                nage = length(age), # of age classes
                r0 = 220,
                selectivity = c(rep(0, 3), rep(1, 62))) # recruitment at unfished biomass
save(rockfish, file = "cleaned_data/rockfish_parms.Rdata")


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
  # length.at.age = matrix(NA, nrow = length(Linf), ncol = length(age))
  # rownames(length.at.age) = c("female", "male")
  # for(i in c("female","male")) {
  #   length.at.age[i,] = calc_lengthxage(Linf[i], k[i], age)
  # }
length.at.age = c(calc_lengthxage(Linf["female"], k["female"], age), calc_lengthxage(Linf["male"], k["male"], age))
names(length.at.age) = c(paste0("LF_", age), paste0("LM_", age))

# Calculate length at age into matrix
weight.at.age = calc_weightxlength(length.at.age, a, b) # weight in kg

# My own estimated vector of maturity for each age class and sex
    #mat.at.age = rbind(c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)),
    #                   c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)))
    #rownames(mat.at.age) = c("female", "male")
mat.at.age <- c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14), 0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14))
names(mat.at.age) = c(paste0("LF_", age), paste0("LM_", age))
    
lingcod = list(length.at.age = length.at.age, # vector of length at age
               weight.at.age = weight.at.age, # vector of weight at age
               mat.at.age = mat.at.age, # vector of maturity at age
               nat.mort = c(0.18, 0.32), # Natural mortality estimate
               h = 0.8,# Steepness
               age = age, #age classes
               nage = length(age), # of age classes
               r0 = 4848,  # recruitment at unfished biomass
               selectivity = c(rep(0, 4), rep(1, 16)))
names(lingcod$nat.mort) = c("female", "male") 
save(lingcod, file = "cleaned_data/lingcod_parms.Rdata")

rm(mat.at.age, length.at.age, weight.at.age, age, a, b, k, Linf) #clean up environment










