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
mat.param <- c(-0.437, 38.78)
names(mat.param) <- c('slope', 'len.at.50')
mat.at.age <- sapply(length.at.age, function(l, par) 1/(1+exp(par[1]*l-par[1]*par[2])), par=mat.param)

# Determine selectivity - 22" i.e. 56cm
min.vulnerable.length = 25
selectivity = length.at.age
selectivity[selectivity<min.vulnerable.length] = 0
selectivity[selectivity>min.vulnerable.length] = 1

rockfish = list(length.at.age = length.at.age, # vector of length at age
                weight.at.age = weight.at.age, # vector of weight at age
                mat.at.age = mat.at.age, # vector of maturity at age
                nat.mort = 0.044, # Natural mortality estimate
                h = 0.718, # Steepness
                age = age, # age classes
                nage = length(age), # of age classes
                r0 = 220, # recruitment at unfished biomass
                selectivity = selectivity) 
save(rockfish, file = "cleaned_data/rockfish_parms.Rdata")


# Lingcod Parameters ----------------------------------------------------------
# From 2017 stock assessment
age = 1:20 # age classes
# Length/weight parameters
Linf = c(100.9, 86.3) # Linf values for lingcod (south) from stock assessment
k = c(0.191, 0.214) # growth rate coefficient for lingcod (south) from stock assessment
a = c(3.308*10^-6, 2.179*10^-6) # allometric scaling parameter
b = c(3.248, 3.36) # allometric scaling parameter
names(Linf) = names(k) = names(a) = names(b) = c("female", "male")

# Calculate length at age
length.at.age = c(calc_lengthxage(Linf["female"], k["female"], age), calc_lengthxage(Linf["male"], k["male"], age))
names(length.at.age) = c(paste0("LF_", age), paste0("LM_", age))
# Calculate weight at age into matrix
weight.at.age = c(calc_weightxlength(length.at.age[1:20], a["female"], b["female"]), calc_weightxlength(length.at.age[21:40], a["male"], b["male"])) # weight in kg

# Caclulate maturity at age
alpha <- c(.994, 1.06); beta <- c(4.323, 2.506)
names(alpha) <- names(beta) <- c('female', 'male')
mat.at.age <- sapply(age, function(x) 1/(1+exp(-alpha*(x-beta))))
mat.at.age = c(mat.at.age[1,], mat.at.age[2,])
names(mat.at.age) = c(paste0("LF_", age), paste0("LM_", age))

# Determine selectivity - 22" i.e. 56cm
min.vulnerable.length = 56
selectivity = length.at.age
selectivity[selectivity<min.vulnerable.length] = 0
selectivity[selectivity>min.vulnerable.length] = 1
  
    
lingcod = list(length.at.age = length.at.age, # vector of length at age
               weight.at.age = weight.at.age, # vector of weight at age
               mat.at.age = mat.at.age, # vector of maturity at age
               nat.mort = c(0.18, 0.32), # Natural mortality estimate
               h = 0.8,# Steepness
               age = age, #age classes
               nage = length(age), # of age classes
               r0 = 4848,  # recruitment at unfished biomass
               selectivity = selectivity)
names(lingcod$nat.mort) = c("female", "male") 
save(lingcod, file = "cleaned_data/lingcod_parms.Rdata")

rm(mat.at.age, length.at.age, weight.at.age, age, a, b, k, Linf) #clean up environment










