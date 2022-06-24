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

# Determine selectivity - 25cm from taylor and wetzel
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
                selectivity = selectivity,
                phi.r = with(rockfish, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = FALSE))) 
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

# Calculate length at age in cm
length.at.age = c(calc_lengthxage(Linf["female"], k["female"], age), calc_lengthxage(Linf["male"], k["male"], age))
names(length.at.age) = c(paste0("LF_", age), paste0("LM_", age))
# Calculate weight at age into matrix in kg
weight.at.age = c(calc_weightxlength(length.at.age[1:20], a["female"], b["female"]), calc_weightxlength(length.at.age[21:40], a["male"], b["male"])) # weight in kg

# Calculate maturity at age
alpha <- c(.994, 1.06); beta <- c(4.323, 2.506)
names(alpha) <- names(beta) <- c('female', 'male')
mat.at.age <- sapply(age, function(x) 1/(1+exp(-alpha*(x-beta))))
mat.at.age = c(mat.at.age[1,], mat.at.age[2,])
names(mat.at.age) = c(paste0("LF_", age), paste0("LM_", age))

# Determine selectivity - current minimum size limit = 22" i.e. 56cm
min.vulnerable.length = 56
max.allow.length = 85
# Minimum size limit
min.selectivity = length.at.age
min.selectivity[min.selectivity<min.vulnerable.length] = 0
min.selectivity[min.selectivity>min.vulnerable.length] = 1
# harvest slot
harvest.slot = length.at.age
harvest.slot[harvest.slot<min.vulnerable.length] = 0 ; harvest.slot[harvest.slot>max.allow.length] = 0
harvest.slot[harvest.slot>min.vulnerable.length & harvest.slot<max.allow.length] = 1
# female biased (sex ratio 38% male deeper than 200ft from Fish Bulletin summary of Lingcod Life History)
#female.biased = length.at.age

# Lingcod consumption rate ------------------------------------------------

# From Beaudreau & Essington 2009
consump.a <- c(rep(10, lingcod$nage), rep(11, lingcod$nage))
consump.b <- c(rep(.75, lingcod$nage), rep(.77, lingcod$nage))
q.10 <- exp(10*0.0647) # From Beaudreau 2009
# Convert daily consumption in g to annual consumption in kg, on per capita basis
consump.at.age.9C <- (365*consump.a*(lingcod$weight.at.age)^consump.b/1000)
# Consumption rate calculated in Puget Sound, ~9C. RCAs are ~6C
consump.at.age <- consump.at.age.9C/(q.10^(3/10))
norm.consump.at.age <- c(consump.at.age[1:20]/sum(consump.at.age[1:20]), consump.at.age[21:40]/sum(consump.at.age[21:40]))
# actual.consump.a <- consump.a *365/1000/(q.10^.3)

    
lingcod = list(length.at.age = length.at.age, # vector of length at age
               weight.at.age = weight.at.age, # vector of weight at age
               mat.at.age = mat.at.age, # vector of maturity at age
               nat.mort = c(0.18, 0.32), # Natural mortality estimate
               h = 0.8,# Steepness
               age = age, #age classes
               nage = length(age), # of age classes
               r0 = 4848,  # recruitment at unfished biomass
               min.selectivity = min.selectivity,
               harvest.slot = harvest.slot,
               consump.at.age = consump.at.age,
               norm.consump.at.age = norm.consump.at.age,
               phi.l = with(lingcod, phi_calc(age, nat.mort, weight.at.age, mat.at.age, lingcod = TRUE))[1])
names(lingcod$nat.mort) = c("female", "male") 
save(lingcod, file = "cleaned_data/lingcod_parms.Rdata")

rm(mat.at.age, length.at.age, weight.at.age, age, a, b, k, Linf, consump.a, consump.at.age, consump.b, consump.at.age.9C, norm.consump.at.age) #clean up environment










