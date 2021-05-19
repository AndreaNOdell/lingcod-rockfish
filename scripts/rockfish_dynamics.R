
# Set up ----------------------------------------------------------------------

library(tidyverse)
library(PNWColors)
theme_set(theme_classic()) # set all ggplot graphs to classic theme


# Functions ----------------------------------------------------------------------

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

BevHolt = function(SBL) (alpha*SBL) / (1 + beta*SBL)

# Parameters -------------------------------------------------------------------

age = 1:20 # vector of ages
nage = length(age) # number of ages
tf = 100 # final time
n0 = rep(100, nage) # Initial starting size
M = 0.044 # Natural Mortality estimate

# Weight calculation
Linf = 63.9 # Linf value for rockfish from Kiva's paper
k = 0.049 # growth rate coefficient for rockfish from stock assessment
aL = 0.000007312807 # allometric scaling parameter
bL = 3.242482 # allometric scaling parameter

# Maturity (arbitrary)
L = 0:65
r_length_age_maturity = matrix(NA, nrow = 3, ncol = length(L))
r_length_age_maturity[1,] = L
r_length_age_maturity[2,] = floor((log(1-(L/Linf)))/(-k)) # Calculate age via length
r_length_age_maturity[2,63:length(L)] = 71 # fill in last ages
r_length_age_maturity[3,1:34] = 0 # immature up until age 33
r_length_age_maturity[3,35:50] = seq(0, 1, length.out = 16) 
r_length_age_maturity[3,50:66] = 1

# My own estimated vector of maturity for each age class
ul = c(rep(0, 15), 0.1, 0.2, 0.27, 0.33, 1)

# Calculate length and weight by age -------------------------------------------------------

# Calculate length at age into matrix
length_l = calc_lengthxage(Linf, k, age)

# Calculate weight by length
wl = calc_weightxlength(length_l, aL, bL) # weight in kg


# Calculate phi ----------------------------------------------------------------

# First calculate survival into each age class (L(a) * exp(-M)) for each sex
phi_n0 = 1 # Initial recruit
phi_struct = numeric(nage) # Empty matrix to fill
phi_struct[1] = phi_n0 # Fill in initial recruit
# Calculate survival across ages
for(j in 1:(nage-1)) {
    phi_struct[j+1] = phi_struct[j] * exp(-M)
}












