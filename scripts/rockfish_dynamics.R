
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

BevHolt = function(alpha, beta, SBL) (alpha*SBL) / (1 + beta*SBL)

# Parameters -------------------------------------------------------------------

age = 0:65 # vector of ages
nage = length(age) # number of ages
tf = 100 # final time
n0 = rep(50, nage) # Initial starting size
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
r_length_age_maturity[2,64:length(L)] = 86 # fill in last ages
r_length_age_maturity[3,1:34] = 0 # immature up until age 33
r_length_age_maturity[3,35:50] = seq(0, 1, length.out = 16) 
r_length_age_maturity[3,50:66] = 1

# My own estimated vector of maturity for each age class
ul = c(rep(0, 15), 0.1, 0.2, 0.27, 0.33, 0.4, 0.47, 0.53, 0.6,  0.67, 0.73, 0.8, 0.87, 
       0.9, rep(1, 38))

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

# Then multiply each age class by weights and maturity and inverse
sb_r = phi_struct * wl * ul # Spawners per recruit
phi = 1/sum(sb_r) # Recruits per spawner

# Recruitment parameters
r0 = 220 # Recruitment at unfished biomass
h = 0.718
alpha = 4*h / (phi*(1-h)) # Carrying Capacity
beta = (5*h-1) / (phi*r0*(1-h)) # Steepness
r_sd = 0.6 # standard deviation for lognormal distribution for stochastic recruitment

# Fishing ----------------------------------------------------------------------

f = 0.1 # fishing mortality
vl = c(rep(0, 4), rep(1, 16)) # susceptibility to fishery (0 or 1)
b = 0.2
mort = M + b*f

# model ------------------------------------------------------------------------

# Create empty matrix to fill in individuals per age (row) through time (column)
nmat = matrix(NA, nrow = nage, ncol = tf)
nmat[,1] = n0

for(t in 2:tf) {
  eps_r = rlnorm(1, meanlog = -0.5*r_sd^2, sdlog = r_sd) # lognormal distribution for varying r
  SBL = sum((nmat[,t-1]) * wl * ul) 
  nmat[1,t] = (BevHolt(alpha, beta, SBL))*eps_r # input Bev Holt recruitment into first row of time t
    
  # Then calculate number of individuals in subsequent ages
  nmat[2:nage, t] = nmat[1:(nage-1), t-1] * exp(-mort)
  nmat[nage, t] = (nmat[nage-1, t-1] * exp(-mort)) + (nmat[nage, t-1] * exp(-mort))
}

# ggplot -----------------------------------------------------------------------

# first let's set up our data into a long format
ntot = colSums(nmat)
ntot = as.data.frame(ntot) %>% 
  mutate(time = 1:tf) %>% 
  rename(abundance = ntot)
# now let's graph!
ggplot(ntot, aes(x = time, y = abundance)) +
  geom_line()











