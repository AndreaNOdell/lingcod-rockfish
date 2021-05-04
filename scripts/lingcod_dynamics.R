#############
# Functions #
#############

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

##############
# Parameters #
##############
age = 1:20 # vector of ages
nage = length(age) # number of ages
tf = 100 # final time
n0 = rep(100, nage) # Initial starting size
nsex = 2 # number of sex 
M = c(0.18, 0.32) # Natural Mortality estimate
names(M) = c("female", "male") 

# Weight calculation
Linf = c(100.9, 86.3) # Linf values for lingcod (south) from stock assessment
k = c(0.191, 0.214) # growth rate coefficient for lingcod (south) from stock assessment
names(Linf) = names(k) = c("female", "male")
aL = c(0.000003308, 0.000002179) # allometric scaling parameter
bL = c(3.248, 3.36) # allometric scaling parameter
names(aL) = names(bL) = c("female", "male")

# Maturity (arbitrary)
L = 0:100
length_age_maturity = matrix(NA, nrow = 3, ncol = length(L))
length_age_maturity[1,] = L
length_age_maturity[2,] = floor((log(1-(L/100.9)))/(-0.191))
length_age_maturity[3,1:29] = 0
length_age_maturity[3,30:40] = seq(0, 0.08, length.out = 11)
length_age_maturity[3,41:53] = seq(0.09, .52, length.out = 13)
length_age_maturity[3,54:60] = seq(0.56, .85, length.out = 7)
length_age_maturity[3,61:65] = seq(0.87, .95, length.out = 5)
length_age_maturity[3,66:70] = seq(0.96, 1, length.out = 5)
length_age_maturity[3,71:101] = 1

# My own estimate of maturity for each age class
ul = rbind(c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)),
          c(0, 0, 0.1, 0.4, 0.75, 0.97, rep(1, 14)))
rownames(ul) = c("female", "male")
  
##########################
# Calculate weight - age #
##########################

# Calculate length at age into matrix
length_l = matrix(NA, nrow = length(Linf), ncol = nage)
rownames(length_l) = c("female", "male")
for(i in c("female","male")) {
  length_l[i,] = calc_lengthxage(Linf[i], k[i], age)
}

# Calculate length at age into matrix
wl = calc_weightxlength(length_l, aL, bL) # weight in kg


#################
# Calculate phi #
#################

# First calculate survival into each age class (L(a) * exp(-M)) for each sex
phi_n0 = 1
phi_struct = matrix(NA, nrow = 2, ncol = nage)
phi_struct[,1] = phi_n0
rownames(phi_struct) = c("female", "male")
for(ling.sex in c("female", "male")) { # Calculate survival across ages
  for(j in 1:(nage-1)) {
    phi_struct[ling.sex,j+1] = phi_struct[ling.sex,j] * exp(-M[ling.sex])
  }
}

# Then multiply each age class by weights and maturity and inverse
sb_r = phi_struct * wl * ul # Spawners per recruit
phi = 1/rowSums(sb_r) # Recruits per spawner

# Beverton Holt parameters
r0 = 4848 # Recruitment at unfished biomass
h = 0.8
alpha = 4*h / (phi*(1-h)) # Carrying Capacity
beta = (5*h-1) / (phi*r0*(1-h)) # Steepness



#########
# model #
#########

# Create empty matrix to fill in individuals per age (row) through time (column)
nmat = array(NA, dim = c(nage, tf, 2), 
             dimnames=list(age=NULL, 
                           year=NULL, ling.sex=c('female', 'male')))
nmat[,1,] = n0

for(ling.sex in c("female", "male")){
  for(t in 2:tf) {
# At new time step, first calculate recruitment via spawning biomass and input into first row
    SBLs = numeric(nsex) # create empty SBL vector
    names(SBLs) = c("female", "male") # name the columns
    for(lingcod.sex in c("female", "male")){ # Calculate Spawning Biomass for each sex
      SBLs[lingcod.sex] = sum((0.5*nmat[,t-1,lingcod.sex]) * wl[lingcod.sex,] * ul[lingcod.sex,]) 
    }
  SBL = sum(SBLs) # sum of male and female spawning biomass for total spawning biomass
  nmat[1,t,ling.sex] = 0.5*(BevHolt(SBL)[ling.sex]) # input Bev Holt recruitment into first row of time t

# Then calculate number of individuals in subsequent ages
  nmat[2:nage, t, ling.sex] = nmat[1:(nage-1), t-1, ling.sex] * exp(-M[ling.sex])
  nmat[nage, t, ling.sex] = (nmat[nage-1, t-1, ling.sex] * exp(-M[ling.sex])) / (1-exp(-M[ling.sex]))
  }
}

# Something is wrong with the spawning biomass input



