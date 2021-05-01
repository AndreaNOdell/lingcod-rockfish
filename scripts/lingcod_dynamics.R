##############
# Parameters #
##############
age = 1:20 # vector of ages
nage = length(age) # number of ages
tf = 100 # final time
n0 = rep(100, nage) # Initial starting size
nsex = 2 # number of sex 

# Weight calculation
Linf = c(100.9, 86.3) # Linf values for lingcod (south) from stock assessment
k = c(0.191, 0.214) # growth rate coefficient for lingcod (south) from stock assessment
names(Linf) = names(k) = c("female", "male")
aL = c(0.000003308, 0.000002179) # allometric scaling parameter
bL = c(3.248, 3.36) # allometric scaling parameter
names(aL) = names(bL) = c("female", "male")

# phi calculation


# Beverton Holt
r0 = 4848 # Recruitment at unfished biomass
h = 0.8
alpha = 4*h / (phi*(1-h))
beta = (5*h-1) / (phi*r0*(1-h))


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

BevHolt = function(SBL) (alpha*SBL) / (beta + SBL)

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

#########
# model #
#########

# Create empty matrix to fill in individuals per age (row) through time (column)
nmat = matrix(NA, nrow = nage, ncol = tf)
nmat[,1] = n0


for(t in 2:tf) {
# At new time step, first calculate recruitment via spawning biomass and input into first row
  SBLs = numeric(nsex) # create empty SBL vector
  names(SBLs) = c("female", "male") # name the columns
    for(ling.sex in c("female","male")) { # Calculate Spawning Biomass for each sex
      SBLs[ling.sex] = sum(0.5*nmat[,t-1] * wl[ling.sex] * ul) 
    }
  SBL = sum(SBLs) # sum of male and female spawning biomass for total spawning biomass
  nmat[1,t] = 0.5*BevHolt(SBL) # input Bev Holt recruitment into first row of time t

# Then calculate number of individuals in subsequent ages

    nmat[2:nage, t] = nmat[1:(nage-1), t-1] * exp(-M)
    
    nmat[nage, t] = (nmat[nage-1, t-1] * exp(-M)) / (1-exp(-M))
}

# make nmat into an array - one dimension be sex



