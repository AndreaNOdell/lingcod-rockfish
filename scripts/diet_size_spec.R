# Lingcod diet size spectrum ----------------------------------------------
# Code via Oken and Essington 2016
# yelloweye rockfish sexes combined because very little contrast in previous assessments (parsimony)

library(tidyverse)
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")

# from Beaudreau & Essington 2007
get_gamma_pars <- function(pars, quant5, quant95) {
  alpha <- exp(pars[1]); beta <- exp(pars[2])
  inv.cdf.5 <- qgamma(.05, shape=alpha, scale=beta)
  inv.cdf.95 <- qgamma(.95, shape=alpha, scale=beta)
  return((inv.cdf.5-quant5)^2+(inv.cdf.95-quant95)^2)
}

optim.result <- exp(optim(c(0,0), get_gamma_pars, quant5=.05, quant95=.29)$par)

size.spec.al <- optim.result[1] # alpha = 3.93
size.spec.be <- optim.result[2] # beta = .038 * lingcod length (in cm)

# Parameters from Echeverria & Lenarz, 1984: Table 1 for S. ruberrimus (yelloweye)
tl.to.sl.int <- -0.1717
tl.to.sl.slope <- 0.826
R_Linf = 63.9
rockfish.sl <- rockfish$length.at.age * tl.to.sl.slope + tl.to.sl.int # Why the conversion? In terms of nutrition?
rockfish.Linf.sl <- R_Linf * tl.to.sl.slope + tl.to.sl.int
  # rockfish bins are in standard length
rockfish.bins <- c(rockfish.sl[1] - (rockfish.sl[2]-rockfish.sl[1])/2, 
                    # beginning of 1st size bin is 
                    # mean length at age 1 - (growth from age 1 to age 2)/2
                    # i.e., force mean length at age 1 as midpoint of bin
                    (rockfish.sl[2:length(rockfish.sl)] + rockfish.sl[1:(length(rockfish.sl)-1)])/2,
                    rockfish.Linf.sl)
                    # end of last size bin is L infinity

binned.size.spec <- matrix(data = NA, nrow = rockfish$nage, ncol = 2*lingcod$nage)
col_name <- c(with(lingcod,c(paste0("LF_", age), paste0("LM_", age)))) # names for column (lingcod sex x ages)
row_name <- with(rockfish, paste0("R_", age)) # names for rows (rockfish ages)
rownames(binned.size.spec) <- row_name ; colnames(binned.size.spec) <-  col_name


for(ling.size in 1:length(lingcod$length.at.age)) {
  # calculate area of gamma distribution/diet spectrum within each prey size bin
  gamma.cdf <- pgamma(rockfish.bins, shape=size.spec.al, 
  scale=size.spec.be*lingcod$length.at.age[ling.size]) 
  unnorm.spec <- gamma.cdf[2:length(gamma.cdf)] - 
  gamma.cdf[1:(length(gamma.cdf)-1)] 
  spec.by.num <- unnorm.spec/sum(unnorm.spec)
  # Convert size spectrum by number into size spectrum by mass!
  binned.size.spec[,ling.size] <- spec.by.num # *rockfish$weight.at.age / sum(spec.by.num*rockfish$weight.at.age)
  }

# Check to see if all columns sum to 1
colSums(binned.size.spec)

save(binned.size.spec, file = "cleaned_data/binned.size.spec.Rdata")


# size spectra graph
size.spec.graph = matrix(data = NA, nrow = rockfish$nage, ncol = 100)
for(ling.size in 20:120) {
  # calculate area of gamma distribution/diet spectrum within each prey size bin
  gamma.cdf <- pgamma(rockfish.bins, shape=size.spec.al, 
                      scale=size.spec.be*ling.size) 
  unnorm.spec <- gamma.cdf[2:length(gamma.cdf)] - 
    gamma.cdf[1:(length(gamma.cdf)-1)] 
  spec.by.num <- unnorm.spec/sum(unnorm.spec)
  # Convert size spectrum by number into size spectrum by mass!
  size.spec.graph[,(ling.size-20)] <- spec.by.num*rockfish$weight.at.age / 
    sum(spec.by.num*rockfish$weight.at.age)
}






