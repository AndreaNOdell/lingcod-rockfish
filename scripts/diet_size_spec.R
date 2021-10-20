# Lingcod diet size spectrum ----------------------------------------------

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

# Paramters from Echeverria & Lenarz, 1984: Table 1 for S. ruberrimus (yelloweye)
tl.to.sl.int <- -0.1717
tl.to.sl.slope <- 0.826
yelloweye.sl <- yelloweye$len.at.age[,'female'] * tl.to.sl.slope + tl.to.sl.int
yelloweye.Linf.sl <- Linf['female'] * tl.to.sl.slope + tl.to.sl.int
yelloweye.bins <- c(yelloweye.sl[1] - (yelloweye.sl[2]-yelloweye.sl[1])/2,
                    # beginning of 1st size bin is 
                    # mean length at age 1 - (growth from age 1 to age 2)/2
                    # i.e., force mean length at age 1 as midpoint of bin
                    (yelloweye.sl[2:length(yelloweye.sl)] + yelloweye.sl[1:(length(yelloweye.sl)-1)])/2,
                    yelloweye.Linf.sl)
# end of last size bin is L infinity

binned.size.spec <- array(0, dim=c(nrow(yelloweye$len.at.age), nrow(lingcod$len.at.age), 2), 
                          dimnames=list(prey.size=NULL, 
                                        lingcod.size=NULL, lingcod.sex=c('female', 'male')))

for(ling.sex in c('female', 'male')) {
  for(ling.size in 1:nrow(lingcod$len.at.age)) {
    # calculate area of gamma distribution/diet spectrum within each prey size bin
    gamma.cdf <- pgamma(yelloweye.bins, shape=size.spec.al, 
                        scale=size.spec.be*lingcod$len.at.age[ling.size, ling.sex]) 
    unnorm.spec <- gamma.cdf[2:length(gamma.cdf)] - 
      gamma.cdf[1:(length(gamma.cdf)-1)] 
    spec.by.num <- unnorm.spec/sum(unnorm.spec)
    # Convert size spectrum by number into size spectrum by mass!
    binned.size.spec[,ling.size, ling.sex] <- spec.by.num*yelloweye$wt.at.age$female / 
      sum(spec.by.num*yelloweye$wt.at.age$female)
  }
}