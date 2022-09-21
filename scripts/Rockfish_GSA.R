library(parallel)
library(randomForest) 
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/functions_for_GSA.R")
library(purrr)


# sample from densities for parameters
m = 2000

f = runif(m, min = 0, max = 0.3)
b = runif(m, min = 0, max = 0.2)
handling = runif(m, min = 0.1, max = 1)
corr = runif(m, min = -0.8, max = 0.8)
quant95 = runif(m, min = 0.27, max = 0.33)
rockfish.prop = runif(m, min = 0.01, max = 0.05)
yelloweye.prop = runif(m, min = 0.01, max = 0.05)
cv = runif(m, min = 0.3, max = 0.7)
autocorr_L = runif(m, min = 0.1, max = 0.7)
autocorr_R = runif(m, min = 0.1, max = 0.7)

args <- commandArgs(TRUE)

#Run Monte Carlo simulations
sampled_pop <- mcmapply(get_pop_stoch, f = f[(1:100) + (args[1] * 100)],
                                       b = b[(1:100) + (args[1] * 100)], 
                                       corr = corr[(1:100) + (args[1] * 100)], 
                                       handling = handling[(1:100) + (args[1] * 100)],
                                       autocorr_L = autocorr_L[(1:100) + (args[1] * 100)],
                                       autocorr_R = autocorr_R[(1:100) + (args[1] * 100)],
                                       cv = cv[(1:100) + (args[1] * 100)],
                                       quant95 = quant95[(1:100) + (args[1] * 100)],
                                       rockfish.prop = rockfish.prop[(1:100) + (args[1] * 100)],
                                       yelloweye.prop = yelloweye.prop[(1:100) + (args[1] * 100)],
         MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 100, tf = 350), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)

save(sampled_pop, file = paste0("data/sampled_pop_", args[1],".Rdata"))

# 1000 simulations 350tf: rockfish avg = 3831  rockfish cv = 8.536
# 100 simulations 350tf: rockfish avg = 3821  rockfish cv = 8.312
# 2 simulations 700tf: rockfish avg = 3860  rockfish cv = 10.50
# 2 simulations 1000tf: rockfish avg = 3962  rockfish cv = 10.33
# 100 simuations 500tf: rockfish avg = 3828  rockfish cv = 8.98
# 100 simulations 450tf: rockfish avg = 3831  rockfish cv = 8.85
# 50 simulations 700tf: rockfish avg = 3834  rockfish cv = 9.32


# save.image(file='MonteCarloSimulations.RData')

rockfish_avg = unname(unlist(lapply(sampled_pop, `[`, "rockfish_avg")))
rockfish_avg[is.na(rockfish_avg)] <- 0

rockfish_cv = unname(unlist(lapply(sampled_pop, `[`, "rockfish_cv")))
rockfish_cv[is.na(rockfish_cv)] <- 0


# Run simulations over a random forest
RF_SBeq_data <- data.frame(rockfish_avg=rockfish_avg,
                           f = f[(1:100) + (args[1] * 100)],
                           b = b[(1:100) + (args[1] * 100)], 
                           corr = corr[(1:100) + (args[1] * 100)], 
                           handling = handling[(1:100) + (args[1] * 100)],
                           autocorr_L = autocorr_L[(1:100) + (args[1] * 100)],
                           autocorr_R = autocorr_R[(1:100) + (args[1] * 100)],
                           cv = cv[(1:100) + (args[1] * 100)],
                           quant95 = quant95[(1:100) + (args[1] * 100)],
                           rockfish.prop = rockfish.prop[(1:100) + (args[1] * 100)],
                           yelloweye.prop = yelloweye.prop[(1:100) + (args[1] * 100)])

RF_cv_data <- data.frame(rockfish_cv=rockfish_cv,
                           f = f[(1:100) + (args[1] * 100)],
                           b = b[(1:100) + (args[1] * 100)], 
                           corr = corr[(1:100) + (args[1] * 100)], 
                           handling = handling[(1:100) + (args[1] * 100)],
                           autocorr_L = autocorr_L[(1:100) + (args[1] * 100)],
                           autocorr_R = autocorr_R[(1:100) + (args[1] * 100)],
                           cv = cv[(1:100) + (args[1] * 100)],
                           quant95 = quant95[(1:100) + (args[1] * 100)],
                           rockfish.prop = rockfish.prop[(1:100) + (args[1] * 100)],
                           yelloweye.prop = yelloweye.prop[(1:100) + (args[1] * 100)])

RF <- randomForest(rockfish_avg~.,data=RF_SBeq_data,importance=TRUE,proximity=TRUE)
RF_cv <- randomForest(rockfish_cv~.,data=RF_cv_data,importance=TRUE,proximity=TRUE)

save(RF, file = paste0("data/RF_", args[1],".Rdata"))
save(RF_cv, file = paste0("data/RF_cv_", args[1],".Rdata"))


#get_pop_stoch(rockfish, lingcod, nsim = 10, corr = corr[6], autocorr = c(autocorr[6],autocorr[6]), cv = cv[6],
#             tf = (30+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
#             f = f[6], b = b[6], quant95 = quant95[6], handling = handling[6], min.selectivity = TRUE,
#             min.age.consumed = 4, rockfish.prop = rockfish.prop[6], yelloweye.prop = yelloweye.prop[6])
