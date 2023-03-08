library(parallel)
library("randomForest", lib.loc = "../R_Packages/R3.6.3")
library("SimDesign", lib.loc = "../R_Packages/R3.6.3")
library("purrr", lib.loc = "../R_Packages/R3.6.3")
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/functions_for_GSA.R")


# sample from densities for parameters
m = 2000

handling = runif(m, min = 0.1, max = 1)
corr = runif(m, min = -0.8, max = 0.8)
quant95 = runif(m, min = 0.27, max = 0.33)
rockfish.prop = runif(m, min = 0.01, max = 0.05)
yelloweye.prop = runif(m, min = 0.01, max = 0.05)
cv = runif(m, min = 0.3, max = 0.7)
autocorr_L = runif(m, min = 0.1, max = 0.7)
autocorr_R = runif(m, min = 0.1, max = 0.7)

args <- commandArgs(TRUE)
args <- as.numeric(args[1])


#Run Monte Carlo simulations
sampled_pop <- mcmapply(get_pop_stoch, corr = corr[(1:100) + (args * 100)], 
                                       handling = handling[(1:100) + (args * 100)],
                                       autocorr_L = autocorr_L[(1:100) + (args * 100)],
                                       autocorr_R = autocorr_R[(1:100) + (args * 100)],
                                       cv = cv[(1:100) + (args * 100)],
                                       quant95 = quant95[(1:100) + (args * 100)],
                                       rockfish.prop = rockfish.prop[(1:100) + (args * 100)],
                                       yelloweye.prop = yelloweye.prop[(1:100) + (args * 100)],
         MoreArgs = list(f = 0, b = 0, hist.f = 0, rockfish = rockfish, lingcod = lingcod, nsim = 100, tf = 350), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)

#save(sampled_pop, file = paste0("data/sampled_pop_", args,".Rdata"))

rockfish_avg = unname(unlist(lapply(sampled_pop, `[`, "rockfish_avg")))
rockfish_avg[is.na(rockfish_avg)] <- 0

rockfish_cv = unname(unlist(lapply(sampled_pop, `[`, "rockfish_cv")))
rockfish_cv[is.na(rockfish_cv)] <- 0

rockfish_age = unname(unlist(lapply(sampled_pop, `[`, "rockfish_age")))
rockfish_age[is.na(rockfish_age)] <- 0

# Run simulations over a random forest
RF_SBeq_data <- data.frame(rockfish_avg=rockfish_avg,
                           corr = corr[(1:100) + (args * 100)], 
                           handling = handling[(1:100) + (args * 100)],
                           autocorr_L = autocorr_L[(1:100) + (args * 100)],
                           autocorr_R = autocorr_R[(1:100) + (args * 100)],
                           cv = cv[(1:100) + (args * 100)],
                           quant95 = quant95[(1:100) + (args * 100)],
                           rockfish.prop = rockfish.prop[(1:100) + (args * 100)],
                           yelloweye.prop = yelloweye.prop[(1:100) + (args * 100)])

RF_cv_data <- data.frame(rockfish_cv=rockfish_cv,
                           corr = corr[(1:100) + (args * 100)], 
                           handling = handling[(1:100) + (args * 100)],
                           autocorr_L = autocorr_L[(1:100) + (args * 100)],
                           autocorr_R = autocorr_R[(1:100) + (args * 100)],
                           cv = cv[(1:100) + (args * 100)],
                           quant95 = quant95[(1:100) + (args * 100)],
                           rockfish.prop = rockfish.prop[(1:100) + (args * 100)],
                           yelloweye.prop = yelloweye.prop[(1:100) + (args * 100)])


RF_age_data <- data.frame(rockfish_age=rockfish_age,
                           corr = corr[(1:100) + (args * 100)],
                           handling = handling[(1:100) + (args * 100)],
                           autocorr_L = autocorr_L[(1:100) + (args * 100)],
                           autocorr_R = autocorr_R[(1:100) + (args * 100)],
                           cv = cv[(1:100) + (args * 100)],
                           quant95 = quant95[(1:100) + (args * 100)],
                           rockfish.prop = rockfish.prop[(1:100) + (args * 100)],
                           yelloweye.prop = yelloweye.prop[(1:100) + (args * 100)])

save(RF_SBeq_data, file = paste0("data/RF_SBeq_data_", args,".Rdata"))
save(RF_cv_data, file = paste0("data/RF_cv_data_", args,".Rdata"))
save(RF_age_data, file = paste0("data/RF_age_data_", args,".Rdata"))
