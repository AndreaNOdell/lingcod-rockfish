library(parallel)
library("randomForest", lib.loc = "../R_Packages/R3.6.3")
library("SimDesign", lib.loc = "../R_Packages/R3.6.3") 
library("purrr", lib.loc = "../R_Packages/R3.6.3")
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/functions_for_GSA.R")


# sample from densities for parameters
m = 100

f = runif(m, min = 0, max = 0.3)
b = runif(m, min = 0, max = 0.2)
handling = runif(m, min = 0.1, max = 1)
corr = runif(m, min = -0.8, max = 0.8)
quant95 = runif(m, min = 0.27, max = 0.33)
rockfish.prop = runif(m, min = 0.01, max = 0.05)
yelloweye.prop = runif(m, min = 0.01, max = 0.05)
cv = runif(m, min = 0.3, max = 0.7)
autocorr_L = runif(m, min = 0.1, max = 0.4)
autocorr_R = runif(m, min = 0.1, max = 0.4)


#Run Monte Carlo simulations
sampled_pop <- mcmapply(get_pop_stoch, f = f,
                                       b = b, 
                                       corr = corr, 
                                       handling = handling,
                                       autocorr_L = autocorr_L,
                                       autocorr_R = autocorr_R,
                                       cv = cv,
                                       quant95 = quant95,
                                       rockfish.prop = rockfish.prop,
                                       yelloweye.prop = yelloweye.prop,
         MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 100), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)


rockfish_avg = unname(unlist(lapply(sampled_pop, `[`, "rockfish_avg")))
rockfish_avg[is.na(rockfish_avg)] <- 0
rockfish_cv = unname(unlist(lapply(sampled_pop, `[`, "rockfish_cv")))
rockfish_cv[is.na(rockfish_cv)] <- 0

# Run simulations over a random forest
RF_SBeq_data <- data.frame(rockfish_avg=rockfish_avg,
                           f = f,
                           b = b, 
                           corr = corr, 
                           handling = handling,
                           autocorr_L = autocorr_L,
                           autocorr_R = autocorr_R,
                           cv = cv,
                           quant95 = quant95,
                           rockfish.prop = rockfish.prop,
                           yelloweye.prop = yelloweye.prop)

RF_cv_data <- data.frame(rockfish_cv=rockfish_cv,
                           f = f,
                           b = b, 
                           corr = corr, 
                           handling = handling,
                           autocorr_L = autocorr_L,
                           autocorr_R = autocorr_R,
                           cv = cv,
                           quant95 = quant95,
                           rockfish.prop = rockfish.prop,
                           yelloweye.prop = yelloweye.prop)

RF <- randomForest(rockfish_avg~.,data=RF_SBeq_data,importance=TRUE,proximity=TRUE)
RF_cv <- randomForest(rockfish_cv~.,data=RF_cv_data,importance=TRUE,proximity=TRUE)
save(sampled_pop, file='data/sampled_pop.RData')
save(RF, file='data/RF.Rdata')
save(RF_cv, file='data/RF_cv.Rdata')   




