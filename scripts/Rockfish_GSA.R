library(parallel)
library(randomForest) 
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/functions_for_GSA.R")


# sample from densities for parameters
m = 20

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
         MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10), SIMPLIFY = FALSE, mc.cores = parallel::detectCores()-1)

# cores <- detectCores()
# cl <- makeCluster(cores)
# clusterExport(cl, ls())
# clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))
# sampled_pop <- parSapply(cl=cl,1:m,function(i){ 
#   get_pop_stoch(f = f[i],
#                 b = b[i], 
#                 corr = corr[i], 
#                 handling = handling[i],
#                 autocorr_L = autocorr_L[i],
#                 autocorr_R = autocorr_R[i],
#                 cv = cv[i],
#                 quant95 = quant95[i],
#                 rockfish.prop = rockfish.prop[i],
#                 yelloweye.prop = yelloweye.prop[i],
#                 rockfish = rockfish, lingcod = lingcod, nsim = 10, min.selectivity = TRUE,
#                 tf = (30+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, min.age.consumed = 4)
# })

save.image(file='MonteCarloSimulations.RData')

rockfish_avg = unname(unlist(lapply(sampled_pop, `[`, "rockfish_avg")))
rockfish_cv = unname(unlist(lapply(sampled_pop, `[`, "rockfish_cv")))




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
save.image(file='PredPrey_GSA.RData')



#get_pop_stoch(rockfish, lingcod, nsim = 10, corr = corr[6], autocorr = c(autocorr[6],autocorr[6]), cv = cv[6],
#             tf = (30+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
#             f = f[6], b = b[6], quant95 = quant95[6], handling = handling[6], min.selectivity = TRUE,
#             min.age.consumed = 4, rockfish.prop = rockfish.prop[6], yelloweye.prop = yelloweye.prop[6])
