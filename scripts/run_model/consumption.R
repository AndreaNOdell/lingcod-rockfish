source("scripts/scratch_script.R")
library(miceadds)
library(parallel)


# First let's alter the slope of the rockfish size spectra vulnerable to rockfish predation ----------
# Load in rockfish size-spectra in size-specific lingcod diet
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) 
load.Rdata( filename="cleaned_data/binned.size.spec.27.Rdata", "binned.size.spec.27" )
load.Rdata( filename="cleaned_data/binned.size.spec.31.Rdata", "binned.size.spec.31" )
load.Rdata( filename="cleaned_data/binned.size.spec.33.Rdata", "binned.size.spec.33" )

as.data.frame(binned.size.spec.33) %>% 
  mutate(r_age = 1:65, .before = LF_1) %>% 
  gather(key = "L_age", value = "preference", LF_1:LM_20, factor_key = TRUE) %>%
  ggplot(aes(x = L_age, y = r_age, fill = preference)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0, 0.9))

# Create scenarios
corr = c(0.8, -0.8)
size.spec = list(binned.size.spec.29, binned.size.spec.27, binned.size.spec.31, binned.size.spec.33)
size.spec.scenario = as.data.frame(expand_grid(corr, size.spec))
size.spec.scenario.names = expand_grid(corr,  c("a_ij.29", "a_ij.27", "a_ij.31", "a_ij.33"))
names(size.spec.scenario.names) = c("corr", "size.spec")

# run model
    # 100 simulations across 4 different size-spectra scenarios
size.spec.run = mcmapply(get_pop_stoch, corr = size.spec.scenario[,1], a_ij = size.spec.scenario[,2],
                MoreArgs = list(rockfish = rockfish, lingcod = lingcod, f = 0.1, b = 0.1, nsim = 100, handling = 0.25), SIMPLIFY = FALSE, mc.cores = 4)

consumption_avg = unname(unlist(lapply(size.spec.run, `[`, "consumption_avg")))
consumption_cv = unname(unlist(lapply(size.spec.run, `[`, "consumption_cv")))
rockfish_avg = unname(unlist(lapply(size.spec.run, `[`, "rockfish_avg")))
rockfish_cv = unname(unlist(lapply(size.spec.run, `[`, "rockfish_cv")))
size.spec.out = cbind(size.spec.scenario.names, consumption_avg, consumption_cv)
size.spec.rockfish.out = cbind(size.spec.scenario.names, rockfish_avg, rockfish_cv)

save(size.spec.out, file = "results/size.spec.out.Rdata")
save(size.spec.rockfish.out, file = "results/size.spec.rockfish.out.Rdata")

#jpeg("plots/consumptionSB.sizespec.jpeg")
size.spec.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, size.spec, consumption_avg) %>% 
  ggplot(aes(x = size.spec, y = consumption_avg, group = correlation, col = correlation)) +
  geom_point() +
  labs(x = "slope of size spectra", y = "avg consumption", 
       title = "avg rockfish consumed") +
  theme_classic()
#dev.off()

#jpeg("plots/consumptionCV.sizespec.jpeg")
size.spec.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, size.spec, consumption_cv) %>% 
  ggplot(aes(x = size.spec, y = consumption_cv, group = correlation, col = correlation)) +
  geom_point() +
  labs(x = "slope of size spectra", y = "avg consumption cv", 
       title = "cv rockfish consumption") +
  theme_classic()
#dev.off()

jpeg("plots/rockfishSB.sizespec.jpeg")
size.spec.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, size.spec, rockfish_avg) %>% 
  ggplot(aes(x = size.spec, y = rockfish_avg, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "slope of size spectra", y = "avg SB", 
       title = "Rockfish SB") +
  theme_classic()
dev.off()

jpeg("plots/rockfishCV.sizespec.jpeg")
size.spec.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, size.spec, rockfish_cv) %>% 
  ggplot(aes(x = size.spec, y = rockfish_cv, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "slope of size spectra", y = "avg CV", 
       title = "Rockfish cv") +
  theme_classic()
dev.off()



# Alter handling time ------------------------------------------------------

handling_trial = c(0.1, 0.3, 0.5, 0.7)
handling_scenarios = as.data.frame(expand_grid(handling_trial, corr))

handling.run = mcmapply(get_pop_stoch, handling = handling_scenarios[,1], corr = handling_scenarios[,2],
         MoreArgs = list(rockfish = rockfish, lingcod = lingcod, f = 0.1, b = 0.1, nsim = 100), SIMPLIFY = FALSE, mc.cores = 4)


consumption_avg = unname(unlist(lapply(handling.run, `[`, "consumption_avg")))
consumption_cv = unname(unlist(lapply(handling.run, `[`, "consumption_cv")))
rockfish_avg = unname(unlist(lapply(handling.run, `[`, "rockfish_avg")))
rockfish_cv = unname(unlist(lapply(handling.run, `[`, "rockfish_cv")))
handling.out = cbind(handling_scenarios, consumption_avg, consumption_cv)
handling.rockfish.out = cbind(handling_scenarios, rockfish_avg, rockfish_cv)

save(handling.out, file = "results/handling.out.Rdata")
save(handling.rockfish.out, file = "results/handling.rockfish.out.Rdata")


jpeg("plots/consumptionSB.handling.jpeg")
handling.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, handling_trial, consumption_avg) %>% 
  ggplot(aes(x = handling_trial, y = consumption_avg, group = correlation, col = correlation)) +
  geom_point(size = 3) +
  labs(x = "handling rate", y = "avg consumption", 
       title = "Rockfish avg consumption") +
  theme_classic()
dev.off()

jpeg("plots/consumptionCV.handling.jpeg")
handling.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, handling_trial, consumption_cv) %>% 
  ggplot(aes(x = handling_trial, y = consumption_cv, group = correlation, col = correlation)) +
  geom_point(size = 3) +
  labs(x = "handling rate", y = "avg consumption cv", 
       title = "Rockfish consumption cv") +
  theme_classic()
dev.off()

jpeg("plots/rockfishSB.handling.jpeg")
handling.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, handling_trial, rockfish_avg) %>% 
  ggplot(aes(x = handling_trial, y = rockfish_avg, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "handling rate", y = "avg SB", 
       title = "Rockfish SB") +
  theme_classic()
dev.off()

jpeg("plots/rockfishCV.handling.jpeg")
handling.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, handling_trial, rockfish_cv) %>% 
  ggplot(aes(x = handling_trial, y = rockfish_cv, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "handling rate", y = "avg CV", 
       title = "Rockfish cv") +
  theme_classic()
dev.off()


# yelloweye proportion ------------------------------------------------------

yelloweye_prop_trials = c(0.04, 0.05, 0.06)
yelloweye_prop_scenarios = as.data.frame(expand_grid(yelloweye_prop_trials, corr))

ye.prop.run = mcmapply(get_pop_stoch, yelloweye.prop = yelloweye_prop_scenarios[,1], corr = yelloweye_prop_scenarios[,2],
         MoreArgs = list(rockfish = rockfish, lingcod = lingcod, f = 0.1, b = 0.1, nsim = 100, 
                         handling = 0.25), SIMPLIFY = FALSE, mc.cores = 4)


consumption_avg = unname(unlist(lapply(ye.prop.run, `[`, "consumption_avg")))
consumption_cv = unname(unlist(lapply(ye.prop.run, `[`, "consumption_cv")))
rockfish_avg = unname(unlist(lapply(ye.prop.run, `[`, "rockfish_avg")))
rockfish_cv = unname(unlist(lapply(ye.prop.run, `[`, "rockfish_cv")))
ye.prop.out = cbind(yelloweye_prop_scenarios, consumption_avg, consumption_cv)
ye.prop.rockfish.out = cbind(yelloweye_prop_scenarios, rockfish_avg, rockfish_cv)

save(ye.prop.out, file = "results/ye.prop.out.Rdata")
save(ye.prop.rockfish.out, file = "results/ye.prop.rockfish.out.Rdata")


jpeg("plots/consumptionSB.yeprop.jpeg")
ye.prop.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, yelloweye_prop_trials, consumption_avg) %>% 
  ggplot(aes(x = yelloweye_prop_trials, y = consumption_avg, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "proportion yelloweye", y = "avg consumption", 
       title = "Rockfish avg consumption") +
  theme_classic()
dev.off()

jpeg("plots/consumptionCV.yeprop.jpeg")
ye.prop.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, yelloweye_prop_trials, consumption_cv) %>% 
  ggplot(aes(x = yelloweye_prop_trials, y = consumption_cv, group = correlation, col = correlation)) +
  geom_point(size = 3) +
  labs(x = "proportion yelloweye", y = "avg consumption cv", 
       title = "Rockfish consumption cv") +
  theme_classic()
dev.off()

jpeg("plots/rockfishSB.yeprop.jpeg")
ye.prop.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, yelloweye_prop_trials, rockfish_avg) %>% 
  ggplot(aes(x = yelloweye_prop_trials, y = rockfish_avg, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "proportion yelloweye", y = "avg SB", 
       title = "Rockfish SB") +
  theme_classic()
dev.off()

jpeg("plots/rockfishCV.yeprop.jpeg")
ye.prop.rockfish.out %>% 
  mutate(correlation = recode(corr, '-0.8' = 'neg', '0.8' = 'pos')) %>% 
  select(correlation, yelloweye_prop_trials, rockfish_cv) %>% 
  ggplot(aes(x = yelloweye_prop_trials, y = rockfish_cv, group = correlation, col = correlation)) +
  geom_point(size = 5) +
  labs(x = "proportion yelloweye", y = "avg CV", 
       title = "Rockfish CV") +
  theme_classic()
dev.off()





