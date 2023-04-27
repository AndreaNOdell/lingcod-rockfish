# For this section, all the parameter values need to match across scenarios except
# the single ecological parameter that's varying and the fishing pressure
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/model.R")

# Fishing focused -----------------------------------
##  Select scenarios --------------------------------------------
SBt = 0.4 # Specified SBt/SB0eq value of interest. In most cases, 0.4
nsims = 100 # number of simulations
cv = 0.5
autocorr_R = 0.23
mpa.year = 0
rec.year = 25 # the max year to calculate average recruitment
rock.prop = 0.15 # 15% of lingcod diet is rockfish

## No Fishing ----------------------------------------------------

trial_no_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                              cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = 0.1, 
                              quant95 = 0.29, handling = 0.23, min.selectivity = TRUE, min.age.consumed = 4, 
                              rockfish.prop = rock.prop, yelloweye.prop = 0.05)
#save(trial_no_fishing, file = "results/fishing_scenarios/trial_no_fishing_highpref.RData")

jpeg("plots/fishing_scenarios/example_timeseries_nofishing.jpeg")
matplot(t(trial_no_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")
dev.off()

# now lets calculate the spawning biomass at time t relative to the average equilibrium spawning biomass
relative_sb_nofishing = trial_no_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_nofishing), type = "l", ylab = "Bt/B0eq", xlab = "time")
abline(h = 1, lty = 2)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value of B0eq
nofishing_SB40 = as.data.frame(which(relative_sb_nofishing >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40$t_SB40); sd(nofishing_SB40$t_SB40)

recruit_rebuild_df_no = cbind(recruitment = trial_no_fishing$early_recruit, 
                              rebuilding_time = nofishing_SB40$t_SB40)
recruit_rebuild_df_no = as.data.frame(recruit_rebuild_df_no) %>% 
  mutate(fishing = "none")


## Moderate Fishing --------------------------------------------------------------

trial_moderate_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                    cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0.3, b = 0.1, 
                                    quant95 = 0.29, handling = 0.23, min.selectivity = TRUE, min.age.consumed = 4, 
                                    rockfish.prop = rock.prop, yelloweye.prop = 0.05)
#save(trial_moderate_fishing, file = "results/fishing_scenarios/trial_moderate_fishing_highpref.RData")


matplot(t(trial_moderate_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_moderatefishing = trial_moderate_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_moderatefishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
max_time = 100 # setting a maximum year for rebuildingto prevent the large standard deviations
moderatefishing_SB40 = as.data.frame(which(relative_sb_moderatefishing >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40_old = min(col)) %>% 
  mutate(t_SB40 = replace(t_SB40_old, t_SB40_old > max_time, max_time)) %>% 
  mutate(fishing = "moderate")
mean(moderatefishing_SB40$t_SB40); sd(moderatefishing_SB40$t_SB40)
sum(moderatefishing_SB40$t_SB40_old > 100)/100 # What fraction of simulations have rebuilding time greater than max_time 

# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_moderate = cbind(recruitment = trial_moderate_fishing$early_recruit, 
                                    rebuilding_time = moderatefishing_SB40$t_SB40)
recruit_rebuild_df_moderate = as.data.frame(recruit_rebuild_df_moderate)  %>% 
  mutate(fishing = "moderate")

# No bycatch pressure, so just looking at how removal of predation alone affects rebuilding
trial_moderate_fishing_nobycatch = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                              cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0.3, b = 0, 
                                              quant95 = 0.29, handling = 0.23, min.selectivity = TRUE, min.age.consumed = 4, 
                                              rockfish.prop = rock.prop, yelloweye.prop = 0.05)
#save(trial_moderate_fishing_nobycatch, file = "results/fishing_scenarios/trial_moderate_fishing_nobycatch_highpref.RData")

relative_sb_moderatefishing_nobycatch = trial_moderate_fishing_nobycatch$rockfish_biomass_ts / apply(trial_moderate_fishing_nobycatch$rockfish_biomass_ts[,-(1:150)], 1, mean)
moderatefishing_nobycatch_SB40 = as.data.frame(which(relative_sb_moderatefishing_nobycatch >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40_old = min(col)) %>% 
  mutate(fishing = "moderate") %>% 
  mutate(t_SB40 = replace(t_SB40_old, t_SB40_old > max_time, max_time))

recruit_rebuild_df_moderate_nobycatch = cbind(recruitment = trial_moderate_fishing_nobycatch$early_recruit, 
                                              rebuilding_time = moderatefishing_nobycatch_SB40$t_SB40_old)
recruit_rebuild_df_moderate_nobycatch = as.data.frame(recruit_rebuild_df_moderate_nobycatch)  %>% 
  mutate(fishing = "moderate no bycatch")
mean(moderatefishing_nobycatch_SB40$t_SB40); sd(moderatefishing_nobycatch_SB40$t_SB40)


## Light Fishing ----------------------------------------------------------------

trial_light_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                 cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = 0.1, 
                                 quant95 = 0.29, handling = 0.23, min.selectivity = TRUE, min.age.consumed = 4, 
                                 rockfish.prop = rock.prop, yelloweye.prop = 0.05)
#save(trial_light_fishing, file = "results/fishing_scenarios/trial_light_fishing_highpref.RData")

matplot(t(trial_light_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_lightfishing = trial_light_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_lightfishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
lightfishing_SB40 = as.data.frame(which(relative_sb_lightfishing >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40$t_SB40); sd(lightfishing_SB40$t_SB40)


# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_light = cbind(recruitment = trial_light_fishing$early_recruit, 
                                 rebuilding_time = lightfishing_SB40$t_SB40)
recruit_rebuild_df_light = as.data.frame(recruit_rebuild_df_light) %>% 
  mutate(fishing = "light")


## plots -----------
recruit_rebuild_df = rbind(recruit_rebuild_df_no, recruit_rebuild_df_light, recruit_rebuild_df_moderate, recruit_rebuild_df_moderate_nobycatch)
ggplot(recruit_rebuild_df, aes(x = recruitment, y = rebuilding_time, color = fishing)) +
  geom_point() +
  labs(x = "Recruitment Variability in early years") +
  theme_classic()

recruit_rebuild_df %>% 
  group_by(fishing) %>% 
  summarise(mean = mean(rebuilding_time), sd = sd(rebuilding_time))

pred_removal_affect <- c(cv = trial_moderate_fishing_nobycatch$rockfish_cv, 
                         sb = trial_moderate_fishing_nobycatch$rockfish_avg,
                         age = trial_moderate_fishing_nobycatch$rockfish_age)
no_fishing <- c(cv = trial_no_fishing$rockfish_cv, 
                sb = trial_no_fishing$rockfish_avg,
                age = trial_no_fishing$rockfish_age)






