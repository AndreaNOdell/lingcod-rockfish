# For this section, all the parameter values need to match across scenarios except
# the single ecological parameter that's varying and the fishing pressure
load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/model.R")

# Fishing focused -----------------------------------
##  Select scenarios --------------------------------------------
SBt = 0.4 # Specified SBt/SB0eq value of interest. In most cases, 0.4
nsims = 150 # number of simulations
cv = 0.5
autocorr_R = 0.23
mpa.year = 0
rec.year = 25
handl = 0.3
bycatch = 0.05
rock.prop = 0.1
yellow.prop = 1

## No Fishing ----------------------------------------------------

no_fishing_highpref = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                              cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = bycatch, 
                              quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                              rockfish.prop = rock.prop, yelloweye.prop = yellow.prop)
#save(no_fishing_highpref, file = "results/fishing_scenarios/no_fishing_highpref.RData")

jpeg("plots/fishing_scenarios/example_timeseries_nofishing.jpeg")
matplot(t(no_fishing_highpref$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")
dev.off()

# now lets calculate the spawning biomass at time t relative to the average equilibrium spawning biomass
relative_sb_nofishing_highpref = no_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_nofishing_highpref), type = "l", ylab = "Bt/B0eq", xlab = "time")
abline(h = 1, lty = 2)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value of B0eq
nofishing_SB40_high_pref = as.data.frame(which(relative_sb_nofishing_highpref >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40_high_pref$t_SB40); sd(nofishing_SB40_high_pref$t_SB40)

recruit_rebuild_df_no_highpref = cbind(recruitment = no_fishing_highpref$early_recruit, 
                              rebuilding_time = nofishing_SB40_high_pref$t_SB40)
recruit_rebuild_df_no_highpref  = as.data.frame(recruit_rebuild_df_no_highpref ) %>% 
  mutate(fishing = "none")


## Moderate Fishing --------------------------------------------------------------

moderate_fishing_highpref = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                    cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0.3, b = bycatch, 
                                    quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                                    rockfish.prop = rock.prop, yelloweye.prop = yellow.prop)
#save(moderate_fishing_highpref, file = "results/fishing_scenarios/moderate_fishing_highpref.RData")

matplot(t(moderate_fishing_highpref$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_moderatefishing_highpref = moderate_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_moderatefishing_highpref), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
max_time = 100 # setting a maximum year for rebuildingto prevent the large standard deviations
moderatefishing_SB40_highpref = as.data.frame(which(relative_sb_moderatefishing_highpref >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40_old = min(col)) %>% 
  mutate(t_SB40 = replace(t_SB40_old, t_SB40_old > max_time, max_time)) %>% 
  mutate(fishing = "moderate")
mean(moderatefishing_SB40_highpref$t_SB40); sd(moderatefishing_SB40_highpref$t_SB40)
sum(moderatefishing_SB40_highpref$t_SB40_old > 100)/100 # What fraction of simulations have rebuilding time greater than max_time 

# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_moderate_highpref = cbind(recruitment = moderate_fishing_highpref$early_recruit, 
                                    rebuilding_time = moderatefishing_SB40_highpref$t_SB40)
recruit_rebuild_df_moderate_highpref = as.data.frame(recruit_rebuild_df_moderate_highpref)  %>% 
  mutate(fishing = "moderate")


## Light Fishing ----------------------------------------------------------------

light_fishing_highpref= get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                 cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = bycatch, 
                                 quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                                 rockfish.prop = rock.prop, yelloweye.prop = yellow.prop)
#save(light_fishing_highpref, file = "results/fishing_scenarios/light_fishing_highpref.RData")

matplot(t(light_fishing_highpref$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_lightfishing_highpref = light_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_lightfishing_highpref), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
lightfishing_SB40_highpref = as.data.frame(which(relative_sb_lightfishing_highpref >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40_highpref$t_SB40); sd(lightfishing_SB40_highpref$t_SB40)


# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_light_highpref = cbind(recruitment = light_fishing_highpref$early_recruit, 
                                 rebuilding_time = lightfishing_SB40_highpref$t_SB40)
recruit_rebuild_df_light_highpref = as.data.frame(recruit_rebuild_df_light_highpref) %>% 
  mutate(fishing = "light")


## plots -----------
recruit_rebuild_df_highpref = rbind(recruit_rebuild_df_no_highpref, recruit_rebuild_df_light_highpref, recruit_rebuild_df_moderate_highpref)
ggplot(recruit_rebuild_df_highpref, aes(x = recruitment, y = rebuilding_time, color = fishing)) +
  geom_point() +
  labs(x = "Recruitment Variability in early years") +
  theme_classic()

recruit_rebuild_df %>% 
  group_by(fishing) %>% 
  summarise(mean = mean(rebuilding_time), sd = sd(rebuilding_time))

pred_removal_affect <- c(cv = trial_moderate_fishing_nobycatch$rockfish_cv, 
                         sb = trial_moderate_fishing_nobycatch$rockfish_avg,
                         age = trial_moderate_fishing_nobycatch$rockfish_age)
no_fishing <- c(cv = no_fishing_highpref$rockfish_cv, 
                sb = trial_no_fishing$rockfish_avg,
                age = trial_no_fishing$rockfish_age)






