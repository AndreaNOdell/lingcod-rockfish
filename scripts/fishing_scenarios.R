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
rec.year = 25

## No Fishing ----------------------------------------------------

trial_no_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                   cv = cv, tf = (5+mpa.year+300), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = 0.1, 
                   quant95 = 0.29, handling = 0.23, min.selectivity = TRUE, min.age.consumed = 4, 
                   rockfish.prop = 0.05, yelloweye.prop = 0.05)
  #save(trial_no_fishing, file = "results/fishing_scenarios/trial_no_fishing.RData")

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
                              rockfish.prop = 0.05, yelloweye.prop = 0.05)
  #save(trial_moderate_fishing, file = "results/fishing_scenarios/trial_moderate_fishing.RData")


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
                                    rockfish.prop = 0.05, yelloweye.prop = 0.05)
  #save(trial_moderate_fishing_nobycatch, file = "results/fishing_scenarios/trial_moderate_fishing_nobycatch.RData")

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
                                    rockfish.prop = 0.05, yelloweye.prop = 0.05)
  #save(trial_light_fishing, file = "results/fishing_scenarios/trial_light_fishing.RData")

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
trial_light_fishing$early_recruit

# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_light = cbind(recruitment = trial_light_fishing$early_recruit, 
                                    rebuilding_time = lightfishing_SB40$t_SB40)
recruit_rebuild_df_light = as.data.frame(recruit_rebuild_df_light) %>% 
  mutate(fishing = "light")


## Compare and Contrast plots ---------------------------------------------------------------------
# plotting relationship between recruitment and rebuilding time across fishing scenarios
recruit_rebuild_df = rbind(recruit_rebuild_df_no, recruit_rebuild_df_light, recruit_rebuild_df_moderate, recruit_rebuild_df_moderate_nobycatch)
jpeg("plots/fishing_scenarios/recruit_rebuild.jpeg", units="in", width=5, height=4, res = 300)
ggplot(recruit_rebuild_df, aes(x = recruitment, y = rebuilding_time, color = fishing)) +
  geom_point() +
  labs(x = "Recruitment Variability in early years") +
  theme_classic()
dev.off()
    # This plot includes both predation and bycatch mortality on rockfish. 
    # Might be worth looking at no bycatch, to disentangle the difference
    # between bycatch mortality and predation effects.

# plotting timeseries of relative spawning biomass across fishing scenarios
jpeg("plots/fishing_scenarios/relativeSB_timeseries.jpeg", units="in", width=5, height=4, res = 300)
matplot(t(relative_sb_nofishing), type = "l", col = "blue", ylab = "SBt/SB0eq", xlab = "time", main = "Relative SB timeseries")
matlines(t(relative_sb_moderatefishing), col = "red")
matlines(t(relative_sb_lightfishing), col = "orange")
dev.off()


# combine minimum time to SB40 dataframes
SB40_df = as.data.frame(rbind(nofishing_SB40, lightfishing_SB40, moderatefishing_SB40[-2]))
SB40_summarydf = as.data.frame(cbind(fishing = c("none", "light", "moderate"), 
                       mean = c(mean(nofishing_SB40$t_SB40), mean(lightfishing_SB40$t_SB40), mean(moderatefishing_SB40$t_SB40)),
                       sd = c(sd(nofishing_SB40$t_SB40), sd(lightfishing_SB40$t_SB40), sd(moderatefishing_SB40$t_SB40))))
SB40_summarydf$mean = as.numeric(SB40_summarydf$mean); SB40_summarydf$sd = as.numeric(SB40_summarydf$sd)

# Rebuilding time to 40% unfished equilibrium (violin plot)
jpeg("plots/fishing_scenarios/rebuildingtime_SB40_violin.jpeg", units="in", width=5, height=4, res = 300)
ggplot(SB40_df, aes(x=reorder(fishing, t_SB40), y=t_SB40)) + 
  geom_violin() +
  ylim(0,100) +
  theme_classic() +
  labs(x = "Fishing Intensity", y = "Years until 0.4B0")
dev.off()

# rebuilding time to 40% unfished equilibrium (error bar plot)
jpeg("plots/fishing_scenarios/rebuildingtime_SB40.jpeg", units="in", width=5, height=4, res = 300)
ggplot(SB40_summarydf, aes(x = reorder(fishing, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  theme_classic() +
  labs(title = "Rebuilding time", x = "fishing pressure", y = "years till 0.4 B0eq")
dev.off()

# What biological attribute might be more prone to longer rebuilding time under moderate fishing pressure


# dataframe containing the spawning biomass, stability, and age-structure across scenarios
bio_out_df = as.data.frame(cbind(fishing = c("none", "light", "moderate"),
                   sb = c(trial_no_fishing$rockfish_avg/trial_no_fishing$rockfish_avg, trial_light_fishing$rockfish_avg/trial_no_fishing$rockfish_avg, trial_moderate_fishing$rockfish_avg/trial_no_fishing$rockfish_avg),
                   cv = c(trial_no_fishing$rockfish_cv/trial_no_fishing$rockfish_cv, trial_light_fishing$rockfish_cv/trial_no_fishing$rockfish_cv, trial_moderate_fishing$rockfish_cv/trial_no_fishing$rockfish_cv),
                   age = c(trial_no_fishing$rockfish_age/trial_no_fishing$rockfish_age, trial_light_fishing$rockfish_age/trial_no_fishing$rockfish_age, trial_moderate_fishing$rockfish_age/trial_no_fishing$rockfish_age)
                   ))
bio_out_df$sb = as.numeric(bio_out_df$sb); bio_out_df$cv = as.numeric(bio_out_df$cv); bio_out_df$age = as.numeric(bio_out_df$age)

# Plot showing change in biological outcomes (age structure, variability, and biomass) across fishing scenarios relative to unfished
jpeg("plots/fishing_scenarios/relative_impact_bio_fishing.jpeg", units="in", width=4, height=4, res = 300)
bio_out_df %>% 
  pivot_longer(!fishing, names_to = "outcome", values_to = "value") %>% 
  ggplot(aes(x = value, y = outcome, color = fishing)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_point() +
  theme_classic() +
  lims(x = c(0,2)) +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) +
  labs(x = "change relative to unfished", title = "Fishing impact on biological outcomes") +
  scale_color_manual(values=c("dodgerblue", "blue", "blue4"), breaks=c('none', 'light', 'moderate'))
dev.off()


## Rebuilding time associations with age-structure and stability -------
  # merge rebuilding time data with the age-structure and stability
  # plot


# Handling focused (under construction do not run) ---------------------- 
##  Select scenarios --------------------------------------------
SBt = 0.4 # SBt/SB0eq value
nsims = 10 # number of simulations
cv = 0.5
autocorr_R = 0.23
handling = 0.1

## No Fishing ----------------------------------------------------

trial_no_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                              cv = cv, tf = (5+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, f = 0, b = 0.1, 
                              quant95 = 0.29, handling = handling, min.selectivity = TRUE, min.age.consumed = 4, 
                              rockfish.prop = 0.05, yelloweye.prop = 0.05)

jpeg("plots/fishing_scenarios/example_timeseries_nofishing.jpeg")
matplot(t(trial_no_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")
dev.off()

# now lets calculate the spawning biomass at time t relative to the average equilibrium spawning biomass
relative_sb_nofishing = trial_no_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_nofishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value of B0eq
nofishing_SB40 = as.data.frame(which(relative_sb_nofishing >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col))
mean(nofishing_SB40$t_SB40); sd(nofishing_SB40$t_SB40)

## Moderate Fishing --------------------------------------------------------------

trial_moderate_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                    cv = cv, tf = (5+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, f = 0.3, b = 0.1, 
                                    quant95 = 0.29, handling = handling, min.selectivity = TRUE, min.age.consumed = 4, 
                                    rockfish.prop = 0.05, yelloweye.prop = 0.05)

matplot(t(trial_moderate_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_moderatefishing = trial_moderate_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_moderatefishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
moderatefishing_SB40 = as.data.frame(which(relative_sb_moderatefishing >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col))
mean(moderatefishing_SB40$t_SB40); sd(moderatefishing_SB40$t_SB40)


## Light Fishing ----------------------------------------------------------------

trial_light_fishing = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                                 cv = cv, tf = (5+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = 0.1, 
                                 quant95 = 0.29, handling = handling, min.selectivity = TRUE, min.age.consumed = 4, 
                                 rockfish.prop = 0.05, yelloweye.prop = 0.05)

matplot(t(trial_light_fishing$rockfish_biomass_ts), type = "l", ylab = "spawning biomass", xlab = "time", main = "Example time-series")

# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_lightfishing = trial_light_fishing$rockfish_biomass_ts / apply(trial_no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_lightfishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
lightfishing_SB40 = as.data.frame(which(relative_sb_lightfishing >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col))
mean(lightfishing_SB40$t_SB40); sd(lightfishing_SB40$t_SB40)


## Compare and Contrast plots ---------------------------------------------------------------------

# plotting timeseries of relative spawning biomass across fishing scenarios
jpeg("plots/fishing_scenarios/relativeSB_timeseries.jpeg", units="in", width=5, height=4, res = 300)
matplot(t(relative_sb_nofishing), type = "l", col = "blue", ylab = "SBt/SB0eq", xlab = "time", main = "Relative SB timeseries")
matlines(t(relative_sb_moderatefishing), col = "red")
matlines(t(relative_sb_lightfishing), col = "orange")
dev.off()


# combine minimum time to SB40 dataframes
lightfishing_SB40 = lightfishing_SB40 %>% 
  mutate(fishing = "light")
moderatefishing_SB40 = moderatefishing_SB40 %>% 
  mutate(fishing = "moderate")
nofishing_SB40 = nofishing_SB40 %>% 
  mutate(fishing = "none")
SB40_df = as.data.frame(rbind(nofishing_SB40, lightfishing_SB40, moderatefishing_SB40))
SB40_summarydf = as.data.frame(cbind(fishing = c("none", "light", "moderate"), 
                                     mean = c(mean(nofishing_SB40$t_SB40), mean(lightfishing_SB40$t_SB40), mean(moderatefishing_SB40$t_SB40)),
                                     sd = c(sd(nofishing_SB40$t_SB40), sd(lightfishing_SB40$t_SB40), sd(moderatefishing_SB40$t_SB40))))
SB40_summarydf$mean = as.numeric(SB40_summarydf$mean); SB40_summarydf$sd = as.numeric(SB40_summarydf$sd)

ggplot(SB40_summarydf, aes(x = reorder(fishing, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  theme_classic()

