# Load in files --------
library(tidyverse)
load(file = "results/fishing_scenarios/trial_no_fishing.RData")
load(file = "results/fishing_scenarios/trial_light_fishing.RData")
load(file = "results/fishing_scenarios/trial_moderate_fishing.RData")

light_fishing = trial_light_fishing
moderate_fishing = trial_moderate_fishing
no_fishing = trial_no_fishing

load(file = "results/fishing_scenarios/no_fishing_highpref.RData")
load(file = "results/fishing_scenarios/light_fishing_highpref.RData")
load(file = "results/fishing_scenarios/moderate_fishing_highpref.RData")

rm(list=ls(pattern="trial"))


SBt = 0.4 # Specified SBt/SB0eq value of interest. In most cases, 0.4

# Moderate prey selectivity -------
## No Fishing --------
# now lets calculate the spawning biomass at time t relative to the average equilibrium spawning biomass
relative_sb_nofishing = no_fishing$rockfish_biomass_ts / apply(no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_nofishing), type = "l", ylab = "Bt/B0eq", xlab = "time")
abline(h = 1, lty = 2)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value of B0eq
nofishing_SB40 = as.data.frame(which(relative_sb_nofishing >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40$t_SB40); sd(nofishing_SB40$t_SB40)

recruit_rebuild_df_no = cbind(recruitment = no_fishing$early_recruit, 
                              rebuilding_time = nofishing_SB40$t_SB40)
recruit_rebuild_df_no = as.data.frame(recruit_rebuild_df_no) %>% 
  mutate(fishing = "none")


## Moderate fishing --------
# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_moderatefishing = moderate_fishing$rockfish_biomass_ts / apply(no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
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
recruit_rebuild_df_moderate = cbind(recruitment = moderate_fishing$early_recruit, 
                                    rebuilding_time = moderatefishing_SB40$t_SB40)
recruit_rebuild_df_moderate = as.data.frame(recruit_rebuild_df_moderate)  %>% 
  mutate(fishing = "moderate")

# with row info
recruit_rebuild_df_moderate_with_row = cbind(recruitment = moderate_fishing$early_recruit, 
                                    rebuilding_time = moderatefishing_SB40$t_SB40,
                                    row = moderatefishing_SB40$row)
recruit_rebuild_df_moderate_with_row = as.data.frame(recruit_rebuild_df_moderate_with_row)  %>% 
  mutate(fishing = "moderate")


## Light Fishing ---------
# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_lightfishing = light_fishing$rockfish_biomass_ts / apply(no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_lightfishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
lightfishing_SB40 = as.data.frame(which(relative_sb_lightfishing >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40$t_SB40); sd(lightfishing_SB40$t_SB40)

# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_light = cbind(recruitment = light_fishing$early_recruit, 
                                 rebuilding_time = lightfishing_SB40$t_SB40)
recruit_rebuild_df_light = as.data.frame(recruit_rebuild_df_light) %>% 
  mutate(fishing = "light")



# High Prey Selectivity -----------
## No Fishing --------
# now lets calculate the spawning biomass at time t relative to the average equilibrium spawning biomass
relative_sb_nofishing_high = no_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value of B0eq
nofishing_SB40_high_pref = as.data.frame(which(relative_sb_nofishing_high >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40_high_pref$t_SB40); sd(nofishing_SB40_high_pref$t_SB40)

recruit_rebuild_df_no_highpref = cbind(recruitment = no_fishing_highpref$early_recruit, 
                              rebuilding_time = nofishing_SB40_high_pref$t_SB40)
recruit_rebuild_df_no_highpref = as.data.frame(recruit_rebuild_df_no_highpref) %>% 
  mutate(fishing = "none")


## Moderate fishing --------
# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_moderatefishing_high = moderate_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
max_time = 100 # setting a maximum year for rebuildingto prevent the large standard deviations
moderatefishing_SB40_highpref = as.data.frame(which(relative_sb_moderatefishing_high >= SBt, arr.ind = T)) %>% 
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


## Light Fishing ---------
# Calculate spawning biomass at time t relative to the unfished equilibrium spawning biomass
relative_sb_lightfishing_high = light_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)

# Determine which timestep is the first year post MPA establishment for Bt to reach specified value for B0eq
lightfishing_SB40_highpref = as.data.frame(which(relative_sb_lightfishing_high >= SBt, arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40_highpref$t_SB40); sd(lightfishing_SB40_highpref$t_SB40)

# look at relationship between early recruitment and rebuilding time
recruit_rebuild_df_light_highpref = cbind(recruitment = light_fishing_highpref$early_recruit, 
                                 rebuilding_time = lightfishing_SB40_highpref$t_SB40)
recruit_rebuild_df_light_highpref = as.data.frame(recruit_rebuild_df_light_highpref) %>% 
  mutate(fishing = "light")


# Plot ------------------
recruit_rebuild_df_base = rbind(recruit_rebuild_df_no, recruit_rebuild_df_light, recruit_rebuild_df_moderate) %>%
  mutate(selectivity = "base")
recruit_rebuild_df_highpref = rbind(recruit_rebuild_df_no_highpref, recruit_rebuild_df_light_highpref, recruit_rebuild_df_moderate_highpref) %>%
  mutate(selectivity = "high")
recruit_rebuild_df = rbind(recruit_rebuild_df_base,recruit_rebuild_df_highpref)


ggplot(recruit_rebuild_df, aes(x = recruitment, y = rebuilding_time, color = fishing)) +
  facet_wrap(~ selectivity) +
  geom_point() +
  labs(x = "Recruitment Variability in early years") +
  theme_classic()


SB40_df_base = as.data.frame(rbind(nofishing_SB40, lightfishing_SB40, moderatefishing_SB40[-2])) %>% 
  mutate(selectivity = "base")
SB40_df_high = as.data.frame(rbind(nofishing_SB40_high_pref, lightfishing_SB40_highpref, moderatefishing_SB40_highpref[-2])) %>% 
  mutate(selectivity = "high")
SB40_df = as.data.frame(rbind(SB40_df_base,SB40_df_high))


jpeg("plots/fishing_scenarios/rebuildingtime_SB40_violin_selectivity_compare.jpeg", units="in", width=5, height=4, res = 300)
SB40_df %>% 
  mutate(selectivity = fct_relevel(selectivity, 
                               "base", "high")) %>% 
ggplot(aes(x=reorder(fishing, t_SB40), y=t_SB40, fill = fishing, color = fishing, alpha=selectivity)) + 
  geom_violin(position = "dodge", adjust = .9) +
  ylim(0,100) +
  theme_classic() +
  scale_fill_manual(values=c("dodgerblue", "blue", "blue4"), breaks=c('none', 'light', 'moderate')) +
  scale_color_manual(values=c("dodgerblue", "blue", "blue4"), breaks=c('none', 'light', 'moderate')) +
  scale_alpha_manual(values=c(1, 0)) +
  labs(x = "Fishing Intensity", y = "Years until 0.4B0") + 
  guides(color = FALSE, fill = FALSE, alpha = FALSE)
dev.off()


SB40_summarydf = as.data.frame(cbind(fishing = c("none", "light", "moderate"), 
                                     mean = c(mean(nofishing_SB40$t_SB40), mean(lightfishing_SB40$t_SB40), mean(moderatefishing_SB40$t_SB40)),
                                     sd = c(sd(nofishing_SB40$t_SB40), sd(lightfishing_SB40$t_SB40), sd(moderatefishing_SB40$t_SB40))))
SB40_summarydf$mean = as.numeric(SB40_summarydf$mean); SB40_summarydf$sd = as.numeric(SB40_summarydf$sd)


# disentangling the simulations with high recruitment and long rebuilding time
tmp = recruit_rebuild_df_moderate_with_row %>% 
  filter(rebuilding_time == 100) %>% 
  filter(recruitment > 1.1) 

tmp_row = tmp$row

jpeg("plots/fishing_scenarios/highrecruit_slowrebuild_moderatefishing.jpeg", units="in", width=5, height=4, res = 300)
matplot(t(relative_sb_moderatefishing), type = "l", col = "blue", lty = 1, 
        xlab = "years", ylab = "relative spawning biomass")
matlines(t(relative_sb_moderatefishing[tmp_row,]), col = "red", lty = 1)
abline(v = 100, lty = 2, col = "black")
dev.off()

#compare the low and high prey selectivity trajectories
relative_sb_moderatefishing = moderate_fishing$rockfish_biomass_ts / apply(no_fishing$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_moderatefishing), type = "l", ylab = "Bt/B0eq", xlab = "time")

relative_sb_moderatefishing_high = moderate_fishing_highpref$rockfish_biomass_ts / apply(no_fishing_highpref$rockfish_biomass_ts[,-(1:150)], 1, mean)
matplot(t(relative_sb_moderatefishing_high), type = "l", ylab = "Bt/B0eq", xlab = "time")


