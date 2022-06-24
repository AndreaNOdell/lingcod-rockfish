source("scripts/scratch_script.R")
library(miceadds)
library(parallel)
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet


# Unfished equilibrium ---------------------------------------------------------------------------------
unfished_vals = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 50,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
                            f = 0, b = 0, a_ij = binned.size.spec.29, handling = 0.7, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough
unfished_vals_no_int = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
                                   f = 0, b = 0, a_ij = a_ij.0, handling = 0.1, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough


# Check if yield is same between the two harvest strategies -----------------------------------------------
fishing_trials = seq(0.1, 0.2, by = 0.1)
lingcod_yield_check_MSL = mapply(get_pop_det, f = fishing_trials, 
                             MoreArgs = list(b = 0.1, rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE)
lingcod_yield_check_HS = mapply(get_pop_det, f = fishing_trials, 
                                 MoreArgs = list(b = 0.1, rockfish = rockfish, lingcod = lingcod, min.selectivity = FALSE), SIMPLIFY = FALSE)

yield_check = matrix(NA, nrow = length(fishing_trials), ncol = 2)
for(j in 1:2) {
yield_check[j,1] = sum(lingcod_yield_check_MSL[[j]]$yield[1,300:470])
yield_check[j,2] = sum(lingcod_yield_check_HS[[j]]$yield[1,300:470])
}

ling_eq_check = matrix(NA, nrow = length(fishing_trials), ncol = 2)
colnames(ling_eq_check) = c("MSL", "HS")
for(j in 1:2) {
  ling_eq_check[j,1] = lingcod_yield_check_MSL[[j]]$lingcod_SBeq
 ling_eq_check[j,2] = lingcod_yield_check_HS[[j]]$lingcod_SBeq
}

ling_eq_check = cbind(fishing_trials, ling_eq_check)
as.data.frame(ling_eq_check) %>% 
  pivot_longer(cols = c("MSL", "HS"), names_to = "strategy", values_to = "equilibrium_SB") %>% 
  ggplot(aes(x = fishing_trials, y = equilibrium_SB, fill = strategy)) +
  geom_point()

yield_check = as.data.frame(yield_check) %>% 
  pivot_longer(cols = c("MSL", "HS"), names_to = "strategy", values_to = "equilibrium_SB") %>% 
  ggplot(aes(x = f, y = equilibrium_SB, fill = strategy)) +
  geom_point()


# Quantify balance between fishing and bycatch rate ---------------------------------------------------
f = seq(0.01,0.3, by = 0.02)
b = seq(0.01,0.3, by = 0.02)
harvest.scenarios = as.data.frame(expand_grid(f,b))

fxb_balance_MSL = mapply(get_fxb_balance, f = harvest.scenarios[,1], b = harvest.scenarios[,2], 
                     MoreArgs = list(rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE)
fxb_balance_df = cbind(harvest.scenarios, round(unlist(fxb_balance_MSL))/round(unfished_vals$rockfish_SBeq))

fxb_balance_HS = mapply(get_fxb_balance, f = harvest.scenarios[,1], b = harvest.scenarios[,2], 
                     MoreArgs = list(rockfish = rockfish, lingcod = lingcod, min.selectivity = FALSE), SIMPLIFY = FALSE)
  
fxb_balance_df = cbind(fxb_balance_df, round(unlist(fxb_balance_HS))/round(unfished_vals$rockfish_SBeq))                  
names(fxb_balance_df) = c("f", "b", "MSL", "HS")

MSL_balance = fxb_balance_df %>% 
  select(f, b, MSL) %>% 
  rename(rockfish_recovery = MSL) %>% 
  mutate(strategy = "MSL")

HS_balance = fxb_balance_df %>% 
  select(f, b, HS) %>%
  rename(rockfish_recovery = HS) %>% 
  mutate(strategy = "HS")

final_fxb_balance = rbind(MSL_balance,HS_balance) %>% 
    group_by(strategy, f) %>% 
    arrange(abs(rockfish_recovery - 1)) %>% 
    slice(1)

ggplot(final_fxb_balance, aes(x = f, y = b, color = strategy)) +
  geom_point(position=position_jitter(h=0.001, w=0.001)) +
  theme_classic() +
  labs(title = "Balance between fishing and bycatch rate", subtitle = "to achieve unfished rockfish biomass")

fxb_trials = as.data.frame(cbind(final_fxb_balance$f, final_fxb_balance$b)) %>% 
  mutate(final_fxb_balance$b + 0.05) %>% 
  mutate(final_fxb_balance$b - 0.05) %>% 
  mutate(final_fxb_balance$b + 0.025) %>% 
  mutate(final_fxb_balance$b - 0.025) %>% 
  pivot_longer(2:6, names_to = "b_names", values_to = "b") %>% 
  select(V1,b) %>% 
  rename(f = V1)
  

# deterministic model to get expected outcomes of fishery objectives -----------

trial_MSL = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
            f = 0.1, b = 0.1, a_ij = binned.size.spec.29, handling = 2, times = 1:2, min.selectivity = TRUE)

trial_HS = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                        f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = FALSE)


c(trial_HS$lingcod_SBeq, trial_MSL$lingcod_SBeq)/unfished_vals$lingcod_SBeq # relative equilibrium spawning biomass
c(trial_HS$lingcod_oldprop, trial_MSL$lingcod_oldprop)/unfished_vals$lingcod_oldprop # relative equilibrium plus size proportion lingcod
c(trial_HS$rockfish_oldprop, trial_MSL$rockfish_oldprop)/unfished_vals$rockfish_oldprop # relative equilibrium plus size proportion rockfish
sum(trial_HS$yield[1,(length(trial_HS$yield[1,])-300):length(trial_HS$yield[1,])]) # total biomass yield harvest slot
sum(trial_MSL$yield[1,(length(trial_MSL$yield[1,])-300):length(trial_MSL$yield[1,])]) # total biomass yield minimum size limit

f = seq(0.01, 0.3, length.out = 5)
b = seq(0.01,0.2, length.out = 3)
#min.selectivity = c(TRUE, FALSE)
harvest.scenarios.det = as.data.frame(expand_grid(f,b))

det_run = mcmapply(get_pop_det, f = harvest.scenarios.det[,1], b = harvest.scenarios.det[,2],
       MoreArgs = list(rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE, mc.cores = 4)
save(det_run, file = "results/det_run.Rdata")


# lingcod_det_SBeq = unname(unlist(lapply(det_run, `[`, "lingcod_SBeq"))) # extract lingcod data
rockfish_det_SBeq = unname(unlist(lapply(det_run, `[`, "rockfish_SBeq")))
# lingcod_det_oldprop = unname(unlist(lapply(det_run, `[`, "lingcod_oldprop")))
rockfish_det_oldprop = unname(unlist(lapply(det_run, `[`, "rockfish_oldprop")))
rockfish_biomass = matrix(unname(unlist(lapply(det_run, `[`, "biomass_rockfish"))), nrow = (150+20+300))

det_df = cbind(harvest.scenarios.det, "Rockfish_SB" = rockfish_det_SBeq/unfished_vals$rockfish_SBeq, "Rockfish_oldprop" = rockfish_det_oldprop/unfished_vals$rockfish_oldprop)
save(det_df, file = "results/det_df.Rdata")


# Stochastic simulations ---------------------------------------------------------------------------------
f = seq(0.01, 0.3, length.out = 4)
b = seq(0.01,0.3, length.out = 4)
corr = c(0.8, -0.8)
min.selectivity = c(TRUE, FALSE)
harvest.scenarios = as.data.frame(expand_grid(f,b, corr))

cv_run = mapply(get_pop_stoch, f = harvest.scenarios[,1], b = harvest.scenarios[,2], corr = harvest.scenarios[,3], 
                min.selectivity = harvest.scenarios[,4],
                MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10), SIMPLIFY = FALSE)

cv_run_no_int = mcmapply(get_pop_stoch, f = harvest.scenarios[,1], b = harvest.scenarios[,2], corr = harvest.scenarios[,3], 
                       MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10, a_ij = a_ij.0), SIMPLIFY = FALSE, mc.cores = 4)

cv_run_compare = mcmapply(get_pop_stoch, f = harvest.scenarios[,1], b = harvest.scenarios[,2], corr = harvest.scenarios[,3],
                MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10), SIMPLIFY = FALSE, mc.cores = 4)


# Test just the different correlations
corr_cv_run = mapply(get_pop_stoch, corr = corr,
                    MoreArgs = list(f = 0, b = 0, rockfish = rockfish, lingcod = lingcod, nsim = 100), SIMPLIFY = FALSE)

corr_cv_run_no_int = get_pop_stoch(rockfish, lingcod, nsim = 100, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6,
                                   tf = (30+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                                   f = 0, b = 0, a_ij = a_ij.0, handling = 0.01, times = 1:2, min.selectivity = TRUE)


unname(unlist(lapply(corr_cv_run, `[`, "rockfish_cv")))
unname(unlist(lapply(corr_cv_run, `[`, "rockfish_avg")))/corr_cv_run_no_int$rockfish_avg


# Extract data from list
lingcod_cv = unname(unlist(lapply(cv_run, `[`, "lingcod_cv"))) # extract lingcod data
rockfish_cv = unname(unlist(lapply(cv_run, `[`, "rockfish_cv"))) # extract rockfish data
lingcod_avg = unname(unlist(lapply(cv_run, `[`, "lingcod_avg")))
rockfish_avg = unname(unlist(lapply(cv_run, `[`, "rockfish_avg")))
cv_df = cbind(harvest.scenarios, lingcod_cv, rockfish_cv, lingcod_avg, rockfish_avg) # add lingcod and rockfish data to harvest scenario info

rockfish_cv_no_int = unname(unlist(lapply(cv_run_no_int, `[`, "rockfish_cv")))
rockfish_avg_no_int = unname(unlist(lapply(cv_run_no_int, `[`, "rockfish_avg")))
rockfish_cv_comp = unname(unlist(lapply(cv_run_compare, `[`, "rockfish_cv")))
rockfish_avg_comp = unname(unlist(lapply(cv_run_compare, `[`, "rockfish_avg")))
cv_df_comp = cbind(harvest.scenarios, rockfish_cv_no_int, rockfish_cv_comp, rockfish_avg_no_int, rockfish_avg_comp)

save(cv_df, file = "results/cv_df.Rdata")
save(cv_df_comp, file = "results/cv_df_comp.Rdata")

# get sample time series info
timeseries_run = mapply(get_sample_ts, corr = c(0.8, -0.8),
                        MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10, min.selectivity = TRUE), SIMPLIFY = FALSE)
timeseries_run_no_interaction = mapply(get_sample_ts, corr = c(0.8, -0.8),
                        MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 20, min.selectivity = TRUE, a_ij = a_ij.0), SIMPLIFY = FALSE)
save(timeseries_run, file = "results/timeseries_run.Rdata")
save(timeseries_run_no_interaction, file = "results/timeseries_run_no_interaction.Rdata")

pcorr_consumption_ts = timeseries_run[[1]][["lingcod_consumption"]]
ncorr_consumption_ts = timeseries_run[[2]][["lingcod_consumption"]]

cons_vis = cbind(timeseries_run[[1]][["biomass_ts_1"]][2,1:449], timeseries_run[[1]][["lingcod_consumption"]][1,1:449]/timeseries_run[[1]][["biomass_ts_1"]][2, 1:449])
cons_vis_neg = cbind(timeseries_run[[2]][["biomass_ts_1"]][2,1:449], timeseries_run[[2]][["lingcod_consumption"]][1,1:449]/timeseries_run[[2]][["biomass_ts_1"]][2, 1:449])

plot(cons_vis_neg[,1], cons_vis_neg[,2], type = "p", pch = 19, cex = 0.5)#, xlim = c(7000,69000), ylim = c(0, .000065), xlab = "lingcod biomass", ylab = "consumption rate")
points(cons_vis[,1], cons_vis[,2], type = "p", pch = 19, cex = 0.5, col = "red")
legend(50000, 0.00006, legend=c("negative", "positive"), col=c("black", "red"), lty = 1 , cex=0.8)

# saving 1 simulation of age structure through time
age_ts_pcorr = as.data.frame(rbind(timeseries_run[[1]][["age_ts_3"]][1:20,] + timeseries_run[[1]][["age_ts_3"]][21:40,],
                     timeseries_run[[1]][["age_ts_3"]][-(1:40),]))
age_ts_ncorr = as.data.frame(rbind(timeseries_run[[2]][["age_ts_1"]][1:20,] + timeseries_run[[2]][["age_ts_3"]][21:40,],
                     timeseries_run[[2]][["age_ts_3"]][-(1:40),]))
colnames(age_ts_ncorr) = colnames(age_ts_pcorr) = (time = 1:450)

# Saving spawning biomass info 
ling_SB_pcorr = rbind(timeseries_run[[1]][["SB_ts_1"]][1,], 
                      timeseries_run[[1]][["SB_ts_2"]][1,], 
                      timeseries_run[[1]][["SB_ts_3"]][1,])
ling_SB_ncorr = rbind(timeseries_run[[2]][["SB_ts_1"]][1,], 
                      timeseries_run[[2]][["SB_ts_2"]][1,], 
                      timeseries_run[[2]][["SB_ts_3"]][1,])
rock_SB_pcorr = rbind(timeseries_run[[1]][["SB_ts_1"]][2,], 
                      timeseries_run[[1]][["SB_ts_2"]][2,], 
                      timeseries_run[[1]][["SB_ts_3"]][2,])
rock_SB_ncorr = rbind(timeseries_run[[2]][["SB_ts_1"]][2,], 
                      timeseries_run[[2]][["SB_ts_2"]][2,], 
                      timeseries_run[[2]][["SB_ts_3"]][2,])
rock_SB_no_int = rbind(timeseries_run_no_interaction[[2]][["SB_ts_1"]][2,], 
                       timeseries_run_no_interaction[[2]][["SB_ts_2"]][2,], 
                       timeseries_run_no_interaction[[2]][["SB_ts_3"]][2,])


jpeg("plots/ling_SB_ncorr.jpeg")
matplot(t(ling_SB_ncorr), type = "l", ylab = "spawning biomass", main = "lingcod SB ncorr", ylim = c(0, 51000), lty = 1)
dev.off()
jpeg("plots/ling_SB_pcorr.jpeg")
matplot(t(ling_SB_pcorr), type = "l", ylab = "spawning biomass", main = "lingcod SB pcorr", ylim = c(0, 51000), lty = 1)
dev.off()
jpeg("plots/rock_SB_pcorr.jpeg")
matplot(t(rock_SB_pcorr), type = "l", ylab = "spawning biomass", main = "rockfish SB pcorr", ylim = c(0, 2700), lty = 1)
dev.off()
jpeg("plots/rock_SB_ncorr.jpeg")
matplot(t(rock_SB_ncorr), type = "l", ylab = "spawning biomass", main = "rockfish SB ncorr", ylim = c(0, 2700), lty = 1)
dev.off()
jpeg("plots/rock_SB_no_int.jpeg")
matplot(t(rock_SB_no_int), type = "l", ylab = "spawning biomass", main = "rockfish SB no interaction", ylim = c(0, 5500), lty = 1)
dev.off()

# Saving total biomass info
rock_B_pcorr = rbind(timeseries_run[[1]][["biomass_ts_1"]][2,], 
                      timeseries_run[[1]][["biomass_ts_2"]][2,], 
                      timeseries_run[[1]][["biomass_ts_3"]][2,])
rock_B_ncorr = rbind(timeseries_run[[2]][["biomass_ts_1"]][2,], 
                      timeseries_run[[2]][["biomass_ts_2"]][2,], 
                      timeseries_run[[2]][["biomass_ts_3"]][2,])


jpeg("plots/consumption_pcorr.jpeg")
matplot(t(pcorr_consumption_ts), type = "l", lty = 1, main = "positive correlation", ylab = "Consumption in biomass", ylim = c(1,14))
dev.off()

jpeg("plots/consumption_ncorr.jpeg")
matplot(t(ncorr_consumption_ts), type = "l", lty = 1, main = "negative correlation", ylab = "Consumption in biomass", ylim = c(1,14))
dev.off()

matplot(t(timeseries_run[[1]]$ts_2), type = "l", lty = 1) # positive
matlines(t(timeseries_run[[2]]$ts_2), lty = 2) # negative

# make plots ----------------------------------------------------------------------

library(PNWColors)
lingcod_variability = cv_df %>% 
  select(f, corr, min.selectivity, lingcod_cv) %>% 
  group_by(min.selectivity, f) %>% 
  summarise(mean_cv = mean(lingcod_cv)) %>% 
  ggplot(aes(x = f, y = mean_cv, colour = min.selectivity)) +
  geom_line() +
  theme_classic() +
  labs(colour = "harvest scenario", y = "Mean CV", x = "Fishing Rate", title = "Lingcod long term variability", subtitle = "across harvest scenarios") +
  scale_color_manual(labels = c("Harvest Slot", "Min. Size Limit"), values = pnw_palette("Sunset2", 2)) 

lingcod_avg_det = det_df %>% 
  select(f, min.selectivity, lingcod_det_SBeq) %>% 
  group_by(f, min.selectivity) %>% 
  summarise(det_SBeq = mean(lingcod_det_SBeq)/unfished_vals$lingcod_SBeq)

lingcod_avg = cv_df %>% 
  select(f, corr, min.selectivity, lingcod_avg) %>% 
  group_by(min.selectivity, f) %>% 
  summarise(mean_avg = mean(lingcod_avg)/unfished_vals$lingcod_SBeq) %>% 
  ggplot(aes(x = f, y = mean_avg, colour = min.selectivity)) +
  geom_line() +
  geom_line(data = lingcod_avg_det, aes(x = f, y = det_SBeq, colour = min.selectivity), linetype= "dashed") + 
  theme_classic() +
  labs(colour = "harvest scenario", y = "Spawning Biomass", x = "Fishing Rate", title = "Lingcod long term average SB", subtitle = "across harvest scenarios") +
  scale_color_manual(labels = c("Harvest Slot", "Min. Size Limit"), values = pnw_palette("Sunset2", 4)) 


rockfish_pcorr_msl_cv = cv_df %>% 
  filter(corr == "positive") %>% 
  filter(min.selectivity == "MSL") %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(18, 2), trans = 'reverse') +
  labs(title = "positive correlation MSL")

rockfish_ncorr_msl_cv = cv_df %>% 
  filter(corr == "negative") %>% 
  filter(min.selectivity == "MSL") %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(18, 2), trans = 'reverse') +
  labs(title = "negative correlation MSL")


rockfish_pcorr_HS_cv = cv_df %>% 
  filter(corr == "positive") %>% 
  filter(min.selectivity == "HS") %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(18, 2), trans = 'reverse') +
  labs(title = "positive correlation HS")

rockfish_ncorr_HS_cv = cv_df %>% 
  filter(corr == "negative") %>% 
  filter(min.selectivity == "HS") %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(18, 2), trans = 'reverse') +
  labs(title = "negative correlation HS")
  

# Looking at how rockfish variability changes with lingcod variability 
library("viridis")
cv_df[cv_df == -0.8] <- "negative" 
cv_df[cv_df == 0.8] <- "positive" 
cv_df[cv_df == TRUE] <- "MSL" 
cv_df[cv_df == FALSE] <- "HS" 
cv_df %>% 
  select(lingcod_cv, rockfish_cv, corr, b, min.selectivity) %>%
  ggplot(aes(x = lingcod_cv, y = rockfish_cv, colour = corr, group=min.selectivity)) +
  geom_point(aes(shape=min.selectivity)) +
  theme_classic() +
  scale_color_manual(values = pnw_palette("Sunset2", 2)) 

cv_df %>% 
  select(lingcod_cv, rockfish_cv, corr, min.selectivity) %>%
  ggplot(aes(x = lingcod_cv, y = rockfish_cv, colour = min.selectivity)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = pnw_palette("Sunset2", 2)) 

cv_df %>% 
  select(f, b, rockfish_avg, corr, min.selectivity) %>%
  ggplot(aes(x = f, y = b, fill = rockfish_avg/unfished_vals$rockfish_SBeq)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient2(midpoint = 1, low = "red", mid = "white", high = "blue") +
  facet_wrap(~corr + min.selectivity) +
  labs(title = "Average Rockfish Spawning Biomass", fill = "Relative SB")

cv_df %>% 
  select(f, b, rockfish_cv, corr, min.selectivity) %>%
  filter(corr == "positive") %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = "plasma", trans = 'reverse') +
  facet_wrap(~min.selectivity) +
  labs(title = "Average Rockfish CV", fill = "Mean CV")

cv_df %>% 
  select(f, b, rockfish_cv, corr, min.selectivity) %>%
  filter(corr == "negative") %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(trans = "reverse") +
  facet_wrap(~min.selectivity) +
  labs(title = "Average Rockfish CV", fill = "Mean CV")


# how rockfish variability and SB compare to no interaction
cv_df_comp[cv_df_comp == -0.8] <- "negative" 
cv_df_comp[cv_df_comp == 0.8] <- "positive" 

jpeg("plots/interaction_impact_cv.jpeg")
as.data.frame(cv_df_comp) %>% 
  select(f, b, corr, rockfish_cv_no_int, rockfish_cv_comp) %>%
  rename(no.interaction = rockfish_cv_no_int) %>% 
  rename(interaction = rockfish_cv_comp) %>% 
  ggplot(aes(x = f, y = b, fill = interaction/no.interaction)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient2(midpoint = 1, low = "red", mid = "white", high = "blue") +
  facet_wrap(~corr) +
  labs(title = "CV relative to no interaction", fill = "Relative CV")
dev.off()

jpeg("plots/interaction_impact_SB.jpeg")
as.data.frame(cv_df_comp) %>% 
  select(f, b, corr, rockfish_avg_no_int, rockfish_avg_comp) %>%
  rename(no.interaction = rockfish_avg_no_int) %>% 
  rename(interaction = rockfish_avg_comp) %>% 
  ggplot(aes(x = f, y = b, fill = interaction/no.interaction)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient2(midpoint = 1, low = "red", mid = "white", high = "blue") +
  facet_wrap(~corr) +
  labs(title = "SB relative to no interaction", fill = "Relative SB")
dev.off()

# Age structure through time

jpeg("plots/lingcod_age_ts_pcorr.jpeg")
as.data.frame(age_ts_pcorr) %>% 
  mutate(age = c(1:20,1:65)) %>% 
  mutate(species = c(rep("lingcod", lingcod$nage),rep("yelloweye", rockfish$nage))) %>% 
  gather(time, N, 1:450) %>% 
  filter(species == "lingcod") %>% 
  filter(time %in% 350:450) %>% 
  ggplot(aes(x = time, y = age, fill = N)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(trans = "reverse", limits = c(12700, 1)) +
  theme(axis.text.x = element_text(angle=45)) +
  lims(y = c(0, 65)) +
  labs(title = "lingcod positive corr")
dev.off()

jpeg("plots/rockfish_age_ts_pcorr_3.jpeg")
as.data.frame(age_ts_pcorr) %>% 
  mutate(age = c(1:20,1:65)) %>% 
  mutate(species = c(rep("lingcod", lingcod$nage),rep("yelloweye", rockfish$nage))) %>% 
  gather(time, N, 1:450) %>% 
  filter(species == "yelloweye") %>% 
  filter(time %in% 350:450) %>% 
  ggplot(aes(x = time, y = age, fill = N)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(low = "yellow", high = "red", limits = c(0, 680)) +
  theme(axis.text.x = element_text(angle=45)) +
  lims(y = c(1,65)) +
  labs(title = "Rockfish positive corr")
dev.off()

jpeg("plots/lingcod_age_ts_ncorr.jpeg")
as.data.frame(age_ts_ncorr) %>% 
  mutate(age = c(1:20,1:65)) %>% 
  mutate(species = c(rep("lingcod", lingcod$nage),rep("yelloweye", rockfish$nage))) %>% 
  gather(time, N, 1:450) %>% 
  filter(species == "lingcod") %>% 
  filter(time %in% 350:450) %>% 
  ggplot(aes(x = time, y = age, fill = N)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(trans = "reverse", limits = c(12700, 1)) +
  theme(axis.text.x = element_text(angle=45)) +
  lims(y = c(0, 65)) +
  labs(title = "lingcod negative corr")
dev.off()

jpeg("plots/rockfish_age_ts_ncorr_3.jpeg")
as.data.frame(age_ts_ncorr) %>% 
  mutate(age = c(1:20,1:65)) %>% 
  mutate(species = c(rep("lingcod", lingcod$nage),rep("yelloweye", rockfish$nage))) %>% 
  gather(time, N, 1:450) %>% 
  filter(species == "yelloweye") %>% 
  filter(time %in% 350:450) %>% 
  ggplot(aes(x = time, y = age, fill = N)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(low = "yellow", high = "red", limits = c(0, 680)) +
  theme(axis.text.x = element_text(angle=45)) +
  lims(y = c(1,65)) +
  labs(title = "Rockfish negative corr")
dev.off()

# check the correlation between lingcod and rockfish age-1s
as.data.frame(age_ts_ncorr) %>% 
  mutate(age = c(1:20,1:65)) %>% 
  mutate(species = c(rep("lingcod", lingcod$nage),rep("yelloweye", rockfish$nage))) %>% 
  gather(time, N, 1:450) %>% 
  filter(age == 1) %>% 
  filter(time %in% 425:450) %>% 
  ggplot(aes(x = time, y = N, col = species)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45)) +
  labs(title = "Age-1 negative corr")


# Biomass and Spawning Biomass through time

# lingcod




# Consumption variation ----------------------------------------------------------







