source("scripts/scratch_script.R")
library(miceadds)
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish) # with interaction

# Unfished equilibrium ---------------------------------------------------------------------------------
unfished_vals = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
                            f = 0, b = 0, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough

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
            f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE)

trial_HS = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                        f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = FALSE)


c(trial_HS$lingcod_SBeq, trial_MSL$lingcod_SBeq)/unfished_vals$lingcod_SBeq # relative equilibrium spawning biomass
c(trial_HS$lingcod_oldprop, trial_MSL$lingcod_oldprop)/unfished_vals$lingcod_oldprop # relative equilibrium plus size proportion lingcod
c(trial_HS$rockfish_oldprop, trial_MSL$rockfish_oldprop)/unfished_vals$rockfish_oldprop # relative equilibrium plus size proportion rockfish
sum(trial_HS$yield[1,(length(trial_HS$yield[1,])-300):length(trial_HS$yield[1,])]) # total biomass yield harvest slot
sum(trial_MSL$yield[1,(length(trial_MSL$yield[1,])-300):length(trial_MSL$yield[1,])]) # total biomass yield minimum size limit

f = seq(0.01, 0.3, length.out = 10)
b = seq(0.01,0.3, length.out = 10)
min.selectivity = c(TRUE, FALSE)
harvest.scenarios.det = as.data.frame(expand_grid(f,b, min.selectivity))

det_run = mapply(get_pop_det, f = harvest.scenarios.det[,1], b = harvest.scenarios.det[,2], 
       min.selectivity = harvest.scenarios.det[,3],
       MoreArgs = list(rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE)
save(det_run, file = "results/det_run.Rdata")


lingcod_det_SBeq = unname(unlist(lapply(det_run, `[`, "lingcod_SBeq"))) # extract lingcod data
rockfish_det_SBeq = unname(unlist(lapply(det_run, `[`, "rockfish_SBeq")))
lingcod_det_oldprop = unname(unlist(lapply(det_run, `[`, "lingcod_oldprop")))
rockfish_det_oldprop = unname(unlist(lapply(det_run, `[`, "rockfish_oldprop")))

det_df = cbind(harvest.scenarios.det, lingcod_det_SBeq, rockfish_det_SBeq, lingcod_det_oldprop, rockfish_det_oldprop)


# Stochastic simulations ---------------------------------------------------------------------------------
f = seq(0.01, 0.3, length.out = 10)
b = seq(0.01,0.3, length.out = 10)
corr = c(0.8, -0.8)
min.selectivity = c(TRUE, FALSE)
harvest.scenarios = as.data.frame(expand_grid(f,b, corr, min.selectivity))

cv_run = mapply(get_pop_stoch, f = harvest.scenarios[,1], b = harvest.scenarios[,2], corr = harvest.scenarios[,3], 
                min.selectivity = harvest.scenarios[,4],
                MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 10), SIMPLIFY = FALSE)

# Extract data from list
lingcod_cv = unname(unlist(lapply(cv_run, `[`, "lingcod_cv"))) # extract lingcod data
rockfish_cv = unname(unlist(lapply(cv_run, `[`, "rockfish_cv"))) # extract rockfish data
lingcod_avg = unname(unlist(lapply(cv_run, `[`, "lingcod_avg")))
rockfish_avg = unname(unlist(lapply(cv_run, `[`, "rockfish_avg")))
cv_df = cbind(harvest.scenarios, lingcod_cv, rockfish_cv, lingcod_avg, rockfish_avg) # add lingcod and rockfish data to harvest scenario info

save(cv_df, file = "results/cv_df.Rdata")


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





