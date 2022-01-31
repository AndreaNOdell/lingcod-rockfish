library(miceadds)
source("scripts/scratch_script.R")
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish) # with interaction

# Unfished equilibrium ---------------------------------------------------------------------------------
unfished_vals = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
                            f = 0, b = 0, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough
lingcod_unfished_equil = unfished_vals$lingcod_SBeq
rockfish_unfished_equil = unfished_vals$rockfish_SBeq
lingcod_unfished_oldprop = lingcod_oldprop
rockfish_unfished_oldprop = rockfish_oldprop
lingcod_unfished_SAD = lingcod_SAD
rockfish_unfished_SAD = rockfish_SAD

# Quantify balance between fishing and bycatch rate ---------------------------------------------------

f = seq(0.01,0.3, by = 0.02)
b = seq(0,0.3, by = 0.02)
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
  filter(MSL > 0.97) %>% 
  filter(MSL < 1.02) %>% 
  rename(rockfish_recovery = MSL) %>% 
  mutate(strategy = "MSL")

HS_balance = fxb_balance_df %>% 
  select(f, b, HS) %>% 
  filter(HS > 0.97) %>% 
  filter(HS < 1.02) %>% 
  rename(rockfish_recovery = HS) %>% 
  mutate(strategy = "HS")

final_fxb_balance = rbind(MSL_balance,HS_balance)

ggplot(final_fxb_balance, aes(x = f, y = b, color = strategy)) +
  geom_point() +
  theme_classic() +
  labs(title = "Balance between fishing and bycatch rate", subtitle = "to achieve unfished rockfish biomass")

# deterministic model to get expected outcomes of fishery objectives 

trial_MSL = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
            f = 0.1, b = 0.15, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE)

trial_HS = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                        f = 0.1, b = 0.13, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = FALSE)


c(trial_HS$lingcod_Beq, trial_MSL$lingcod_Beq)/unfished_vals$lingcod_Beq
c(trial_HS$lingcod_oldprop, trial_MSL$lingcod_oldprop)/unfished_vals$lingcod_oldprop
c(trial_HS$rockfish_oldprop, trial_MSL$rockfish_oldprop)/unfished_vals$rockfish_oldprop
sum(trial_HS$yield[2,(length(trial_HS$yield[2,])-300):length(trial_HS$yield[2,])])
sum(trial_MSL$yield[2,(length(trial_MSL$yield[2,])-300):length(trial_MSL$yield[2,])])

# Stochastic simulations ---------------------------------------------------------------------------------

stoch_trial_MSL = get_pop_stoch(rockfish, lingcod, init.l = 200, init.r = 70, nsim = 4, 
                                corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6,
                                tf = (150+20+300), mpa.yr = 20, hist.f = 0.5, hist.by = 0.5, 
                                f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
                                min.selectivity = TRUE)

f = seq(0.01, 0.3, length.out = 10)
b = seq(0.01,0.3, length.out = 10)
corr = c(0.8, -0.8)
min.selectivity = c(TRUE, FALSE)
harvest.scenarios = as.data.frame(expand_grid(f,b, corr, min.selectivity))

cv_run = mapply(get_pop_stoch, f = harvest.scenarios[,1], b = harvest.scenarios[,2], corr = harvest.scenarios[,3], 
                min.selectivity = harvest.scenarios[,4],
                MoreArgs = list(rockfish = rockfish, lingcod = lingcod, nsim = 5), SIMPLIFY = FALSE)

# Extract data from list
lingcod_cv = unname(unlist(lapply(cv_run, `[`, "lingcod_cv"))) # extract lingcod data
rockfish_cv = unname(unlist(lapply(cv_run, `[`, "rockfish_cv"))) # extract rockfish data
cv_df = cbind(harvest.scenarios, lingcod_cv, rockfish_cv) # add lingcod and rockfish data to harvest scenario info

save(cv_df, file = "results/cv_df.Rdata")


# make plots ----------------------------------------------------------------------
lingcod_pcorr_msl_cv = cv_df %>% 
  filter(corr == 0.8) %>% 
  filter(min.selectivity == TRUE) %>% 
  select(f,b,lingcod_cv) %>% 
  ggplot(aes(x = f, y = b, fill = lingcod_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(25,12), trans = 'reverse') +
  labs(title = "positive correlation MSL")

lingcod_ncorr_msl_cv = cv_df %>% 
  filter(corr == -0.8) %>% 
  filter(min.selectivity == TRUE) %>% 
  select(f,b,lingcod_cv) %>% 
ggplot(aes(x = f, y = b, fill = lingcod_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(25,12), trans = 'reverse') +
  labs(title = "negative correlation MSL")


rockfish_pcorr_msl_cv = cv_df %>% 
  filter(corr == 0.8) %>% 
  filter(min.selectivity == TRUE) %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(17, 2), trans = 'reverse') +
  labs(title = "positive correlation MSL")

rockfish_ncorr_msl_cv = cv_df %>% 
  filter(corr == -0.8) %>% 
  filter(min.selectivity == TRUE) %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(17, 2), trans = 'reverse') +
  labs(title = "negative correlation MSL")

lingcod_pcorr_HS_cv = cv_df %>% 
  filter(corr == 0.8) %>% 
  filter(min.selectivity == FALSE) %>% 
  select(f,b,lingcod_cv) %>% 
  ggplot(aes(x = f, y = b, fill = lingcod_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(25,12), trans = 'reverse') +
  labs(title = "positive correlation HS")

lingcod_ncorr_HS_cv = cv_df %>% 
  filter(corr == -0.8) %>% 
  filter(min.selectivity == FALSE) %>% 
  select(f,b,lingcod_cv) %>% 
  ggplot(aes(x = f, y = b, fill = lingcod_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(25,12), trans = 'reverse') +
  labs(title = "negative correlation HS")


rockfish_pcorr_HS_cv = cv_df %>% 
  filter(corr == 0.8) %>% 
  filter(min.selectivity == FALSE) %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(17, 2), trans = 'reverse') +
  labs(title = "positive correlation HS")

rockfish_ncorr_HS_cv = cv_df %>% 
  filter(corr == -0.8) %>% 
  filter(min.selectivity == FALSE) %>% 
  select(f,b,rockfish_cv) %>% 
  ggplot(aes(x = f, y = b, fill = rockfish_cv)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient(limits = c(17, 2), trans = 'reverse') +
  labs(title = "negative correlation HS")
  
  
  
  
  
  
  

