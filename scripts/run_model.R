library(tidyverse)
library(miceadds)
library(abind)
library(magrittr)
library(viridis)
source("scripts/impulse_dynamics.R") # load in population dynamic function

# Key parameters to adjust
burn.in = 150
mpa.yr = 20
tf = burn.in+mpa.yr+300
cv = 0.6
simulations = 5

load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish) # with interaction
a_ij.0 = binned.size.spec.29 * 0 # without interaction

# Check alternative stable states ----------------------------------------------------
init.r = seq(30, 90, length.out = 7)
init.l = seq(200, 800, length.out = 7)
init.vals = expand.grid(init.l, init.r)

# empty matrix
alt.states.check = matrix(NA, nrow(init.vals), 2)

# run the model with unique combinations of initial starting values
for (r in 1:nrow(init.vals)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = init.vals[r,1], init.r = init.vals[r,2], corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = FALSE)
  alt.states.check[r,1] = lingcod_pop[1,470]
  alt.states.check[r,2] = rockfish_pop[1,470]
}

# create dataframe
alt.states.check = cbind(init.vals, alt.states.check)
colnames(alt.states.check) = c("init.l", "init.r", "lingcod", "rockfish")
alt.states.check = pivot_longer(alt.states.check, !c(init.l,init.r), names_to = "species", values_to = "abundance")

jpeg('plots/alt.states.check.jpg')
ggplot(alt.states.check, aes(init.l, init.r)) +
  geom_tile(aes(fill = abundance)) +
  scale_fill_continuous(trans = 'reverse') +
  facet_wrap(~species) +
  theme_classic()
dev.off()


# Run deterministic model with varying size-selective predation --------------------------------------------

load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
load.Rdata( filename="cleaned_data/binned.size.spec.27.Rdata", "binned.size.spec.27" ) 
load.Rdata( filename="cleaned_data/binned.size.spec.31.Rdata", "binned.size.spec.31" ) 
load.Rdata( filename="cleaned_data/binned.size.spec.33.Rdata", "binned.size.spec.33" ) 
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.27 = binned.size.spec.27 %*% diag(diet.frac.rockfish) 
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish)
a_ij.31 = binned.size.spec.31 %*% diag(diet.frac.rockfish)
a_ij.33 = binned.size.spec.33 %*% diag(diet.frac.rockfish)

a_ij_trials = abind(a_ij.27, a_ij.33, along = 3) 
vary_dietspec_rockfish = matrix(NA, nrow = tf, ncol = 2)

for(t in 1:2) {
get_pop_n(rockfish, lingcod, nsim = 1, init.l = 500, init.r = 50, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
          tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = 0.1, b = 0.1, a_ij = a_ij_trials[,,t], handling = 0.01, times = 1:2, 
          stochastic = FALSE)
  vary_dietspec_rockfish[,t] = rockfish_pop
}

# what is the percent difference between the equilibriums of the size spectra scenarios
vary_dietspec_percdiff = (vary_dietspec_rockfish[tf,2] - vary_dietspec_rockfish[tf,1]) / vary_dietspec_rockfish[tf,1] * 100

# Pivot dataset into long format
vary_dietspec_rockfish = as.data.frame(vary_dietspec_rockfish) %>% 
  set_colnames(c("lower", "higher")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "size.spec", values_to = "abundance")

# Plot
jpeg('plots/det.vary.dietspec.jpg')
ggplot(vary_dietspec_rockfish, aes(x=time, y = abundance, color = size.spec)) + 
  geom_line() + 
  labs(y = "Abundance", x = "Years", title = "Varying Size-Spectra") +
  theme_classic() +
  theme(legend.position = c(0.7,0.7)) +
  scale_color_discrete(name = "Slope of 95th quantile", labels = c("0.33", "0.27"))
dev.off()

# Deterministic run - the qualitative results do not change (similar trajectory and 
# slightly different equilibrium)
# I looked just at the rockfish dynamics because the lingcod dynamics are not impacted by the
# varying consumption of rockfish (since their consumption of rockfish only affects rockfish abundance)


# Run deterministic model with varying post-MPA fishing pressure ------------------------------------
fishing_trials = seq(0, 0.5, by = 0.1)
vary_fishing_rock_det = matrix(NA, nrow = tf, ncol = length(fishing_trials))
vary_fishing_ling_det = matrix(NA, nrow = tf, ncol = length(fishing_trials))

for (f in 1:length(fishing_trials)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = 200, init.r = 70, corr = 0.8, autocorr = c(0.23,0.23), cv =log(0.6), mn = log(1-(0.5*cv^2)), 
            tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = fishing_trials[f], b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = FALSE)
  vary_fishing_rock_det[,f] = rockfish_pop
  vary_fishing_ling_det[,f] = lingcod_pop
}

vary_fishing_rock_det_long = as.data.frame(vary_fishing_rock_det) %>% 
  set_colnames(c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "Fishing.Pressure", values_to = "Abundance")

vary_fishing_ling_det_long = as.data.frame(vary_fishing_ling_det) %>% 
  set_colnames(c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "Fishing.Pressure", values_to = "Abundance")

# Plot
jpeg('plots/det.vary.fishing.rockf.jpg')
ggplot(vary_fishing_rock_det_long, aes(x=time, y = Abundance, color = Fishing.Pressure)) + 
  geom_line() + 
  geom_point(aes(x = burn.in, y = 2000), size = 0.75) +
  #geom_point(aes(x = 150+mpa.yr, y = 2000), size = 0.75) +
  geom_segment(aes(x = burn.in, xend = (burn.in+mpa.yr), y = 2000, yend = 2000)) +
  labs(y = "Abundance", x = "Years", title = "Rockfish Time-Series", 
       subtitle = "Varying Fishing Pressure") +
  #lims(y = c(1500, 3500)) +
  theme_classic() +
  theme(legend.position = "right") + 
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1)
dev.off()

jpeg('plots/det.vary.fishing.ling.jpg')
ggplot(vary_fishing_ling_det_long, aes(x=time, y = Abundance, color = Fishing.Pressure)) + 
  geom_line() +
  labs(y = "Abundance", x = "Years", title = "Lingcod Time-Series", 
       subtitle = "Varying Fishing Pressure") +
  theme_classic() +
  theme(legend.position = "right") + 
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1)
dev.off()



# Run deterministic model with varying handling time ---------------------------------------
h_trials = c(0.01, 0.05, 0.15, 0.25)
vary_handling_rock_det = matrix(NA, nrow = tf, ncol = length(h_trials))

for (h in 1:length(h_trials)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = 500, init.r = 50, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = 0.2, b = 0.1, a_ij = a_ij.29, handling = h_trials[h], times = 1:2, 
            stochastic = FALSE)
  vary_handling_rock_det[,h] = rockfish_pop
}

jpeg('plots/det.vary.handling.equil.rockf.jpg')
plot(vary_handling_rock_det[tf,], type = "p", ylab = "Equilibrium Abundance", xlab = "handling time")
dev.off()

vary_handling_percdiff = (vary_handling_rock_det[tf,length(h_trials)] - vary_handling_rock_det[tf,1]) / vary_handling_rock_det[tf,1] * 100
    # 2565 -> 2571 = 0.23% difference between lowest and highest handling scenarios

vary_handling_rock_det_long = as.data.frame(vary_handling_rock_det) %>% 
  set_colnames(c("0.01", "0.05", "0.15", "0.25")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "Handling.time", values_to = "Abundance")

jpeg('plots/det.vary.handling.rockf.jpg')
ggplot(vary_handling_rock_det_long, aes(x=time, y = Abundance, color = Handling.time)) + 
  geom_line() + 
  labs(y = "Abundance", x = "Years", title = "Rockfish Time-Series", 
       subtitle = "Varying Handing Time") +
  theme_classic() +
  theme(legend.position = "right") + 
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1)
dev.off()
# This dynamics of this model are not sensitive to handling time.



# Run deterministic model with varying pre-MPA fishing pressure -------------------------------------
hist.f_trials = seq(0.1, 0.5, by = 0.1)
vary_hist.f_rock_det = matrix(NA, nrow = tf, ncol = length(hist.f_trials))
vary_hist.f_ling_det = matrix(NA, nrow = tf, ncol = length(hist.f_trials))

for (hf in 1:length(hist.f_trials)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = 700, init.r = 70, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = tf, mpa.yr = mpa.yr, hist.f = hist.f_trials[hf], hist.by = 0.5, f = 0.2, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = FALSE)
  vary_hist.f_rock_det[,hf] = rockfish_pop
  vary_hist.f_ling_det[,hf] = lingcod_pop
}

vary_hist.f_rock_det_long = as.data.frame(vary_hist.f_rock_det) %>% 
  set_colnames(c("0.1", "0.2","0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "Hist.Fishing.Pressure", values_to = "Abundance")

vary_hist.f_ling_det_long = as.data.frame(vary_hist.f_ling_det) %>% 
  set_colnames(c("0.1", "0.2","0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:tf) %>% 
  pivot_longer(!time, names_to = "Hist.Fishing.Pressure", values_to = "Abundance")

ggplot(vary_hist.f_rock_det_long, aes(x=time, y = Abundance, color = Hist.Fishing.Pressure)) + 
  geom_line() + 
  #geom_point(aes(x = 150, y = 2000), size = 0.75) +
  #geom_point(aes(x = 150+mpa.yr, y = 2000), size = 0.75) +
  #geom_segment(aes(x = 150, xend = (150+mpa.yr), y = 2000, yend = 2000)) +
  labs(y = "Abundance", x = "Years", title = "Rockfish Time-Series", 
       subtitle = "Varying Historical Fishing Pressure") +
  #lims(y = c(2000, 3500)) +
  theme_classic() +
  theme(legend.position = "right") + 
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1)

ggplot(vary_hist.f_ling_det_long, aes(x=time, y = Abundance, color = Hist.Fishing.Pressure)) + 
  geom_line() + 
  #geom_point(aes(x = 150, y = 2000), size = 0.75) +
  #geom_point(aes(x = 150+mpa.yr, y = 2000), size = 0.75) +
  #geom_segment(aes(x = 150, xend = (150+mpa.yr), y = 2000, yend = 2000)) +
  labs(y = "Abundance", x = "Years", title = "Lingcod Time-Series", 
       subtitle = "Varying Historical Fishing Pressure") +
  #lims(y = c(2000, 3500)) +
  theme_classic() +
  theme(legend.position = "right") + 
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1)


# Compare deterministic models with and without species interaction (no lingcod fishery) ---------------------------
# Look at the effect of the closure when you consider single-species approach to multi-species approach

interaction_effect_rockfish = matrix(NA, nrow = tf, ncol = 2)
colnames(interaction_effect_rockfish) = c("with_interaction", "without_interaction")


get_pop_n(rockfish, lingcod, nsim = 1, init.l = 200, init.r = 70, corr = -0.8, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1 - (0.5*cv^2)), 
          tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = 0, b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
          stochastic = FALSE) # run deterministic model WITH interaction
interaction_effect_rockfish[ ,1] = rockfish_pop # plug in rockfish time-series output into matrix

get_pop_n(rockfish, lingcod, nsim = 1, init.l = 200, init.r = 70, corr = -0.8, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1 - (0.5*cv^2)), 
          tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = 0, b = 0.1, a_ij = a_ij.0, handling = 0.01, times = 1:2, 
          stochastic = FALSE) # run deterministic model WITHOUT interaction
interaction_effect_rockfish[ ,2] = rockfish_pop # plug in rockfish time-series output into matrix

interaction_effect_rockfish = as.data.frame(interaction_effect_rockfish) %>% 
  mutate(tf = 1:tf) %>% 
  pivot_longer(!tf, names_to = "interaction", values_to = "abundance")

jpeg('plots/rock.recovery.interaction.effect.jpg')
ggplot(interaction_effect_rockfish, aes(x = tf, y = abundance, color = interaction)) +
  geom_line() +
  theme_classic() +
  geom_vline(xintercept = (burn.in), linetype = "dotted") +
  labs(title = "Recovery expectations following spatial closure", subtitle = "from single vs multi-species perspective", x = "year") +
  scale_color_discrete(name = "Approach", labels = c("multi-species", "single species"))
dev.off()


# Check best time frame to subset for calculating CV --------------------------------------------------------------
fishing_trials = seq(0, 0.1, by = 0.1)
t_removed = seq(150, 400, by = 50)
rockfish_cv_moving_window = matrix(NA, nrow = length(fishing_trials), ncol = length(t_removed)) # nrow is number of simulations
lingcod_cv_moving_window = matrix(NA, nrow = length(fishing_trials), ncol = length(t_removed)) # nrow is number of simulations

for(f in 1:length(fishing_trials)) {
    # Run the function over each fishing pressure scenario 
    get_pop_n(rockfish, lingcod, nsim = simulations, init.l = 200, init.r = 70, corr = -0.8, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1 - (0.5*cv^2)), 
              tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = fishing_trials[f], b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
              stochastic = TRUE)
  for(cv in 1:length(t_removed)) {
    # Years to toss
    rockfish_pop_subset = rockfish_pop[, -(1:t_removed[cv])]
    lingcod_pop_subset = lingcod_pop[, -(1:t_removed[cv])]
    
    # Calculate mean CV and print into dataset
    rockfish_cv_moving_window[f, cv] = mean(apply(rockfish_pop_subset, 1, function(x) sd(x) / mean(x) * 100))
    lingcod_cv_moving_window[f, cv] = mean(apply(lingcod_pop_subset, 1, function(x) sd(x) / mean(x) * 100))
  }
}

# Run stochastic model with varying fishing pressure -----------------------------------------
fishing_trials = seq(0, 0.5, by = 0.1)
vary_fishing_rock_stoch_summary = array(NA, dim = c(tf, 4, length(fishing_trials)))
vary_fishing_ling_stoch_summary = array(NA, dim = c(tf, 4, length(fishing_trials)))
vary_fishing_rock_stock_cv = matrix(NA, nrow = simulations, ncol = length(fishing_trials)) # nrow is number of simulations
vary_fishing_ling_stock_cv = matrix(NA, nrow = simulations, ncol = length(fishing_trials))# nrow is number of simulations
fishery_yields_vary_f = array(NA, dim = c(simulations, tf, length(fishing_trials)))

for(f in 1:length(fishing_trials)) {
 # Run the function over each fishing pressure scenario 
  get_pop_n(rockfish, lingcod, nsim = simulations, init.l = 200, init.r = 70, corr = -0.8, autocorr = c(0.23,0.23), cv = log(cv), mn = log(1 - (0.5*cv^2)), 
            tf = tf, mpa.yr = mpa.yr, hist.f = 0.5, hist.by = 0.5, f = fishing_trials[f], b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = TRUE)
  # Fill in array with data
  vary_fishing_rock_stoch_summary[,1,f] = 1:tf
  vary_fishing_ling_stoch_summary[,1,f] = 1:tf
  vary_fishing_rock_stoch_summary[,2,f] = colMeans(rockfish_pop)
  vary_fishing_ling_stoch_summary[,2,f] = colMeans(lingcod_pop)
  vary_fishing_rock_stoch_summary[,c(3,4),f] = t(apply(rockfish_pop, 2, function(x) quantile(x, c(0.1, 0.9))))
  vary_fishing_ling_stoch_summary[,c(3,4),f] = t(apply(lingcod_pop, 2, function(x) quantile(x, c(0.1, 0.9))))
  
 # Toss out the initial 300 years as a burn-in phase
  rockfish_pop_subset = rockfish_pop[, -(1:250)]
  lingcod_pop_subset = lingcod_pop[, -(1:250)]
 # Matrix as output - but in the forloop, all we want is the coefficient of variation 
    # (and maybe mean)
  # Which will look like matrix with fishing_trials columns and nsims rows for CV 
  # and similarly for mean
  vary_fishing_rock_stock_cv[, f] = apply(rockfish_pop_subset, 1, function(x) sd(x) / mean(x) * 100)
  vary_fishing_ling_stock_cv[, f] = apply(lingcod_pop_subset, 1, function(x) sd(x) / mean(x) * 100)
  
  # Add fishery yield data to array
  fishery_yields_vary_f[,,f] = fishery_yields
}

colnames(vary_fishing_rock_stoch_summary) = colnames(vary_fishing_ling_stoch_summary) = c("time", "mean", "lower", "upper")

# saving the data - adjust the names to fit the correct run (i.e. whether is positive and negative correlation)
# vary_fishing_rock_stoch_summary_ncorr = vary_fishing_rock_stoch_summary
# vary_fishing_ling_stoch_summary_ncorr = vary_fishing_ling_stoch_summary
# vary_fishing_rock_stock_cv_ncorr = vary_fishing_rock_stock_cv
# vary_fishing_ling_stock_cv_ncorr = vary_fishing_ling_stock_cv
# vary_fishing_fishery_yield_ncorr = fishery_yields_vary_f
# save(vary_fishing_rock_stock_cv_ncorr, file = "cleaned_data/vary_fishing_rock_stock_cv_ncorr.Rdata")  
# save(vary_fishing_ling_stock_cv_ncorr, file = "cleaned_data/vary_fishing_ling_stock_cv_ncorr.Rdata")   
# save(vary_fishing_rock_stoch_summary_ncorr, file = "cleaned_data/vary_fishing_rock_stoch_summary_ncorr.Rdata")  
# save(vary_fishing_ling_stoch_summary_ncorr, file = "cleaned_data/vary_fishing_ling_stoch_summary_ncorr.Rdata")  
# save(vary_fishing_fishery_yield_ncorr, file = "cleaned_data/vary_fishing_fishery_yield_ncorr")
 
# rockfish
    # CV plots
jpeg('plots/cv.stoch.varyfishpress.rock.jpg')  
plot(fishing_trials,colMeans(vary_fishing_rock_stock_cv_pcorr), type = "l", ylab = "mean CV", xlab = "Fishing Pressure", main = "Rockfish", ylim = c(10, 20))  
lines(fishing_trials,colMeans(vary_fishing_rock_stock_cv_ncorr), col = "blue")
legend(0, 20, legend=c("positive", "negative"),
       col=c("black", "blue"), lty=1)
dev.off()

    # ribbon plots
jpeg('plots/stoch.det.varyfishpress.rock.jpg') 
ggplot(as.data.frame(vary_fishing_rock_stoch_summary_ncorr[,,2])) +
  geom_line(aes(x = time, y = mean, color = "mean")) +
  labs(title = "rockfish") +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = time, fill = "mid 80th quant"), alpha = 0.3) + # Central 80% of randomness in abundance
  theme_classic() +
  geom_line(subset(vary_fishing_rock_det_long, Fishing.Pressure == 0.1), mapping = aes(x=time, y = Abundance))
dev.off()

# lingcod
jpeg('plots/cv.stoch.varyfishpress.ling.jpg') 
plot(fishing_trials,colMeans(vary_fishing_ling_stock_cv_pcorr), type = "l", ylab = "mean CV", xlab = "Fishing Pressure", main = "Lingcod")  
lines(fishing_trials,colMeans(vary_fishing_ling_stock_cv_ncorr), col = "blue")
legend(0.33, 23.6, legend=c("positive", "negative"),
       col=c("black", "blue"), lty=1)
dev.off()

    # ribbon plot
jpeg('plots/stoch.det.varyfishpress.ling.jpg') 
ggplot(as.data.frame(vary_fishing_ling_stoch_summary_ncorr[,,2])) +
  geom_line(aes(x = time, y = mean, color = "mean")) +
  labs(title = "lingcod") +
  geom_ribbon(aes(ymin = lower, ymax = upper, x = time, fill = "mid 80th quant"), alpha = 0.3) + # Central 80% of randomness in abundance
  theme_classic() +
  geom_line(subset(vary_fishing_ling_det_long, Fishing.Pressure == 0.1), mapping = aes(x=time, y = Abundance))
dev.off()  

# Determine the change in CV over fishing scenarios
rockfish_cv_range = c((mean(vary_fishing_rock_stock_cv_pcorr[,1]) - mean(vary_fishing_rock_stock_cv_pcorr[,length(fishing_trials)])), # Positive 
                   (mean(vary_fishing_rock_stock_cv_ncorr[,1]) - mean(vary_fishing_rock_stock_cv_ncorr[,length(fishing_trials)]))) # Negative

lingcod_cv_range = c((mean(vary_fishing_ling_stock_cv_pcorr[,1]) - mean(vary_fishing_ling_stock_cv_pcorr[,length(fishing_trials)])), # rockfish 
                     (mean(vary_fishing_ling_stock_cv_ncorr[,1]) - mean(vary_fishing_ling_stock_cv_ncorr[,length(fishing_trials)]))) # lingcod
names(rockfish_cv_range) = names(lingcod_cv_range) =c("positive", "negative")
cv_range = as.data.frame(abs(cbind(rockfish_cv_range,lingcod_cv_range))) %>% 
  set_colnames(c("rockfish", "lingcod")) %>% 
  mutate(correlation = c("positive", "negative")) %>% 
  pivot_longer(!correlation, names_to = "species", values_to = "mean.CV")
  
jpeg('plots/cv.range.varyfishpress.jpg')
ggplot(cv_range, aes(reorder(species, -mean.CV), mean.CV, fill = reorder(correlation, -mean.CV))) +
  geom_col(position = "dodge") +
  labs(y = "Difference in CV", x = "", title = "Change in CV", subtitle= "across fishing pressure scenarios") +
  theme_classic() +
  scale_fill_discrete(name = "Correlation")
dev.off()  
  






  
  
  