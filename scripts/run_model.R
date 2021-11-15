library(tidyverse)
library(miceadds)
library(abind)
library(magrittr)
library(viridis)

# Run deterministic model with varying size-selective predation --------------------------------------------

load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
load.Rdata( filename="cleaned_data/binned.size.spec.27.Rdata", "binned.size.spec.27" ) 
load.Rdata( filename="cleaned_data/binned.size.spec.31.Rdata", "binned.size.spec.31" ) 
load.Rdata( filename="cleaned_data/binned.size.spec.33.Rdata", "binned.size.spec.33" ) 
diet.frac.rockfish <- c(rep(0, 4), rep(0.001, (lingcod$nage-4)), rep(0, 4), rep(0.001, (lingcod$nage-4)))
a_ij.27 = binned.size.spec.27 %*% diag(diet.frac.rockfish) 
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish)
a_ij.31 = binned.size.spec.31 %*% diag(diet.frac.rockfish)
a_ij.33 = binned.size.spec.33 %*% diag(diet.frac.rockfish)

a_ij_trials = abind(a_ij.27, a_ij.33, along = 3) 
vary_dietspec_rockfish = matrix(NA, nrow = 300, ncol = 2)

for(t in 1:2) {
get_pop_n(rockfish, lingcod, nsim = 1, init.l = 10, init.r = 20, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
          tf = 300, fish.mort = 0.1, b = 0.1, a_ij = a_ij_trials[,,t], handling = 0.01, times = 1:2, 
          stochastic = FALSE)
  vary_dietspec_rockfish[,t] = rockfish_pop
}

# Pivot dataset into long format
vary_dietspec_rockfish = as.data.frame(vary_dietspec_rockfish) %>% 
  set_colnames(c("lower", "higher")) %>% 
  mutate(time = 1:300) %>% 
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




# Run deterministic model with varying fishing pressure ------------------------------------
fishing_trials = seq(0, 0.5, by = 0.1)
vary_fishing_rock_det = matrix(NA, nrow = 300, ncol = length(fishing_trials))
vary_fishing_ling_det = matrix(NA, nrow = 300, ncol = length(fishing_trials))

for (f in 1:length(fishing_trials)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = 10, init.r = 20, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = 300, fish.mort = fishing_trials[f], b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = FALSE)
  vary_fishing_rock_det[,f] = rockfish_pop
  vary_fishing_ling_det[,f] = lingcod_pop
}

vary_fishing_rock_det_long = as.data.frame(vary_fishing_rock_det) %>% 
  set_colnames(c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:300) %>% 
  pivot_longer(!time, names_to = "Fishing.Pressure", values_to = "Abundance")

vary_fishing_ling_det_long = as.data.frame(vary_fishing_ling_det) %>% 
  set_colnames(c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) %>% 
  mutate(time = 1:300) %>% 
  pivot_longer(!time, names_to = "Fishing.Pressure", values_to = "Abundance")

# Plot
jpeg('plots/det.vary.fishing.rockf.jpg')
ggplot(vary_fishing_rock_det_long, aes(x=time, y = Abundance, color = Fishing.Pressure)) + 
  geom_line() + 
  labs(y = "Abundance", x = "Years", title = "Rockfish Time-Series", 
       subtitle = "Varying Fishing Pressure") +
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
vary_handling_rock_det = matrix(NA, nrow = 300, ncol = length(h_trials))

for (h in 1:length(h_trials)) {
  get_pop_n(rockfish, lingcod, nsim = 1, init.l = 10, init.r = 20, corr = 0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = 300, fish.mort = 0.2, b = 0.1, a_ij = a_ij.29, handling = h_trials[h], times = 1:2, 
            stochastic = FALSE)
  vary_handling_rock_det[,h] = rockfish_pop
}

vary_handling_rock_det_long = as.data.frame(vary_handling_rock_det) %>% 
  set_colnames(c("0.01", "0.05", "0.15", "0.25")) %>% 
  mutate(time = 1:300) %>% 
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



# Run stochastic model with varying fishing pressure -----------------------------------------
fishing_trials = seq(0, 0.5, by = 0.1)
vary_fishing_rock_stoch = matrix(NA, nrow = 300, ncol = length(fishing_trials))
vary_fishing_ling_stoch = matrix(NA, nrow = 300, ncol = length(fishing_trials))
vary_fishing_rock_stock_cv = matrix(NA, nrow = 10, ncol = length(fishing_trials)) # nrow is number of simulations
vary_fishing_ling_stock_cv = matrix(NA, nrow = 10, ncol = length(fishing_trials)) # nrow is number of simulations

for(f in 1:length(fishing_trials)) {
 # Run the function over each fishing pressure scenario 
  get_pop_n(rockfish, lingcod, nsim = 10, init.l = 10, init.r = 20, corr = -0.8, autocorr = c(0.23,0.23), cv = 0.6, mn = 0.5*cv^2, 
            tf = 300, fish.mort = fishing_trials[f], b = 0.1, a_ij = a_ij.29, handling = 0.01, times = 1:2, 
            stochastic = TRUE)
  #vary_fishing_rock_stoch[,f] = rockfish_pop
  #vary_fishing_ling_stoch[,f] = lingcod_pop
 # Toss out the initial 25 years as a burn-in phase
  rockfish_pop_subset = rockfish_pop[, -(1:25)]
  lingcod_pop_subset = lingcod_pop[, -(1:25)]
 # Matrix as output - but in the forloop, all we want is the coefficient of variation 
    # (and maybe mean)
  # Which will look like matrix with fishing_trials columns and nsims rows for CV 
  # and similarly for mean
  vary_fishing_rock_stock_cv[, f] = apply(rockfish_pop_subset, 1, function(x) sd(x) / mean(x) * 100)
  vary_fishing_ling_stock_cv[, f] = apply(lingcod_pop_subset, 1, function(x) sd(x) / mean(x) * 100)
}

# saving the data - adjust the names to fit the correct run (i.e. whether is positive and negative correlation)
# vary_fishing_rock_stock_cv_ncorr = vary_fishing_rock_stock_cv
# vary_fishing_ling_stock_cv_ncorr = vary_fishing_ling_stock_cv
# save(vary_fishing_rock_stock_cv_ncorr, file = "cleaned_data/vary_fishing_rock_stock_cv_ncorr.Rdata")  
# save(vary_fishing_ling_stock_cv_ncorr, file = "cleaned_data/vary_fishing_ling_stock_cv_ncorr.Rdata")   

# rockfish
jpeg('plots/cv.stoch.varyfishpress.rock.jpg')  
plot(fishing_trials,colMeans(vary_fishing_rock_stock_cv_pcorr), type = "l", ylab = "mean CV", xlab = "Fishing Pressure", main = "Rockfish")  
lines(fishing_trials,colMeans(vary_fishing_rock_stock_cv_ncorr), col = "blue")
legend(0.3, 75, legend=c("positive", "negative"),
       col=c("black", "blue"), lty=1)
dev.off()

# lingcod
jpeg('plots/cv.stoch.varyfishpress.ling.jpg') 
plot(fishing_trials,colMeans(vary_fishing_ling_stock_cv_pcorr), type = "l", ylab = "mean CV", xlab = "Fishing Pressure", main = "Lingcod")  
lines(fishing_trials,colMeans(vary_fishing_ling_stock_cv_ncorr), col = "blue")
legend(0.3, 27, legend=c("positive", "negative"),
       col=c("black", "blue"), lty=1)
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
  
  
  
  