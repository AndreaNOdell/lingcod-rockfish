library(tidyverse)
library(miceadds)
library(abind)
library(magrittr)

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

vary_dietspec_rockfish = as.data.frame(vary_dietspec_rockfish) %>% 
  set_colnames(c("lower", "higher")) %>% 
  mutate(time = 1:300) %>% 
  pivot_longer(!time, names_to = "size.spec", values_to = "abundance")

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



# 







# 

