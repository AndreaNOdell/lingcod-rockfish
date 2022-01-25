library(miceadds)
source("scripts/scratch_script.R")
load.Rdata( filename="cleaned_data/binned.size.spec.Rdata", "binned.size.spec.29" ) # Load in rockfish size-spectra in size-specific lingod diet
diet.frac.rockfish <- c(rep(0, 4), rep(0.135*0.001, (lingcod$nage-4)), rep(0, 4), rep(0.165*0.001, (lingcod$nage-4)))
a_ij.29 = binned.size.spec.29 %*% diag(diet.frac.rockfish) # with interaction

# Unfished equilibrium ---------------------------------------------------------------------------------
unfished_vals = get_pop_det(rockfish, lingcod, init.l = 200, init.r = 70,  tf = 300, mpa.yr = 20, hist.f = 0, hist.by = 0.5, 
                            f = 0, b = 0, a_ij = a_ij.29, handling = 0.01, times = 1:2, min.selectivity = TRUE) # I checked and 300 yrs long enough
lingcod_unfished_equil = unfished_vals$lingcod_Beq
rockfish_unfished_equil = unfished_vals$rockfish_Beq
lingcod_unfished_oldprop = lingcod_oldprop
rockfish_unfished_oldprop = rockfish_oldprop
lingcod_unfished_SAD = lingcod_SAD
rockfish_unfished_SAD = rockfish_SAD

# Quantify balance between fishing and bycatch rate ---------------------------------------------------

f = seq(0.01,0.3, by = 0.03)
b = seq(0,0.3, by = 0.01)
harvest.scenarios = as.data.frame(expand_grid(f,b))

fxb_balance = mapply(get_fxb_balance, f = harvest.scenarios[,1], b = harvest.scenarios[,2], 
                     MoreArgs = list(rockfish = rockfish, lingcod = lingcod), SIMPLIFY = FALSE)

fxb_balance_HS = mapply(get_fxb_balance, f = harvest.scenarios[,1], b = harvest.scenarios[,2], 
                     MoreArgs = list(rockfish = rockfish, lingcod = lingcod, min.selectivity = FALSE), SIMPLIFY = FALSE)
                     

fxb_balance_df = cbind(harvest.scenarios, round(unlist(fxb_balance))/round(rockfish_unfished_equil))
names(fxb_balance_df) = c("f", "b", "rockfish_recovery")

fxbplot = fxb_balance_df %>% 
  filter(rockfish_recovery > 0.97) %>% 
  filter(rockfish_recovery < 1.03) %>% 
  ggplot(aes(x = f, y = b)) +
  geom_point() +
  theme_classic()











