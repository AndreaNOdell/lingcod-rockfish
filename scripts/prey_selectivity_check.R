# prey selectivity effect no fishing

load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/rockfish_parms.Rdata")
source("scripts/model.R")

SBt = 0.4 # Specified SBt/SB0eq value of interest. In most cases, 0.4
nsims = 150 # number of simulations
cv = 0.5
autocorr_R = 0.23
mpa.year = 0
rec.year = 25
handl = 0.3
bycatch = 0.05

# Run model ----------
no_prey = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                     cv = cv, tf = (5+mpa.year+350), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = bycatch, 
                     quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                     rockfish.prop = 0.000001, yelloweye.prop = 0.02)
no_prey$prop_crashed = sum(no_prey$rockfish_biomass_ts[,ncol(no_prey$rockfish_biomass_ts)] == 0)/nsims
#save(no_prey, file = "results/prey_selectivity/no_prey.RData")


low_prey = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                              cv = cv, tf = (5+mpa.year+350), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = bycatch, 
                              quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                              rockfish.prop = 0.05, yelloweye.prop = 0.02)
low_prey$prop_crashed = sum(low_prey$rockfish_biomass_ts[,ncol(low_prey$rockfish_biomass_ts)] == 0)/nsims
#save(low_prey, file = "results/prey_selectivity/low_prey.RData")

high_prey = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                       cv = cv, tf = (5+mpa.year+350), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = bycatch, 
                       quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                       rockfish.prop = 0.2, yelloweye.prop = 0.07)
high_prey$prop_crashed = sum(high_prey$rockfish_biomass_ts[,ncol(high_prey$rockfish_biomass_ts)] == 0)/nsims
#save(high_prey, file = "results/prey_selectivity/high_prey.RData")

xtra_high_prey = get_pop_ts(rockfish, lingcod, nsim = nsims, corr = 0, autocorr_L = 0.23, autocorr_R = autocorr_R, 
                             cv = cv, tf = (5+mpa.year+350), rec.yr = rec.year, mpa.yr = mpa.year, hist.f = 0.5, hist.by = 0.5, f = 0, b = bycatch, 
                             quant95 = 0.29, handling = handl, min.selectivity = TRUE, min.age.consumed = 4, 
                             rockfish.prop = 0.1, yelloweye.prop = 1)
#save(xtra_high_prey, file = "results/prey_selectivity/xtra_high_prey.RData")


# Plot ----------
load(file = "results/prey_selectivity/no_prey.RData")
load(file = "results/prey_selectivity/low_prey.RData")
load(file = "results/prey_selectivity/high_prey.RData")
load(file = "results/prey_selectivity/xtra_high_prey.RData")

prey_selectivity_df = as.data.frame(cbind("spawning biomass" = c(no_prey$rockfish_avg/no_prey$rockfish_avg, low_prey$rockfish_avg/no_prey$rockfish_avg, high_prey$rockfish_avg/no_prey$rockfish_avg, xtra_high_prey$rockfish_avg/no_prey$rockfish_avg),
                        "variability" = c(no_prey$rockfish_cv/no_prey$rockfish_cv, low_prey$rockfish_cv/no_prey$rockfish_cv, high_prey$rockfish_cv/no_prey$rockfish_cv, xtra_high_prey$rockfish_cv/no_prey$rockfish_cv),
                        "age-structure" = c(no_prey$rockfish_age/no_prey$rockfish_age, low_prey$rockfish_age/no_prey$rockfish_age, high_prey$rockfish_age/no_prey$rockfish_age, xtra_high_prey$rockfish_age/no_prey$rockfish_age))) %>% 
  mutate(selectivity = c("no", "low", "high", "extra high")) %>% 
  pivot_longer(cols = c("spawning biomass", "variability", "age-structure"), names_to = "outcome", values_to = "value") %>%
  mutate(selectivity = fct_relevel(selectivity, c("no", "low", "high", "extra high"))) %>% 
  mutate(outcome = fct_relevel(outcome, c("spawning biomass", "variability", "age-structure")))


jpeg("plots/prey_selectivity/prey_selectivity_on_dynamics.jpeg", units="in", width=6, height=4, res = 300)
prey_selectivity_df %>% 
ggplot(aes(x = selectivity, y = value, col = selectivity)) +
  theme_classic() +
  geom_point( size = 2.25, position = position_dodge(0.7)) +
  facet_wrap( ~ outcome, strip.position = "bottom", scales = "free_x") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "gray90",
                                        colour = "gray90",
                                        size = 0.5, linewidth = "solid")) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey44") +
  labs(title = "Effect of Prey Selectivity on Dynamics") +
  scale_color_manual(values=c("dodgerblue", "blue", "blue4", "navyblue"), breaks=c('no', 'low', 'high', 'extra high'), name = "Prey Selectivity")
dev.off()

# look to see if high prey selectivity effect on variability is due 
# to high cv of lingcod






# FIXED - Collapsed pops ----------
example_ts = as.data.frame(rbind(high_prey$rockfish_biomass_ts[high_prey$rockfish_biomass_ts[,350] == 0,], 
      high_prey$rockfish_biomass_ts[1:10,]))
write_csv(example_ts, file ="results/prey_selectivity/example_ts.csv")

matplot(t(as.matrix(example_ts)), type = "l")

collapsed_pops = as.data.frame(high_prey$rockfish_biomass_ts) %>% 
  mutate(row = 1:nsims) %>% 
  filter(V350 == 0)
collapsed_pops$row

example_ts = as.data.frame(rbind(high_prey$rockfish_biomass_ts[collapsed_pops$row[1:5],], 
                    high_prey$rockfish_biomass_ts[1:5,])) %>% 
  mutate(persist = c(rep("collapsed", length(collapsed_pops$row[1:5])), rep("persist", 5))) %>%
  mutate(sim = 1:(length(collapsed_pops$row[1:5])+5)) %>% 
  pivot_longer(!c(persist,sim), names_to = "time", values_to = "SB") %>% 
  mutate(species = "rockfish")
example_ts$time = gsub("V","",as.character(example_ts$time))
example_ts$time = as.numeric(example_ts$time)
example_ts$sim = as.factor(example_ts$sim)

example_ts_lingcod = as.data.frame(rbind(high_prey$lingcod_biomass_ts[collapsed_pops$row[1:5],], 
                    high_prey$lingcod_biomass_ts[1:5,])) %>% 
  mutate(persist = c(rep("collapsed", length(collapsed_pops$row[1:5])), rep("persist", 5))) %>%
  mutate(sim = 1:(length(collapsed_pops$row[1:5])+5)) %>% 
  pivot_longer(!c(persist,sim), names_to = "time", values_to = "SB") %>% 
  mutate(species = "lingcod")
example_ts_lingcod$time = gsub("V","",as.character(example_ts_lingcod$time))
example_ts_lingcod$time = as.numeric(example_ts_lingcod$time)
example_ts_lingcod$sim = as.factor(example_ts_lingcod$sim)

example_ts_df = bind_rows(example_ts, example_ts_lingcod)

library(viridis)

jpeg("plots/prey_selectivity/collapse_rockfish_lingcod.jpeg", units="in", width=7, height=4, res = 300)
ggplot(example_ts_df, aes(x = time, y = SB)) +
  geom_line(aes(linetype = persist, color = sim)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  facet_wrap( ~ species, strip.position = "bottom", scales = "free_y") +
  scale_color_viridis(discrete = TRUE)
dev.off()

high_prey$lingcod_biomass_ts[collapsed_pops$row,]
