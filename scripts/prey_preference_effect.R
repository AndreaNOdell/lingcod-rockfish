load(file = "results/fishing_scenarios/trial_no_fishing.RData")
load(file = "results/fishing_scenarios/trial_light_fishing.RData")
load(file = "results/fishing_scenarios/trial_moderate_fishing.RData")

light_fishing = trial_light_fishing
moderate_fishing = trial_moderate_fishing
no_fishing = trial_no_fishing

load(file = "results/fishing_scenarios/trial_no_fishing_highpref.RData")
load(file = "results/fishing_scenarios/trial_light_fishing_highpref.RData")
load(file = "results/fishing_scenarios/trial_moderate_fishing_highpref.RData")
load(file = "results/fishing_scenarios/trial_moderate_fishing_nobycatch_highpref.RData")

rm(list=ls(pattern="trial"))

# Biological Outcomes ---------

cv = cbind(c(no_fishing$rockfish_cv/no_fishing$rockfish_cv, light_fishing$rockfish_cv/no_fishing$rockfish_cv, moderate_fishing$rockfish_cv/no_fishing$rockfish_cv,
             no_fishing_highpref$rockfish_cv/no_fishing_highpref$rockfish_cv, light_fishing_highpref$rockfish_cv/no_fishing_highpref$rockfish_cv, moderate_fishing_highpref$rockfish_cv/no_fishing_highpref$rockfish_cv),
           rep(c("no", "light", "moderate"),2),
           rep("variability", 6),
           c(rep("base",3), rep("high",3)))


age = cbind(c(no_fishing$rockfish_age/no_fishing$rockfish_age, light_fishing$rockfish_age/no_fishing$rockfish_age, moderate_fishing$rockfish_age/no_fishing$rockfish_age,
              no_fishing_highpref$rockfish_age/no_fishing_highpref$rockfish_age, light_fishing_highpref$rockfish_age/no_fishing_highpref$rockfish_age, moderate_fishing_highpref$rockfish_age/no_fishing_highpref$rockfish_age),
            rep(c("no", "light", "moderate"),2),
            rep("age-structure", 6),
            c(rep("base",3), rep("high",3)))


sb = cbind(c(no_fishing$rockfish_avg/no_fishing$rockfish_avg, light_fishing$rockfish_avg/no_fishing$rockfish_avg, moderate_fishing$rockfish_avg/no_fishing$rockfish_avg,
             no_fishing_highpref$rockfish_avg/no_fishing_highpref$rockfish_avg, light_fishing_highpref$rockfish_avg/no_fishing_highpref$rockfish_avg, moderate_fishing_highpref$rockfish_avg/no_fishing_highpref$rockfish_avg),
           rep(c("no", "light", "moderate"),2),
           rep("spawning biomass", 6),
           c(rep("base",3), rep("high",3)))


df = as.data.frame(rbind(cv, age, sb))
colnames(df) = c("value", "fishing", "variable", "prey_selectivity")
df$value = as.numeric(df$value)

jpeg("plots/fishing_scenarios/prey_preference_effect_onlybase.jpeg", units="in", width=6, height=4, res = 300)
df %>% 
  mutate(fishing = fct_relevel(fishing, 
                            "no", "light", "moderate")) %>% 
  mutate(variable = fct_relevel(variable,
                                "spawning biomass", "variability", "age-structure")) %>%
  mutate(prey_selectivity = fct_relevel(prey_selectivity,
                                "base", "high")) %>% 
ggplot(aes(x = fishing, y = value, color = fishing, alpha = prey_selectivity)) +
  theme_classic() +
  geom_point( size = 2.25, position = position_dodge(0.7)) +
  facet_wrap( ~ variable, strip.position = "bottom", scales = "free_x") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "gray90",
                                        colour = "gray90",
                                        size = 0.5, linewidth = "solid")) +
  scale_alpha_discrete(range = c(1, 0), name = "prey selectivity") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey44") +
  scale_color_manual(values=c("dodgerblue", "blue", "blue4"), breaks=c('no', 'light', 'moderate'), name = "fishing pressure")
dev.off()

jpeg("plots/fishing_scenarios/prey_preference_effect.jpeg", units="in", width=6, height=4, res = 300)
df %>% 
  mutate(fishing = fct_relevel(fishing, 
                               "no", "light", "moderate")) %>% 
  mutate(variable = fct_relevel(variable,
                                "spawning biomass", "variability", "age-structure")) %>%
  mutate(prey_selectivity = fct_relevel(prey_selectivity,
                                        "base", "high")) %>% 
  ggplot(aes(x = fishing, y = value, color = fishing, shape = prey_selectivity)) +
  theme_classic() +
  geom_point( size = 2.25, position = position_dodge(0.7)) +
  facet_wrap( ~ variable, strip.position = "bottom", scales = "free_x") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "gray90",
                                        colour = "gray90",
                                        size = 0.5, linewidth = "solid")) +
  scale_shape_manual(values = c(19,10), name = "prey selectivity") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey44") +
  scale_color_manual(values=c("dodgerblue", "blue", "blue4"), breaks=c('no', 'light', 'moderate'), name = "fishing pressure")
dev.off()



