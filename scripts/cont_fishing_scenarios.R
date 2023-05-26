library(tidyverse)
library(RColorBrewer)
library(latex2exp)

# Set up ---------

load("results/fishing_scenarios_cont/fishing_base_selec.Rdata")
load("results/fishing_scenarios_cont/fishing_extreme_selec.Rdata")
load("results/fishing_scenarios_cont/fishing.Rdata")

# extract data 
base_avg = unname(unlist(lapply(fishing_base_selec, `[`, "rockfish_avg")))
base_cv = unname(unlist(lapply(fishing_base_selec, `[`, "rockfish_cv")))
base_age = unname(unlist(lapply(fishing_base_selec, `[`, "rockfish_age")))

extreme_avg = unname(unlist(lapply(fishing_extreme_selec, `[`, "rockfish_avg")))
extreme_cv = unname(unlist(lapply(fishing_extreme_selec, `[`, "rockfish_cv")))
extreme_age = unname(unlist(lapply(fishing_extreme_selec, `[`, "rockfish_age")))

# F40% calculated below
f40 = round(round(fishing,2) / 0.27, 2)

# make datasets
df_base_avg = cbind(base_avg/base_avg[1], fishing, rep("spawning biomass", length(fishing)), rep("base"), f40)
df_base_cv = cbind(base_cv/base_cv[1], fishing, rep("variability", length(fishing)), rep("base"), f40)
df_base_age = cbind(base_age/base_age[1], fishing, rep("age-structure", length(fishing)), rep("base"), f40)

df_extreme_avg = cbind(extreme_avg/extreme_avg[1], fishing, rep("spawning biomass", length(fishing)), rep("extreme"), f40)
df_extreme_cv = cbind(extreme_cv/extreme_cv[1], fishing, rep("variability", length(fishing)), rep("extreme"), f40)
df_extreme_age = cbind(extreme_age/extreme_age[1], fishing, rep("age-structure", length(fishing)), rep("extreme"), f40)

# combine data into one master dataset
df = as.data.frame(rbind(df_base_avg, df_base_cv, df_base_age, df_extreme_avg, df_extreme_cv, df_extreme_age))
colnames(df) <- c("value", "fishing", "variable", "prey_selectivity", "f40")
df$value = as.numeric(df$value)
df$fishing = round(as.numeric(df$fishing),3)
df$f40 = as.numeric(df$f40)


fishing_vals = c(0, 0.48, 1)

df_point = df %>% 
  filter(f40 %in% fishing_vals) %>% 
  mutate("variable" = fct_relevel(variable,
                                  "spawning biomass", "variability", "age-structure"))


# Prey preference on fishing impact -------------------------------
  
# Make plot
jpeg("plots/fishing_scenarios/prey_preference_effect_contF.jpeg", units="in", width=6, height=4, res = 300)
df %>% 
  mutate("variable" = fct_relevel(variable,
                                "spawning biomass", "variability", "age-structure")) %>%
  mutate(prey_selectivity = fct_relevel(prey_selectivity,
                                        "base", "extreme")) %>% 
  ggplot(aes(x = f40, y = value, group = prey_selectivity)) +
  theme_classic() +
  geom_line(aes(linetype = prey_selectivity)) +
  facet_wrap( ~ variable, strip.position = "bottom", scales = "free_x") +
  geom_point(data = df_point, mapping = aes(x = f40, y = value, 
                                            col = fishing, 
                                            shape = prey_selectivity), 
             size = 3, stroke = 1) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "gray90",
                                        colour = "gray90",
                                        size = 0.5, linewidth = "solid")) +
  scale_shape_manual(values = c(19,1), name = "prey selectivity") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey44") +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  labs(x  = "F/F40%", y = "Relative Change") +
  guides(shape = FALSE, color = FALSE)
dev.off()

jpeg("plots/fishing_scenarios/prey_preference_effect_contF_baseonly.jpeg", units="in", width=6, height=4, res = 300)
df %>% 
  mutate("variable" = fct_relevel(variable,
                                  "spawning biomass", "variability", "age-structure")) %>%
  filter(prey_selectivity == "base") %>% 
  ggplot(aes(x = f40, y = value, group = prey_selectivity)) +
  theme_classic() +
  geom_line(aes(linetype = prey_selectivity)) +
  facet_wrap( ~ variable, strip.position = "bottom", scales = "free_x") +
  geom_point(data = df_point[df_point$prey_selectivity == "base",], mapping = aes(x = f40, y = value, 
                                            col = fishing, 
                                            shape = prey_selectivity), 
             size = 3, stroke = 1) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.background = element_rect(fill = "gray90",
                                        colour = "gray90",
                                        size = 0.5, linewidth = "solid")) +
  scale_shape_manual(values = c(19,1), name = "prey selectivity") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey44") +
  scale_colour_gradient(low = "#56B1F7", high = "#132B43",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour") +
  labs(x  = "F/F40%", y = "Relative Change")  +
  guides(shape = FALSE, color = FALSE)
dev.off()


# Rebuilding time --------
SBt = 0.4
f_rows = which(f40 %in% fishing_vals)

# base selectivity
rockfish_base_avg_df = lapply(fishing_base_selec, `[`, "rockfish_biomass_ts")
no_fishing_base_df <- rockfish_base_avg_df[[f_rows[1]]]$rockfish_biomass_ts/apply(rockfish_base_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)
light_fishing_base_df <- rockfish_base_avg_df[[f_rows[2]]]$rockfish_biomass_ts/apply(rockfish_base_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)
mod_fishing_base_df<- rockfish_base_avg_df[[f_rows[3]]]$rockfish_biomass_ts/apply(rockfish_base_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)

nofishing_SB40_base = as.data.frame(which(no_fishing_base_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40_base$t_SB40); sd(nofishing_SB40_base$t_SB40)

lightfishing_SB40_base = as.data.frame(which(light_fishing_base_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40_base$t_SB40); sd(lightfishing_SB40_base$t_SB40)

modfishing_SB40_base = as.data.frame(which(mod_fishing_base_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "moderate")
mean(modfishing_SB40_base$t_SB40); sd(modfishing_SB40_base$t_SB40)

# Extreme Selectivity
rockfish_extreme_avg_df = lapply(fishing_extreme_selec, `[`, "rockfish_biomass_ts")
no_fishing_extreme_df <- rockfish_extreme_avg_df[[f_rows[1]]]$rockfish_biomass_ts/apply(rockfish_extreme_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)
light_fishing_extreme_df <- rockfish_extreme_avg_df[[f_rows[2]]]$rockfish_biomass_ts/apply(rockfish_extreme_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)
mod_fishing_extreme_df<- rockfish_extreme_avg_df[[f_rows[3]]]$rockfish_biomass_ts/apply(rockfish_extreme_avg_df[[f_rows[1]]]$rockfish_biomass_ts[,-(1:150)], 1, mean)

nofishing_SB40_extreme = as.data.frame(which(no_fishing_extreme_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "none")
mean(nofishing_SB40_extreme$t_SB40); sd(nofishing_SB40_extreme$t_SB40)

lightfishing_SB40_extreme = as.data.frame(which(light_fishing_extreme_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "light")
mean(lightfishing_SB40_extreme$t_SB40); sd(lightfishing_SB40_extreme$t_SB40)

modfishing_SB40_extreme = as.data.frame(which(mod_fishing_extreme_df >= SBt,arr.ind = T)) %>% 
  group_by(row) %>% 
  summarise(t_SB40 = min(col)) %>% 
  mutate(fishing = "moderate")
mean(modfishing_SB40_extreme$t_SB40); sd(modfishing_SB40_extreme$t_SB40)


# combine datasets
SB40_df_base = as.data.frame(rbind(nofishing_SB40_base, lightfishing_SB40_base, modfishing_SB40_base)) %>% 
  mutate(selectivity = "base")
SB40_df_extreme = as.data.frame(rbind(nofishing_SB40_extreme, lightfishing_SB40_extreme, modfishing_SB40_extreme)) %>% 
  mutate(selectivity = "extreme")
SB40_df = as.data.frame(rbind(SB40_df_base,SB40_df_extreme))

jpeg("plots/fishing_scenarios/rebuildingtime_SB40_violin_selectivity_compare.jpeg", units="in", width=5, height=4, res = 300)
SB40_df %>% 
  mutate(selectivity = fct_relevel(selectivity, 
                                   "base", "extreme")) %>% 
  ggplot(aes(x=reorder(fishing, t_SB40), y=t_SB40, fill = fishing, color = fishing, alpha=selectivity)) + 
  geom_violin(position = "dodge", adjust = 1.1) +
  ylim(0,100) +
  theme_classic() +
  scale_fill_manual(values=c("#56B1F7", "#01579B", "#132B43"), breaks=c('none', 'light', 'moderate')) +
  scale_color_manual(values=c("#56B1F7", "#01579B", "#132B43"), breaks=c('none', 'light', 'moderate')) +
  scale_alpha_manual(values=c(1, 0)) +
  labs(x = "Lingcod Fishing Intensity", y = TeX(r'(Years till Rockfish 40% $SB_0$)')) + 
  guides(color = FALSE, fill = FALSE, alpha = FALSE) +
  scale_x_discrete(labels = c("none", TeX(r'(50% $F_{40}$)'), TeX(r'($F_{40}$)')))
dev.off()


# Check F/Fmsy for lingcod ---------
# In 2017 lingcod stock assessment for South population, SB0 = 20,260 mt and
# SBMSY = 5,265mt. So, BMSY = SB26%

lingcod_base_avg_df = (lapply(fishing_base_selec, `[`, "lingcod_biomass_ts"))
lincod_base_avg = numeric(length(tmp))
for(i in 1:length(lincod_base_avg)) {
  lincod_base_avg[i] = mean(apply(lingcod_base_avg_df[[i]]$lingcod_biomass_ts[,-(1:150)], 1, mean))
}

cbind(SB = round(lincod_base_avg/lincod_base_avg[1],2), fishing = round(fishing,2))

f40 = round(round(fishing,2) / 0.27, 2)

# it looks like SB40% is achieved by F40% = 0.27






