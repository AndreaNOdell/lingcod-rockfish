library(tidyverse)
library(randomForest)
library(rpart)
library(rpart.plot)
library(caret)
library(plotrix)
library(ggpubr)
require(grid) 

load("GSA_results_final/GSA_total_df.Rdata")
load("GSA_results_final/RF_SBeq_total_data.Rdata")
load("GSA_results_final/RF_CVeq_total_data.Rdata")
load("GSA_results_final/RF_age_total_data.Rdata")

# Data exploration -------------------------------------

#### histograms of outcomes
cv_cv = (sd(RF_CVeq_total_data[,1])/mean(RF_CVeq_total_data[,1]))*100
cv_hist = ggplot(GSA_total_df, aes(rockfish_cv)) +
  geom_histogram(bins = 20, fill = "gray24", col = "white") +
  theme_classic() +
  labs(title = "Stability", subtitle = paste("cv =", round(cv_cv, 3)), x = "CV of Spawning Biomass")


cv_sb = (sd(RF_SBeq_total_data[,1])/mean(RF_SBeq_total_data[,1]))*100
avg_hist = ggplot(GSA_total_df, aes(rockfish_avg)) +
  geom_histogram(bins = 20, fill = "gray24", col = "white") +
  theme_classic() +
  labs(title = "Spawning Biomass", subtitle = paste("cv =", round(cv_sb, 3)), x = "Avg of Spawning Biomass")


cv_age = (sd(RF_age_total_data[,1])/mean(RF_age_total_data[,1]))*100
age_hist = ggplot(GSA_total_df, aes(rockfish_age)) +
  geom_histogram(bins = 20, fill = "gray24", col = "white") +
  theme_classic() +
  labs(title = "Age-Structure", subtitle = paste("cv =", round(cv_age, 3)), x = "Avg of Spawning Biomass")

jpeg(filename = "plots/GSA/outcome_histograms.jpeg", units="in", width=8, height=4, res = 300)
figure = ggarrange(avg_hist + rremove("ylab"), cv_hist + rremove("ylab"), age_hist + rremove("ylab"), ncol = 3)
annotate_figure(figure, left = textGrob("count", rot = 90, vjust = 1, gp = gpar(cex = 1)))
dev.off()

# Random forest -----------------------------------------

RF_avg <- randomForest(rockfish_avg~.,data=RF_SBeq_total_data, mtry = 4, importance = TRUE, proximity=TRUE)
RF_cv <- randomForest(rockfish_cv~.,data=RF_CVeq_total_data, mtry = 4, importance=TRUE,proximity=TRUE)
RF_age <- randomForest(rockfish_age~.,data=RF_age_total_data, mtry = 4, importance=TRUE,proximity=TRUE)

# Importance Plots ---------------------------------------

### from Jorge - barplots with unnormalized 

SBeq_GSA <-
  data.frame(Parameter = rownames(RF_avg$importance),
             Importance = RF_avg$importance[, 1])
SBeq_GSA$Parameter <-
  factor(SBeq_GSA$Parameter, levels = SBeq_GSA$Parameter[order(SBeq_GSA$Importance,decreasing = TRUE)])
jpeg(filename = "plots/GSA/importance_plots_SBeq.jpeg", units="in", width=4, height=3, res = 300)
SBeq_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  ) +
  labs(subtitle  = "Long-term Spawning Biomass")
dev.off()

cv_GSA <-
  data.frame(Parameter = rownames(RF_cv$importance),
             Importance = RF_cv$importance[, 1])
cv_GSA$Parameter <-
  factor(cv_GSA$Parameter, levels = cv_GSA$Parameter[order(cv_GSA$Importance,decreasing = TRUE)])
jpeg(filename = "plots/GSA/importance_plots_CVeq.jpeg", units="in", width=4, height=3, res = 300)
cv_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  ) +
  labs(subtitle  = "Long-term Stability")
dev.off()


age_GSA <-
  data.frame(Parameter = rownames(RF_age$importance),
             Importance = RF_age$importance[, 1])
age_GSA$Parameter <-
  factor(age_GSA$Parameter, levels = age_GSA$Parameter[order(age_GSA$Importance,decreasing = TRUE)])
jpeg(filename = "plots/GSA/importance_plots_age.jpeg", units="in", width=4, height=3, res = 300)
age_GSA %>% ggplot(aes(Parameter, Importance)) + geom_col() + theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=16),
    axis.text.x = element_text(angle = 90)
  ) +
  labs(subtitle  = "Long-term Age-structure")
dev.off()

### Lineplots of normalized %IncMSE

ImpData <- as.data.frame(importance(RF_age))
ImpData$Var.Names <- row.names(ImpData)

age_lineplot = ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  ylim(-4, 400) +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=10)
  ) +
  labs(title = "Age-Structure")

ImpData <- as.data.frame(importance(RF_avg))
ImpData$Var.Names <- row.names(ImpData)

avg_lineplot = ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  ylim(-4, 400) +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=10)
  ) +
  labs(title = "Spawning Biomass")


ImpData <- as.data.frame(importance(RF_cv))
ImpData$Var.Names <- row.names(ImpData)

cv_lineplot = ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  ylim(-4, 400) +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size=10)
  ) +
  labs(title = "Stability")


jpeg(filename = "plots/GSA/importance_lineplots_aggregated.jpeg", units="in", width=4, height=8, res = 300)
figure = ggarrange(avg_lineplot + rremove("ylab") + rremove("xlab"), 
                   cv_lineplot + rremove("ylab") + rremove("xlab"),
                   age_lineplot + rremove("ylab") + rremove("xlab"), 
                   ncol = 1, 
                   common.legend = TRUE, legend = "right")
annotate_figure(figure, left = textGrob("parameter", rot = 90, vjust = 1, gp = gpar(cex = 0.9)),
                bottom = textGrob("%IncMSE", gp = gpar(cex = 0.9)))
dev.off()


# Important Parameter influence ------------------------------------

#### Spawning Biomass ####

# The effect of handling time on spawning biomass
jpeg(filename = "plots/GSA/SB_handling.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = handling, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass")
dev.off()
# for the handling time, are these values reasonable? The effect on spawning biomass
# seems to occur in such a narrow range of handling values which seems
# to be driving their importance indicated from the random forest.

# The effect of yelloweye and rockfish proportion on spawning biomass
jpeg(filename = "plots/GSA/SB_yelloweye.prop.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = rockfish.prop, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "proportion rockfish")
dev.off()
# The effect of yelloweye.prop (the proportion of rockfish that is yelloweye) and 
# rockfish proportion (the proportion of lingcod diet that is rockfish)

jpeg(filename = "plots/GSA/SB_cv.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = cv, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "recruitment variability")
dev.off()
# there is not a strong effect of recruitment variability and 
# autocorrelation on spawning biomass. The influence between 
# the two drivers are very similar

#The effect of both rockfish proportion in diet and handling time on spawning biomass
ggplot(GSA_total_df, aes(x = rockfish.prop, y = rockfish_avg, color = handling)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "proportion rockfish") +
  scale_color_continuous(name = "handling")

# the decrease in long term spawning biomass is driven in part by handling and yelloweye rockfish proportion


#### Stability #####

# the effect of recruitment variability and autocorrelation on stability
jpeg(filename = "plots/GSA/CV_cv.x.autocorrR.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = cv, y = rockfish_cv, col = autocorr_R)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "recruitment variability") +
  scale_color_continuous(name = "autocorr")
dev.off()
# Increases in recruitment variability and autocorrelation reduce stability

jpeg(filename = "plots/GSA/CV_autocorrR.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = autocorr_R, y = rockfish_cv)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "recruitment autocorrelation")
dev.off()


#### Age structure ####

# the effect of recruitment variability and autocorrelation on age-structure
jpeg(filename = "plots/GSA/AGE_cv.x.autocorrR.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = cv, y = rockfish_age, col = autocorr_R)) +
  geom_point() +
  theme_classic() +
  labs(y = "Proportion in plus group", x = "recruitment variability") +
  scale_color_continuous(name = "autocorr")
dev.off()
# Increases in recruitment variability and autocorrelation increase proportion 
# of individuals in plus group, but the effect is small

jpeg(filename = "plots/GSA/AGE_autocorrR.jpeg", units="in", width=5, height=4, res = 300)
ggplot(GSA_total_df, aes(x = autocorr_R, y = rockfish_age)) +
  geom_point() +
  theme_classic() +
  labs(y = "Proportion in plus group", x = "recruitment autocorrelation")
dev.off()

# Other Parameter influence -----------------------------------------------

#### Spawning Biomass ####
ggplot(GSA_total_df, aes(x = corr, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "recruitment correlation")

ggplot(GSA_total_df, aes(x = autocorr_L, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "lingcod autocorrelation")

ggplot(GSA_total_df, aes(x = quant95, y = rockfish_avg)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass", x = "size-selectivity")

#### Stability ####
ggplot(GSA_total_df, aes(x = corr, y = rockfish_cv)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "recruitment correlation")

ggplot(GSA_total_df, aes(x = autocorr_L, y = rockfish_cv)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "lingcod autocorrelation")

ggplot(GSA_total_df, aes(x = quant95, y = rockfish_cv)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "size-selectivity")

ggplot(GSA_total_df, aes(x = handling, y = rockfish_cv)) +
  geom_point() +
  theme_classic() +
  labs(y = "spawning biomass cv", x = "handling time")

#### Age-structure ####
