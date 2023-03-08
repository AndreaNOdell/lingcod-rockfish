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

# The effect of handling time on spawning biomass
ggplot(GSA_total_df, aes(x = handling, y = rockfish_avg)) +
  geom_point() +
  theme_classic()
# for the handling time, are these values reasonable? The effect on spawning biomass
# seems to occur in such a narrow range of handling values which seems
# to be driving their importance indicated from the random forest.

ggplot(GSA_total_df, aes(x = rockfish.prop, y = rockfish_avg, color = handling)) +
  geom_point() +
  theme_classic()
# the decrease in long term spawning biomass is driven in part by handling and yelloweye rockfish proportion

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

