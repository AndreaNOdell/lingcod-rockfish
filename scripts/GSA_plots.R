load("GSA_results/PredPrey_GSA.RData")
library(tidyverse)

# Importance plots of SB equilibrium
SBeq_GSA <-
  data.frame(Parameter = rownames(RF_avg$importance),
             Importance = RF_avg$importance[, 1])
SBeq_GSA$Parameter <-
  factor(SBeq_GSA$Parameter, levels = SBeq_GSA$Parameter[order(SBeq_GSA$Importance,decreasing = TRUE)])
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



# Importance plots of cv
cv_GSA <-
  data.frame(Parameter = rownames(RF_cv$importance),
             Importance = RF_cv$importance[, 1])
cv_GSA$Parameter <-
  factor(cv_GSA$Parameter, levels = cv_GSA$Parameter[order(cv_GSA$Importance,decreasing = TRUE)])
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
  





