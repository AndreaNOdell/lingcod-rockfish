library(tidyverse)
library(randomForest)
library(rpart)
library(rpart.plot)
library(caret)

# Create dataframe ------------------------------------------

RF_SBeq_total_data = matrix(NA, nrow = 2000, ncol = 11)
RF_CVeq_total_data = matrix(NA, nrow = 2000, ncol = 11)

#load in all data from slurm job array and concatenate results from 
# different slurm jobs
for(i in 0:19) {
load(paste0("GSA_results/SBeq_data/RF_SBeq_data_", i ,".Rdata"))
RF_SBeq_total_data[(1:100) + (i*100), ] = as.matrix(RF_SBeq_data)
}
colnames(RF_SBeq_total_data) = colnames(RF_SBeq_data)

for(i in 0:19) {
  load(paste0("GSA_results/CVeq_data/RF_cv_data_", i ,".Rdata"))
  RF_CVeq_total_data[(1:100) + (i*100), ] = as.matrix(RF_cv_data)
}
colnames(RF_CVeq_total_data) = colnames(RF_cv_data)



# Random forest --------------------------------------------

RF_avg <- randomForest(rockfish_avg~.,data=RF_SBeq_total_data, mtry = 4, importance = TRUE, proximity=TRUE)
#save(RF_avg, file = "GSA_results/random_forest/RF_avg.RData")
RF_cv <- randomForest(rockfish_cv~.,data=RF_CVeq_total_data, mtry = 4, importance=TRUE,proximity=TRUE)
#save(RF_cv, file = "GSA_results/random_forest/RF_cv.RData")

# get best mtry
mtry <- tuneRF(RF_CVeq_total_data[-1],RF_CVeq_total_data$rockfish_cv, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m)

#Evaluate variable importance
importance(RF_cv)
varImpPlot(RF_cv)

# CART ----------------------------------------------------

RF_SBeq_total_data = as.data.frame(RF_SBeq_total_data)
RF_CVeq_total_data = as.data.frame(RF_CVeq_total_data)

# Step1: Begin with a small cp. 
set.seed(123)
tree <- rpart(rockfish_avg ~ ., data = RF_SBeq_total_data, control = rpart.control(cp = 0.008))
printcp(tree)
bestcp <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
tree.pruned <- prune(tree, cp = bestcp)

plot(tree.pruned)
text(tree.pruned, cex = 0.3, use.n = TRUE, xpd = TRUE)
prp(tree.pruned, faclen = 0, cex = 0.3, extra = 1)


# pair wise plots --------------------------------------------------
ggplot(RF_SBeq_total_data) +
  geom_tile(aes(x = f, y = b, fill = rockfish_avg),width=0.012,height=0.006) +
  theme_classic()

ggplot(RF_CVeq_total_data) +
  geom_tile(aes(x = cv, y = autocorr_R, fill = rockfish_cv),width=0.015,height=0.015) +
  theme_classic()


ggplot(RF_SBeq_total_data) +
  geom_point(aes(x = f, y = b, color = rockfish_avg)) +
  theme_classic()


