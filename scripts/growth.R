library(tidyverse)
library(AquaticLifeHistory) # growth package

# This script adds the predicted ages into the lingcod_full dataset and re-saves this dataset.

# Load in full dataset
load(file = "cleaned_data/lingcod_full.Rdata")

# organize what I need
lingcod_full <- lingcod_full %>% 
  mutate(Sex = ifelse(Sex == "U", Sex.1, Sex)) %>% # I changed the U in the Sex column to the value in the Sex.1 columm
  filter(Sex == "F" | Sex == "M") # Then I only want males and females 


# Calculate growth using a package (AquaticLifeHistory) https://github.com/jonathansmart/AquaticLifeHistory
growth_data = lingcod_full %>% 
  select(Ages, TL.cm, Sex)
colnames(growth_data) = c("Age", "Length", "Sex")
growth_data_F = growth_data %>% 
  filter(Sex == "F")
growth_data_M = growth_data %>% 
  filter(Sex == "M")
F_growthparms = Estimate_Growth(growth_data_F)
F_Linf = F_growthparms$VonB[1,1]
F_k = F_growthparms$VonB[2,1]
# Linf = 108.56   se = 2.542
# k = 0.182     se = 0.014
M_growthparms = Estimate_Growth(growth_data_M)
M_Linf = M_growthparms$VonB[1,1]
M_k = M_growthparms$VonB[2,1]
# Linf = 84.557   se = 2.461
# k = 0.192     se = 0.021

# Create a function that predicts age given lengths
#females
predict_age_f = function(length) {
  if(length >= 106) { # I chose this length because that's when the predicted age goes above 20
    age_pred = 20 # set max age to 20
  } else {
    age_pred = log(1-length/F_Linf) / -F_k
  }
  return(age_pred)
}
#males
predict_age_m = function(length) {
  if(length >= 83) { # I chose this length because that's when the predicted age goes above 20
    age_pred = 20 # set max age to 20
  } else {
    age_pred = log(1-length/M_Linf) / -M_k
  }
  return(age_pred)
}

# Vectorize to prevent warning messages
predict_age_m <- Vectorize(predict_age_m)
predict_age_f <- Vectorize(predict_age_f)

# Predict ages for males and females and add to the dataset
lingcod_full_F = lingcod_full %>% #female
  filter(Sex == "F") %>% 
  mutate(age_pred = round(predict_age_f(TL.cm)))
lingcod_full_M = lingcod_full %>% #male
  filter(Sex == "M") %>% 
  mutate(age_pred = round(predict_age_m(TL.cm)))
lingcod_full = bind_rows(lingcod_full_F, lingcod_full_M) %>% # bind rows
  mutate(age_error = Ages - age_pred) # Determine error

ggplot(lingcod_full, aes(TL.cm, age_error)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(70,110), linetype='dotted', col = 'red') +
  labs(y = "Error", x = "Length (cm)", title = "Error in Predicted Lengths")

age_NA <- lingcod_full[is.na(lingcod_full$Ages),]

ggplot(lingcod_full, aes(TL.cm)) +
  geom_histogram() +
  theme_classic() +
  geom_vline(xintercept = c(70,110), linetype='dotted', col = 'red') +
  labs(title = "Lengths to be predicted", x = "Length (cm)")

### Create a new column with complete age (NA's filled in with predicted)

lingcod_full_with_age = lingcod_full %>% 
  mutate(age_useful = Ages)
lingcod_full_with_age$age_useful[is.na(lingcod_full_with_age$age_useful)] <- lingcod_full_with_age$age_pred[is.na(lingcod_full_with_age$age_useful)]


#age_check = lingcod_full_with_age %>% 
#  select(Sex, TL.cm, Ages, age_pred, age_error, age_useful)

save(lingcod_full_with_age, file = "cleaned_data/lingcod_full_with_age.Rdata")


####################################
####################################
####### Predict Age by Length ######
####################################
####################################


# Create a function that predicts age using the francis parameterization of the 
# VBGF. At ages that are too high, the function fails because it attempts to 
# take a log of a negative number which isn't possible. So, an if statement
# is used to set a specific age for lengths that are too large. This max age is 
# is determined by the max age in the observed data.

# Here I use the lingcod_full dataset which has all lingcod with and without stomach contents.

# MALES
ages_francis_male <- c(3,11) # male t1 and t3 age
male_TL_aget1 = lingcod_full %>% # Filter out all male age at t1
  filter(Sex == "M", Ages == ages_francis_male[1]) %>% 
  select(TL.cm)
male_TL_aget2 = lingcod_full %>% # Filter out all male at t2
  filter(Sex == "M", Ages == mean(ages_francis_male)) %>% 
  select(TL.cm)
male_TL_aget3 = lingcod_full %>% # Filter out all male at t3
  filter(Sex == "M", Ages == ages_francis_male[2]) %>% 
  select(TL.cm)
L1_m <- mean(male_TL_aget1$TL.cm) # mean length at relatively young age t1 (age 3 used here)
L2_m <- mean(male_TL_aget2$TL.cm) # mean length at the average of t1 and t2 (age 7 used here)
L3_m <- mean(male_TL_aget3$TL.cm) # mean length at relatively old age t2 (age 11 used here)
r_m <- (L3_m-L2_m)/(L2_m-L1_m)

predict_age_m <- function(length) { # Function that predicts male age given length
  if(length >= 85) { # I chose this length because that's when the predicted age goes above 20
    pred_age <-  20 # set max age to 20
  } else {
    pred_age <- ages_francis_male[1] + (ages_francis_male[2]-ages_francis_male[1])*((log(1-(1-r_m^2)*((length-L1_m)/(L3_m-L1_m))))/(2*log(r_m)))
  }
  return(pred_age) 
}

# FEMALES
ages_francis_female <- c(3,15)
female_TL_aget1 = lingcod_full %>% # Filter out all female age at t1
  filter(Sex == "F", Ages == ages_francis_female[1]) %>% 
  select(TL.cm)
female_TL_aget2 = lingcod_full %>% # Filter out all female age at t2
  filter(Sex == "F", Ages == mean(ages_francis_female)) %>% 
  select(TL.cm)
female_TL_aget3 = lingcod_full %>% # Filter out all female age at t2
  filter(Sex == "F", Ages == ages_francis_female[2]) %>% 
  select(TL.cm)
L1_f <- mean(female_TL_aget1$TL.cm)  # mean length at relatively young age t1 (age 3 used here)
L2_f <- mean(female_TL_aget2$TL.cm)  # mean length at the average of t1 and t2 (age 9 used here)
L3_f <- mean(female_TL_aget3$TL.cm) # mean length at relatively old age t2 (age 15 used here)
r_f <- (L3_f-L2_f)/(L2_f-L1_f)

predict_age_f <- function(length) {
  if(length >= 117) {
    pred_age <-  20
  } else {
    pred_age <- ages_francis_female[1] + (ages_francis_female[2]-ages_francis_female[1])*((log(1-(1-r_f^2)*((length-L1_f)/(L3_f-L1_f))))/(2*log(r_f)))
  }
  return(pred_age) 
}

# this helps prevent the warning messages.
v_predict_age_m <- Vectorize(predict_age_m)
v_predict_age_f <- Vectorize(predict_age_f)


# making a predicted age x length dataframe to graph over the data

lingcod_male <- lingcod_full %>% 
  filter(Sex == "M")
lingcod_female <- lingcod_full %>% 
  filter(Sex == "F")

#males
predicted_age_m = sapply(1:120, predict_age_m)
predicted_age_m = pmax(predicted_age_m, 0)
predicted_age_m = as.data.frame(predicted_age_m) %>% 
  mutate(length = 1:120)
#females
predicted_age_f = sapply(1:120, predict_age_f)
predicted_age_f = pmax(predicted_age_f, 0)
predicted_age_f = as.data.frame(predicted_age_f) %>% 
  mutate(length = 1:120)

plot(lingcod_female$Ages, lingcod_female$TL.cm, type = "p", xlab = "Ages", ylab = "Length (cm)", main = "female", ylim = c(0,120))
lines(predicted_age_f$predicted_age_f, predicted_age_f$length)

plot(lingcod_male$Ages, lingcod_male$TL.cm, type = "p", xlab = "Ages", ylab = "Length (cm)", main = "male", ylim = c(0,100))
lines(predicted_age_m$predicted_age_m, predicted_age_m$length)


# Add a new row with the predicted ages
## Males
lingcod_male <- lingcod_full %>% 
  filter(Sex == "M") %>% 
  mutate(age_pred = v_predict_age_m(TL.cm))

## Females

lingcod_female <- lingcod_full %>% 
  filter(Sex == "F") %>% 
  mutate(age_pred = v_predict_age_f(TL.cm))

#predict_age = function(length, sex) {
#  if (sex == "F") {
#    age = v_predict_age_f(length)
#    print(age)
#  } else if (sex == "M") {
#    age = v_predict_age_m(length)
#    print(age)
#  } else
#    print(NA)
#}

#lingcod_full <- lingcod_full %>% 
#  mutate(age_pred = predict_age(TL.cm, Sex))

# Let's compare how well the predicted ages align with the observed ages
lingcod_full = bind_rows(lingcod_male, lingcod_female)
lingcod_full$age_pred = round(lingcod_full$age_pred)

lingcod_agecheck = lingcod_full %>% 
  select(Sex, TL.cm, Ages, age_pred) %>% 
  mutate(underestimate = Ages > age_pred,
         overestimate = Ages < age_pred,
         correct = Ages == age_pred)








