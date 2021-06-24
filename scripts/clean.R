library(tidyverse)
library(FSA)
library(FSAdata)
library(nlstools)
library(car)


####################################
####################################
######### Setting up Data ##########
####################################
####################################

### Read the data in
lingcod <- read.csv("data/LingcodwithGutContents_All_BB.csv")

## Calculate the sum of all the gut content in weight and create new column
## with total gut content in weight
lingcod_gut_content_wt <- lingcod %>% 
  select(Sample.ID, contains("wt.."))
lingcod_gut_content_wt[is.na(lingcod_gut_content_wt)] <- 0
lingcod_gut_content_wt <- lingcod_gut_content_wt %>% 
  mutate(wt..Total = rowSums(lingcod_gut_content_wt[,2:103])) %>% 
  select(Sample.ID, wt..Total)

## Calculate the sum of all the gut content by number
lingcod_gut_content_X <- lingcod %>% 
  select(Sample.ID, contains("X.."))
lingcod_gut_content_X[is.na(lingcod_gut_content_X)] <- 0
lingcod_gut_content_X <- lingcod_gut_content_X %>% 
  mutate(X..Total = rowSums(lingcod_gut_content_X[,2:103])) %>% 
  select(Sample.ID, X..Total)

# New working dataset with total food in stomach  
lingcod_working <- left_join(lingcod, lingcod_gut_content_wt, by = "Sample.ID")
lingcod_working <- left_join(lingcod_working, lingcod_gut_content_X, by = "Sample.ID")
lingcod_working <- lingcod_working %>% # Remove the NA rows at the end
  slice(1:1321)

# Create working dataset with relevant columns and make a
# new column grouping all Sebastes species into one
lingcod_rockfish <- lingcod_working %>%
  select(Sample.ID, Ages, Date..yymmdd., Port, Location, MPA.1.fish.take.prohibited., Lat, Long,
         Depth.ft, TL.cm, Sex.1, Wt.kg, Gape, Maturity, Stage, 
         Gut.wt.kg, Gut.contents.kg, X..Sebastes, wt..Sebastes, 
         X..Sebastes.jordani, wt..Sebastes.jordani, X..Yellowtail.Rockfish,
         wt..Yellowtail.Rockfish, X..Blue.Rockfish, wt..Blue.Rockfish,
         X..Black.Rockfish, wt..Black.RF, X..Puget.Sound.Rockfish,
         wt..Puget.Sound.Rockfish, wt..Total, X..Total) %>% 
  replace(is.na(.), 0) %>%
  mutate(X..Sebastes.Total = rowSums(select(.,X..Sebastes, X..Sebastes.jordani,
                                            X..Yellowtail.Rockfish, X..Blue.Rockfish,
                                            X..Black.Rockfish, X..Puget.Sound.Rockfish))) %>% 
  mutate(wt..Sebastes.Total = rowSums(select(.,wt..Sebastes, wt..Sebastes.jordani,
                                             wt..Yellowtail.Rockfish, wt..Blue.Rockfish,
                                             wt..Black.RF, wt..Puget.Sound.Rockfish))) %>% 
  filter(Sex.1 != "") ## We realized that rows with NA for Sex were all just 0s

# Replace all 0s in Ages column with NAs
lingcod_rockfish$Ages <- replace(lingcod_rockfish$Ages, lingcod_rockfish$Ages == 0, NA)


####################################
####################################
##### Exploratory Data Analysis ####
####################################
####################################


# Ages x Total Length (colored by sex)
ggplot(lingcod_rockfish, aes(x = Ages, y = TL.cm, color = Sex.1)) +
  geom_point() +
  theme_classic()

# Ages x total Sebastes weight (colored by sex)

SebastesWeightRangeByAge <- lingcod_rockfish %>% 
  group_by(Ages, Sex.1) %>%
  summarise(
    MaxSebastesWByAge = max(wt..Sebastes.Total, na.rm = T),
    MeanSebastesWByAge = mean(wt..Sebastes.Total, na.rm = T),
    SumSebastesWByAge = sum(wt..Sebastes.Total, na.rm = T),
    n()
  )

ggplot(SebastesWeightRangeByAge, aes(Ages, SumSebastesWByAge, fill = Sex.1)) +
  geom_col() +
  theme_classic()

# Age x Depth (colored by sex)
ggplot(lingcod_rockfish, aes(x = Ages, y = Depth.ft, color = Sex.1)) +
  geom_point() +
  theme_classic()

# Age x Lingcod weight (colored by sex)
ggplot(lingcod_rockfish, aes(x = Ages, y = Wt.kg, color = Sex.1)) +
  geom_point() +
  theme_classic()

### How many samples do we have for males and females for each age class
FreqOfSexByAge <- lingcod_rockfish %>% 
  group_by(Ages, Sex.1) %>%
  summarise(
    n = n()
  ) %>%
  arrange(Ages)

# Lets look at how many observations there are that actually contain sebastes in the stomach
lingcodWithSebastes <- lingcod_rockfish %>%
  filter(wt..Sebastes.Total>0) %>% 
  tally() # there are 165 observations with sebastes in stomach



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

# MALES
ages_francis_male <- c(3,11) # male t1 and t3 age
male_TL_aget1 = lingcod_rockfish %>% # Filter out all male age at t1
  filter(Sex.1 == "M", Ages == ages_francis_male[1]) %>% 
  select(TL.cm)
male_TL_aget2 = lingcod_rockfish %>% # Filter out all male at t2
  filter(Sex.1 == "M", Ages == mean(ages_francis_male)) %>% 
  select(TL.cm)
male_TL_aget3 = lingcod_rockfish %>% # Filter out all male at t3
  filter(Sex.1 == "M", Ages == ages_francis_male[2]) %>% 
  select(TL.cm)
L1_m <- mean(male_TL_aget1$TL.cm) # mean length at relatively young age t1 (age 3 used here)
L2_m <- mean(male_TL_aget2$TL.cm) # mean length at the average of t1 and t2 (age 7 used here)
L3_m <- mean(male_TL_aget3$TL.cm) # mean length at relatively old age t2 (age 11 used here)
r_m <- (L3_m-L2_m)/(L2_m-L1_m)

predict_age_m <- function(length) { # Function that predicts male age given length
  if(length >= 83.8) {
    pred_age <-  14
  } else {
    pred_age <- ages_francis_male[1] + (ages_francis_male[2]-ages_francis_male[1])*((log(1-(1-r_m^2)*((length-L1_m)/(L3_m-L1_m))))/(2*log(r_m)))
  }
  return(pred_age) 
}

# FEMALES
ages_francis_female <- c(3,15)
female_TL_aget1 = lingcod_rockfish %>% # Filter out all female age at t1
  filter(Sex.1 == "F", Ages == ages_francis_female[1]) %>% 
  select(TL.cm)
female_TL_aget2 = lingcod_rockfish %>% # Filter out all female age at t2
  filter(Sex.1 == "F", Ages == mean(ages_francis_female)) %>% 
  select(TL.cm)
female_TL_aget3 = lingcod_rockfish %>% # Filter out all female age at t2
  filter(Sex.1 == "F", Ages == ages_francis_female[2]) %>% 
  select(TL.cm)
L1_f <- mean(female_TL_aget1$TL.cm)  # mean length at relatively young age t1 (age 3 used here)
L2_f <- mean(female_TL_aget2$TL.cm)  # mean length at the average of t1 and t2 (age 9 used here)
L3_f <- mean(female_TL_aget3$TL.cm) # mean length at relatively old age t2 (age 15 used here)
r_f <- (L3_f-L2_f)/(L2_f-L1_f)

predict_age_f <- function(length) {
  if(length >= 105) {
    pred_age <-  17
  } else {
    pred_age <- ages_francis_female[1] + (ages_francis_female[2]-ages_francis_female[1])*((log(1-(1-r_f^2)*((length-L1_f)/(L3_f-L1_f))))/(2*log(r_f)))
  }
  return(pred_age) 
}

# this helps prevent the warning messages.
v_predict_age_m <- Vectorize(predict_age_m)
v_predict_age_f <- Vectorize(predict_age_f)


# Add a new row with the predicted ages
## Males
lingcod_male <- lingcod_rockfish %>% 
  filter(Sex.1 == "M")

lingcod_male <- lingcod_male %>% 
  mutate(age_pred = v_predict_age_m(lingcod_male$TL.cm))

## Females

lingcod_female <- lingcod_rockfish %>% 
  filter(Sex.1 == "F")

lingcod_female <- lingcod_female %>% 
  mutate(age_pred = v_predict_age_f(lingcod_female$TL.cm))


####################################
####################################
##### Cleaning Updated  Dataset ####
####################################
####################################

### Let's remove the decimal places for the predicted age
lingcod_male$age_pred <- round(lingcod_male$age_pred)
lingcod_female$age_pred <- round(lingcod_female$age_pred)

### Create a new column with complete age (NA's filled in with predicted)
lingcod_male <- lingcod_male %>% 
  mutate(age_useful = Ages)

lingcod_female <- lingcod_female %>% 
  mutate(age_useful = Ages)

lingcod_male$age_useful[is.na(lingcod_male$age_useful)] <- lingcod_male$age_pred[is.na(lingcod_male$age_useful)]
lingcod_female$age_useful[is.na(lingcod_female$age_useful)] <- lingcod_female$age_pred[is.na(lingcod_female$age_useful)]

### Merge male and female dataset to re-create the lingcod_rockfish dataset
lingcod_rockfish <- rbind(lingcod_male, lingcod_female)

### Create a new column with the proportion of stomach content that is sebastes
lingcod_rockfish <- lingcod_rockfish %>% 
  mutate(gut.ratio.sebastes.wt = wt..Sebastes.Total/wt..Total) %>% 
  mutate(gut.ratio.sebastes.X = X..Sebastes.Total/X..Total)
lingcod_rockfish$gut.ratio.sebastes.wt[is.na(lingcod_rockfish$gut.ratio.sebastes.wt)] <- 0

# Create new dataset with just lingcod caught from California ports
CA_ports <- c("MOR", "BOD", "SDG", "FTB", "SLO", "LOS", "EUR", "HMB", "MON", "EME", "SBA")
lingcod_rockfish_CA <- lingcod_rockfish %>% 
  filter(Port %in% CA_ports)

####################################
####################################
##### Diet Fraction Calculation ####
####################################
####################################

# Dataframe with the average diet fraction by weight that is sebastes for each age and sex 
dietfrac_by_age_wt <- lingcod_rockfish_CA %>% 
  group_by(Sex.1, age_useful) %>% 
  summarise(mean = mean(gut.ratio.sebastes.wt), sd = sd(gut.ratio.sebastes.wt), 
            n = n(), max = max(gut.ratio.sebastes.wt), min = min(gut.ratio.sebastes.wt))
dietfrac_by_age_wt$sd[is.na(dietfrac_by_age_wt$sd)] <- 0

# Dataframe with the average diet fraction by number that is sebastes for each age and sex 
dietfrac_by_age_X <- lingcod_rockfish_CA %>% 
  group_by(Sex.1, age_useful) %>% 
  summarise(mean = mean(gut.ratio.sebastes.X), sd = sd(gut.ratio.sebastes.X), n = n(),
            max = max(gut.ratio.sebastes.X), min = min(gut.ratio.sebastes.X))
dietfrac_by_age_X$sd[is.na(dietfrac_by_age_X$sd)] <- 0

# The proportion of stomachs that contained rockfish of all the stomachs
lingcod_rockfish_CA_F <- lingcod_rockfish_CA %>% 
  filter(Sex.1 == "F")
lingcod_rockfish_CA_M <- lingcod_rockfish_CA %>% 
  filter(Sex.1 == "M")
perc.occurrence <- (sum(lingcod_rockfish_CA$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA))*100
perc.occurrence.F <- (sum(lingcod_rockfish_CA_F$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA_F))*100
perc.occurrence.M <- (sum(lingcod_rockfish_CA_M$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA_M))*100

perc.abundance.number <- sum(lingcod_rockfish_CA$X..Sebastes.Total)/sum(lingcod_rockfish_CA$X..Sebastes.Total > 0)
perc.abundance.weight <- sum(lingcod_rockfish_CA$wt..Sebastes.Total)/sum(lingcod_rockfish_CA$X..Sebastes.Total > 0)

### Proportion of gut content that is Sebastes across Ages
  # males
ggplot(lingcod_rockfish_CA[which(lingcod_rockfish_CA$Sex.1 == "M"),], aes(age_useful, gut.ratio.sebastes.wt)) +
  geom_point(color = "cyan3") +
  theme_classic() +
  labs( x = "Lingcod Age", y = "Fraction of diet that is Sebastes", title = "Male") +
  xlim(0,17)
  # females
ggplot(lingcod_rockfish_CA[which(lingcod_rockfish_CA$Sex.1 == "F"),], aes(age_useful, gut.ratio.sebastes.wt)) +
  geom_point(color = "lightcoral") +
  theme_classic() +
  labs( x = "Lingcod Age", y = "Fraction of diet that is Sebastes", title = "Female") +
  xlim(0,17)

#graphing the mean diet fraction by age
ggplot(dietfrac_by_age_wt, aes(age_useful, mean, col = Sex.1)) +
  geom_point() +
  labs( x = "Lingcod Age", y = "Average fraction of diet that is Sebastes") +
  theme_classic()

####################################
####################################
##### Exploratory Data Analysis ####
####################################
####################################

### Sebastes consumption by length (colored by sex)
    # Including the zeros
ggplot(lingcod_rockfish, aes(TL.cm, wt..Sebastes.Total)) +
  geom_point(aes(color = Sex.1)) +
  theme_classic() +
  labs( x = "Lingcod Length (cm)", y = "Sebastes consumption") 
    # Excluding the zeros
ggplot(lingcod_rockfish[which(lingcod_rockfish$wt..Sebastes.Total > 0),], aes(TL.cm, wt..Sebastes.Total)) +
  geom_point(aes(color = Sex.1)) +
  theme_classic() +
  labs( x = "Lingcod Length (cm)", y = "Sebastes consumption") 


### Proportion of gut content that is Sebastes across Ages
ggplot(lingcod_rockfish[which(lingcod_rockfish$wt..Total>0),], aes(age_useful, gut.ratio.sebastes.wt)) +
  geom_point(aes(color = Sex.1)) +
  theme_classic() +
  labs( x = "Lingcod Age", y = "Proportion consumed that is Sebastes") 
    # Just males
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "M"),], aes(age_useful, gut.ratio.sebastes.wt)) +
  geom_point(color = "cyan3") +
  theme_classic() +
  labs( x = "Lingcod Age", y = "Proportion consumed that is Sebastes") +
  xlim(0,17)
    # Just females
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "F"),], aes(age_useful, gut.ratio.sebastes.wt)) +
  geom_point(color = "lightcoral") +
  theme_classic() +
  labs( x = "Lingcod Age", y = "Proportion consumed that is Sebastes") +
  xlim(0,17)

### Proportion of gut content that is Sebastes across Sizes
ggplot(lingcod_rockfish[which(lingcod_rockfish$wt..Total>0),], aes(TL.cm, gut.ratio.sebastes.wt)) +
  geom_point(aes(color = Sex.1)) +
  theme_classic() +
  labs( x = "Lingcod Size", y = "Proportion consumed that is Sebastes")
    #Just males
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "M"),], aes(TL.cm, gut.ratio.sebastes.wt)) +
  geom_point(color = "cyan3") +
  theme_classic() +
  labs( x = "Lingcod Size", y = "Proportion consumed that is Sebastes") +
  xlim(25,120)
    #Just females
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "F"),], aes(TL.cm, gut.ratio.sebastes.wt)) +
  geom_point(color = "lightcoral") +
  theme_classic() +
  labs( x = "Lingcod Size", y = "Proportion consumed that is Sebastes") +
  xlim(25,120)


### Proportion of gut content that is Sebastes across depths (Sex done separately)
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "M"),], aes(Depth.ft, gut.ratio.sebastes.wt)) +
  geom_point(color = "cyan3") +
  theme_classic() +
  labs( x = "Depth", y = "Proportion consumed that is Sebastes") +
  xlim(0,600)

ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "F"),], aes(Depth.ft, gut.ratio.sebastes.wt)) +
  geom_point(color = "lightcoral") +
  theme_classic() +
  labs( x = "Depth", y = "Proportion consumed that is Sebastes") +
  xlim(0,600)


### Sebastes consumption in MPA vs fishing ground
ggplot(lingcod_rockfish, aes(MPA.1.fish.take.prohibited., wt..Sebastes.Total)) +
  geom_point() +
  theme_classic()

### Total consumption in MPA vs fishing ground
ggplot(lingcod_rockfish, aes(MPA.1.fish.take.prohibited., wt..Total)) +
  geom_point() +
  theme_classic()

### Proportion of gut that is sebastes in MPA vs fishing ground
ggplot(lingcod_rockfish, aes(MPA.1.fish.take.prohibited., gut.ratio.sebastes.wt)) +
  geom_point(aes(color = Sex.1)) +
  theme_classic()
    # With Female
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "F"),], aes(MPA.1.fish.take.prohibited., gut.ratio.sebastes.wt)) +
  geom_point(color = "lightcoral") +
  theme_classic()
    # With Males
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "M"),], aes(MPA.1.fish.take.prohibited., gut.ratio.sebastes.wt)) +
  geom_point(color = "cyan3") +
  theme_classic()


### How does length range vary inside and outside of an MPA?
lingcod_rockfish$MPA.1.fish.take.prohibited. <- factor(lingcod_rockfish$MPA.1.fish.take.prohibited.)
ggplot(lingcod_rockfish[which(lingcod_rockfish$Sex.1 == "M"),], aes(MPA.1.fish.take.prohibited., TL.cm, fill = MPA.1.fish.take.prohibited.)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none")
  



