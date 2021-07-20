library(tidyverse)
library(FSA)
library(FSAdata)
library(nlstools)
library(car)
library(AquaticLifeHistory)


# lingcod dataset: raw data from Bonnie et al. that only inclues stomach contents
# lingcod_empty dataset: raw data from Bonnie et al. that only includes lingcod with empty stomachs
# lingcod_full dataset: The binded dataset with lingcod and lingcod_empty
# lingcod_full_with_age dataset: The lingcod_full dataset with complete age information but only with identified sexes
# lingcod_rockfish dataset: lingcod_full_with_age dataset that focuses only on rockfish diet data, and include the 
#       total weight/number of stomach content, total weight/number of stomach content that is sebastes,
#       and fraction of total weight/number that is sebastes
# lingcod_rockfish_CA dataset: lingcod_rockfish dataset subsetted for ports off of CA


#################################
######### Load in Data ##########
#################################

### Read the data in
# With stomach contents
lingcod <- read.csv("raw_data/LingcodwithGutContents_All_BB.csv")
  # gotta fix some classifications...
  lingcod$Gape <- as.numeric(lingcod$Gape)
  lingcod$Gut.contents.kg <- as.numeric(lingcod$Gut.contents.kg)
  lingcod$Body.wt.minus.gut.contents <- as.numeric(lingcod$Body.wt.minus.gut.contents)
# Without stomach contents
lingcod_empty <- read.csv("raw_data/Lingcod_empty_stomachs_All_BB.csv")
  #gotta fix some classifications...
  lingcod_empty$logTL.cm <- as.numeric(lingcod_empty$logTL.cm)
  lingcod_empty$logWt.kg <- as.numeric(lingcod_empty$logWt.kg)
  lingcod_empty$Gape <- as.numeric(lingcod_empty$Gape)
# Bind the two datasets
lingcod_full <- bind_rows(lingcod, lingcod_empty)
  # save(lingcod_full, file = "cleaned_data/lingcod_full.Rdata")


#########################################
##### LETS ADD predicted ages!!! ########
#########################################

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


# Look at population model growth curve against this estimate growth curve from diet data
# first go to lingcod_dynamics script and run the length_l vector
plot(age, F_Linf*(1-exp(-F_k*age)), main = "Female growth curve", ylab = "length", type = "l", lty = 2) #estimated from diet data
lines(length_l[1,]) # population model
legend(x = "bottomright", legend = c("diet data", "stock assessment"), lty = c(2, 1))

plot(age, M_Linf*(1-exp(-M_k*age)), main = "Male growth curve", ylab = "length", type = "l", lty = 2) # estimated from diet data
lines(length_l[2,]) # population model
legend(x = "bottomright", legend = c("diet data", "stock assessment"), lty = c(2, 1))

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




#######################################
#### Let's add gut content info!! #####
#######################################

## Calculate the sum of all the gut content in weight and create new column
## with total gut content in weight
lingcod_gut_content_wt <- lingcod_full_with_age %>% 
  select(Sample.ID, contains("wt..")) %>% 
  filter(Sample.ID != "") # remove lingcod with missing sample.ID since these are empty stomachs anyways
lingcod_gut_content_wt[is.na(lingcod_gut_content_wt)] <- 0
lingcod_gut_content_wt <- lingcod_gut_content_wt %>% 
  mutate(wt..Total = rowSums(lingcod_gut_content_wt[,2:103])) %>% 
  select(Sample.ID, wt..Total)

## Calculate the sum of all the gut content by number (same as above but for X)
lingcod_gut_content_X <- lingcod_full_with_age %>% 
  select(Sample.ID, contains("X..")) %>% 
  filter(Sample.ID != "")
lingcod_gut_content_X[is.na(lingcod_gut_content_X)] <- 0
lingcod_gut_content_X <- lingcod_gut_content_X %>% 
  mutate(X..Total = rowSums(lingcod_gut_content_X[,2:103])) %>% 
  select(Sample.ID, X..Total)

# New working dataset with total food in stomach  
lingcod_working <- left_join(lingcod_full_with_age, lingcod_gut_content_wt, by = "Sample.ID")
lingcod_working <- left_join(lingcod_working, lingcod_gut_content_X, by = "Sample.ID")

# Create working dataset with relevant columns and make a
# new column grouping all Sebastes species into one
# Note - we are using the lingcod_full_with_age which subsets ONLY the identified male and female
#        but I checked the Sex == "" and they are all empty stomachs which shouldn't impact the diet composition stuff
lingcod_rockfish <- lingcod_working %>%
  select(Sample.ID, age_useful, Date..yymmdd., Port, Location, MPA.1.fish.take.prohibited., Lat, Long,
         Depth.ft, TL.cm, Sex, Wt.kg, Gape, Maturity, Stage, 
         Gut.wt.kg, Gut.contents.kg, X..Sebastes, wt..Sebastes, 
         X..Sebastes.jordani, wt..Sebastes.jordani, X..Yellowtail.Rockfish,
         wt..Yellowtail.Rockfish, X..Blue.Rockfish, wt..Blue.Rockfish,
         X..Black.Rockfish, wt..Black.RF, X..Puget.Sound.Rockfish,
         wt..Puget.Sound.Rockfish, wt..Total, X..Total) %>% 
  replace(is.na(.), 0) %>%
  mutate(X..Sebastes.Total = rowSums(select(.,X..Sebastes, X..Sebastes.jordani, # total number of sebastes contents
                                            X..Yellowtail.Rockfish, X..Blue.Rockfish,
                                            X..Black.Rockfish, X..Puget.Sound.Rockfish)),
         wt..Sebastes.Total = rowSums(select(.,wt..Sebastes, wt..Sebastes.jordani, # total weight of sebastes contents
                                             wt..Yellowtail.Rockfish, wt..Blue.Rockfish,
                                             wt..Black.RF, wt..Puget.Sound.Rockfish)),
         gut.ratio.sebastes.wt = wt..Sebastes.Total/wt..Total, # fraction of total weight in stomach that is sebastes
         gut.ratio.sebastes.X = X..Sebastes.Total/X..Total) # fraction of total number in stomach that is sebastes 

lingcod_rockfish$gut.ratio.sebastes.wt[is.na(lingcod_rockfish$gut.ratio.sebastes.wt)] <- 0
lingcod_rockfish$gut.ratio.sebastes.X[is.na(lingcod_rockfish$gut.ratio.sebastes.X)] <- 0

save(lingcod_rockfish, file = "cleaned_data/lingcod_rockfish.Rdata")

lingcod_rockfish_w_gutcontent = lingcod_rockfish %>% 
  filter(wt..Total > 0)

save(lingcod_rockfish_w_gutcontent, file = "cleaned_data/lingcod_rockfish_w_gutcontent.Rdata")

# Create new dataset with just lingcod caught from California ports
CA_ports <- c("MOR", "BOD", "SDG", "FTB", "SLO", "LOS", "EUR", "HMB", "MON", "EME", "SBA")
lingcod_rockfish_CA <- lingcod_rockfish %>% 
  filter(Port %in% CA_ports)

#save lingcod_rockfish_CA as data
save(lingcod_rockfish_CA, file = "cleaned_data/lingcod_rockfish_CA.Rdata")

  



