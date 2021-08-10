#############################
####### Age-Length Key ######
#############################

# These age-length keys are using only the observed ages (and not the estimated ages)
# This shows the number of fish within a length bin for each age. You can see the variation 
# of lengths for a given age as well as which lengths are most common for each age
library(tidyverse)
library(FSA)
load(file = "cleaned_data/lingcod_full.Rdata")

# organize what I need by changing the sex U to F and filtering dataset to only have F and M
lingcod_full <- lingcod_full %>% 
  mutate(Sex = ifelse(Sex == "U", Sex.1, Sex)) %>% # I changed the U in the Sex column to the value in the Sex.1 columm
  filter(Sex == "F" | Sex == "M")

# Create length bins
bin_size = 5
round_to = 5
bins <- seq(plyr::round_any(min(lingcod_full$TL.cm), round_to, f = floor), plyr::round_any(max(lingcod_full$TL.cm), round_to, f = ceiling), by= bin_size)
TL.cm <- as.numeric(lingcod_full$TL.cm)
rangelabels <- paste(head(bins,-1), tail(bins,-1), sep="-")
lingcod_full$Bin <- cut(TL.cm, bins, rangelabels)


age_length_F = lingcod_full %>% 
  group_by(Bin, Ages, Sex) %>% 
  filter(Sex == "F") %>% 
  summarise(n=n()) %>% 
  dplyr::select(Bin, Ages, n) %>% 
  drop_na() %>% 
  spread(Ages, n)

age_length_M = lingcod_full %>% 
  group_by(Bin, Ages, Sex) %>% 
  filter(Sex == "M") %>% 
  summarise(n=n()) %>% 
  dplyr::select(Bin, Ages, n) %>% 
  drop_na() %>% 
  spread(Ages, n)

# Let's visualize!
#females
age_length_F$Bin = gsub("\\-.*","\\", age_length_F$Bin)
age_length_F = as.data.frame(age_length_F)
age_length_F_2 <- age_length_F[,-1]
rownames(age_length_F_2) <- age_length_F[,1]
age_length_F = age_length_F_2
rm(age_length_F_2)
age_length_F[is.na(age_length_F)] <- 0
age_length_F = as.matrix(age_length_F)
age_length_F_proportion = prop.table(age_length_F, 2)
alkPlot(age_length_F_proportion,"splines")
#males
age_length_M$Bin = gsub("\\-.*","\\", age_length_M$Bin)
age_length_M = as.data.frame(age_length_M)
age_length_M_2 <- age_length_M[,-1]
rownames(age_length_M_2) <- age_length_M[,1]
age_length_M = age_length_M_2
rm(age_length_M_2)
age_length_M[is.na(age_length_M)] <- 0
age_length_M = as.matrix(age_length_M)
age_length_M_proportion = prop.table(age_length_M, 2)
alkPlot(age_length_M_proportion,"splines")

# I want to see the range of lengths for each age. So, I extract the first and last length bin
# with a non-zero value for each age class
# Couldn't figure out the code to extract the first and last non-zero index, and so I just did
# it by hand.
F_length_intervals = bind_rows(initial = c(15,30,30,35,55,50,65,65,55,85,75,85,100,95,105,110), 
                               final = c(30,55,75,90,95,95,105,100,105,110,110,115,105,115,115,115)) %>% 
  mutate(age = c(1:13,15:17))
# The data is limited after age 12, so it might be worth aggregating those ages
M_length_intervals = bind_rows(initial = c(15,30,30,40,45,35,55,55,65,70,65,70,70,70), 
                               final = c(35,70,80,90,85,85,85,90,90,95,90,85,95,75)) %>% 
  mutate(age = 1:14)
# The data is limited after age 11, so it might be worth aggregating those ages
plot(F_length_intervals$age, F_length_intervals$final, type = "p", ylim = c(0,120))
points(F_length_intervals$age, F_length_intervals$initial)
plot(M_length_intervals$age, M_length_intervals$final, type = "p", ylim = c(0,120))
points(M_length_intervals$age, M_length_intervals$initial)

# Let's fit a gam into these upper and lower interval length bins to calculate the variation
# in lengths at each age
F_gam_length_UI = gam(final ~ s(age), data = F_length_intervals, family = gaussian(), method = "ML")
plot(F_gam_length_UI,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
F_gam_length_LI = gam(initial ~ s(age), data = F_length_intervals, family = gaussian(), method = "ML")
plot(F_gam_length_LI,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)

M_gam_length_UI = gam(final ~ s(age), data = M_length_intervals, family = gaussian(), method = "ML")
plot(M_gam_length_UI,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
M_gam_length_LI = gam(initial ~ s(age), data = M_length_intervals, family = gaussian(), method = "ML")
plot(M_gam_length_LI,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)


# Now, let's use the gam to predict the upper and lower lengths for each age, then 
# take the difference of them to find the length variation for each age
F_age_to_predict = as.data.frame(1:12) # 12 is the max age after which we aggregate due to limited data
colnames(F_age_to_predict) = "age"
F_UI_prediction = predict.gam(F_gam_length_UI, F_age_to_predict)
F_LI_prediction = predict.gam(F_gam_length_LI, F_age_to_predict)
F_length_variation = (F_UI_prediction - F_LI_prediction)/2 # Divide by 2 to get variation only going 1 way

M_age_to_predict = as.data.frame(1:11) # 11 is the max age after which we aggregate due to limited data
colnames(M_age_to_predict) = "age"
M_UI_prediction = predict.gam(M_gam_length_UI, M_age_to_predict)
M_LI_prediction = predict.gam(M_gam_length_LI, M_age_to_predict)
M_length_variation = (M_UI_prediction - M_LI_prediction)/2


# save(M_length_variation, file = "cleaned_data/M_length_variation.Rdata")
# save(F_length_variation, file = "cleaned_data/F_length_variation.Rdata")
# save(age_length_F_proportion, file = "cleaned_data/age_length_F_proportion.Rdata")
# save(age_length_M_proportion, file = "cleaned_data/age_length_M_proportion.Rdata")
# save(age_length_F, file = "cleaned_data/age_length_F.Rdata")
# save(age_length_M, file = "cleaned_data/age_length_M.Rdata")



