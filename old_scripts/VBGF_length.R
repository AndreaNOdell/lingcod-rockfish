library(tidyverse)
library(FSA)
library(FSAdata)
library(nlstools)
library(car)


# Males
lingcoddiets_male <- lingcod_rockfish %>% 
  filter(Sex.1 == "M") %>% 
  select(Sample.ID, Ages, TL.cm) %>%  ### 690 observations prior to omitting NA's
  na.omit()   ### 551 complete observations (NA's removed)

count(lingcoddiets_male)

# Females
lingcoddiets_female <- lingcod_rockfish %>% 
  filter(Sex.1 == "F") %>% 
  select(Sample.ID, Ages, TL.cm) %>%  ### 630 observations prior to omitting NA's
  na.omit()   ### 580 complete observations (NA's removed)

count(lingcoddiets_female)


# Visualize the data
ggplot(lingcoddiets_female, aes(x = Ages, y = TL.cm)) +
  geom_point() +
  theme_classic()

ggplot(lingcoddiets_male, aes(x = Ages, y = TL.cm)) +
  geom_point() +
  theme_classic()


###
### Getting starting values for parameters
###

  # Males
sv_francis_l_males <- vbStarts(Ages~TL.cm,data=lingcoddiets_male,type="Francis")
names(sv_francis_l_males) <- c("t1", "t2", "t3")
unlist(sv_francis_l_males)
    ### T1 = 1, T2 = 7, T3 = 10

  # Females
sv_francis_l_females <- vbStarts(Ages~TL.cm,data=lingcoddiets_female,type="Francis",tFrancis=c(2,17))
names(sv_francis_l_females) <- c("t1", "t2", "t3")
unlist(sv_francis_l_females)
    ### T1 = 2, T2 = 6.833333, T3 = 15.000000


###
### Call in the Growth function
###

vbfrancis <- vbFuns("Francis")
ages_francis_male <- c(2,12)
ages_francis_female <- c(2,17)

###
### Fit growth model to the data
###

# male

fittedgrowth_francis_male <- nls(TL.cm~vbfrancis(Ages,L1,L2,L3,t1=ages_francis_male[1],t3=ages_francis_male[2]), 
                                 data=lingcoddiets_male,start=sv_francis_males)
overview(fittedgrowth_francis_male)
bootfrancis_male <- nlsBoot(fittedgrowth_francis_male,niter=1000)
confint(bootfrancis_male,plot=TRUE)
plot(bootfrancis_male)
  # a little better but still somewhat correlated

#female
fittedgrowth_francis_female <- nls(TL.cm~vbfrancis(Ages,L1,L2,L3,t1=ages_francis_female[1],t3=ages_francis_female[2]), 
                                   data=lingcoddiets_female,start=sv_francis_females)
overview(fittedgrowth_francis_female)
bootfrancis_female <- nlsBoot(fittedgrowth_francis_female,niter=500)
confint(bootfrancis_female,plot=TRUE)
plot(bootfrancis_female)
# A little better, but still somewhat correlated


###
### Quick visualization
###

fitPlot(fittedgrowth_francis_male,xlab="Age",ylab="Total Length (mm)",main="male")
fitPlot(fittedgrowth_francis_female,xlab="Age",ylab="Total Length (mm)",main="female")




### Welp, the equation is predicting length, not age. So let's do some algebra 
### and have it predict Age

#male

##### The males length ranges from 27 through 90.4, so I will create a dataframe that
##### predicts from 32.9 (because anything smaller is predicting a negative number) through 83.8 (because the oldest recorded age 13 is predicted to be 83.7
##### and the predictions for anything above 83.7 may not be biologically sound)
LengthPredicted_male <- data.frame(TL.cm = 62.3)
predict(fittedgrowth_francis_male, LengthPredicted_male)

ests_m <- bootfrancis_male$coefboot
L1_m <- 47.1 #ests_m[,"L1"]
L2_m <- 69.6 #ests_m[,"L2"]
L3_m <- 81 #ests_m[,"L3"]
r_m <- (L3_m-L2_m)/(L2_m-L1_m)

      # Predicting
male_lengthrange <- seq(from = 32.9, to = 83.7, by = 0.1)
pva_m <- numeric(length(male_lengthrange))
for (i in 1:length(male_lengthrange)) {
  pva_m[i] <- ages_francis_male[1] + (ages_francis_male[2]-ages_francis_male[1])*((log(1-(1-r_m^2)*((male_lengthrange[i]-L1_m)/(L3_m-L1_m))))/(2*log(r_m)))
}

## convert to dataframe and add a column with associated lengths
pva_m <- as.data.frame(pva_m)
pva_m <- mutate(pva_m, length = seq(from = 32.9, to = 83.7, by = 0.1))

## Since we only predicted until 83.7, but there is length's up to 90.4, I will create a new dataframe
## to fill in the rest of the data with 83.8-90.4 all "predicting" age 14
seq_83 <- seq(from = 83.8, to = 90.4, by = 0.1)
pva_m_83up <- cbind(rep(14, times = 67), seq_83)
colnames(pva_m_83up) <- c("pva_m", "length")

###
### This is the final dataframe with predicted vals for length
###
pva_m <- rbind(pva_m, pva_m_83up)

## graph of predicted x length
ggplot(data = pva_m, aes(pva_m, length)) +
  geom_point() +
  geom_point(data = lingcoddiets_male, aes(x = Ages, y = TL.cm), colour = 'blue')

#FEMALE

  ##### The females length ranges from 33.3 through 115.0, so I will create a dataframe that
  ##### predicts from 33.3 through 105 (because the oldest recorded age 17 is predicted to be 105
  ##### and the predictions for anything above 105 may not be biologically sound)

LengthPredicted_female <- data.frame(TL.cm = 89.5)
predict(fittedgrowth_francis_female, LengthPredicted_female)

ests_f <- bootfrancis_female$coefboot
L1_f <- 48    #ests_f[,"L1"]
L2_f <- 92   #ests_f[,"L2"]
L3_f <- 105   #ests_f[,"L3"]
r_f <- (L3_f-L2_f)/(L2_f-L1_f)

female_lengthrange <- seq(from = 33.3, to = 105, by = 0.1)
pva_f <- numeric(length(female_lengthrange))
for (i in 1:length(female_lengthrange)) {
  pva_f[i] <- ages_francis_female[1] + (ages_francis_female[2]-ages_francis_female[1])*((log(1-(1-r_f^2)*((female_lengthrange[i]-L1_f)/(L3_f-L1_f))))/(2*log(r_f)))
}

## convert to dataframe and add a column with associated lengths
pva_f <- as.data.frame(pva_f)
pva_f <- mutate(pva_f, length = seq(from = 33.3, to = 105, by = 0.1))

## Since we only predicted until 105, but there is length's up to 115, I will create a new dataframe
## to fill in the rest of the data with 105 - 115 all "predicting" age 17
seq_105 <- seq(from = 105.1, to = 115, by = 0.1)
pva_f_105up <- cbind(rep(17, times = 100), seq_105)
colnames(pva_f_105up) <- c("pva_f", "length")

###
### This is the final dataframe with predicted vals for length
###
pva_f <- rbind(pva_f, pva_f_105up)

## graph of predicted x length
ggplot(data = pva_f, aes(pva_f, length)) +
  geom_point() +
  geom_point(data = lingcoddiets_female, aes(x = Ages, y = TL.cm), colour = 'blue')



### Okay, now lets add these predicted ages back into a complete dataset with 
### males and females together


## MALES
# Start by looking at just males from the working mastersheet
lingcodworking_male <- lingcod_rockfish %>% 
  filter(Sex.1 == "M")

# rename columns in pva_m to match lingcodworking_male for merging
pva_m <- pva_m %>% 
  rename(
    age_p = pva_m,
    TL.cm = length
  )

# merge the datasets
lingcodworking_male <- left_join(x = lingcodworking_male, y = pva_m, by = "TL.cm")
      ## there are 21 NAs because they have two numbers following the decimal


## FEMALES
# Start by looking at just females from the working mastersheet
lingcodworking_female <- lingcod_rockfish %>% 
  filter(Sex.1 == "F")

# rename columns in pva_m to match lingcodworking_male for merging
pva_f <- pva_f %>% 
  rename(
    age_p = pva_f,
    TL.cm = length)

# merge the datasets
lingcodworking_female <- left_join(x = lingcodworking_female, y = pva_f, by = "TL.cm")

    ## Produced random NA's so lets try to fix those
lingcodworking_female_NAs <- subset(lingcodworking_female,is.na(age_p))
lingcodworking_female_NAs$age_p <- NULL
lingcodworking_female_NAs$age_p <- pva_f$age_p[match(lingcodworking_female_NAs$TL.cm,pva_f$TL.cm)]
    ## Still producing NA's when it shouldn't....


## Let's go ahead and just merge the male and female datasets together to create our working dataset


#### ALTERNATIVELY, I can create a function with an if statement

predict_age_m <- function(length) {
if(length >= 83.8) {
  pred_age <-  14
} else {
  pred_age <- ages_francis_male[1] + (ages_francis_male[2]-ages_francis_male[1])*((log(1-(1-r_m^2)*((length-L1_m)/(L3_m-L1_m))))/(2*log(r_m)))
}
return(pred_age) 
}

predict_age_f <- function(length) {
  if(length >= 105) {
    pred_age <-  17
  } else {
    pred_age <- ages_francis_female[1] + (ages_francis_female[2]-ages_francis_female[1])*((log(1-(1-r_f^2)*((length-L1_f)/(L3_f-L1_f))))/(2*log(r_f)))
  }
  return(pred_age) 
}

#this helps prevent the warning messages.
v_predict_age_m <- Vectorize(predict_age_m)
v_predict_age_f <- Vectorize(predict_age_f)

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



### Now let's remove the decimal places for the predicted age
lingcod_male$age_pred <- floor(lingcod_male$age_pred)
lingcod_female$age_pred <- floor(lingcod_female$age_pred)

### Replace NA ages with the predicted ages
lingcod_male$Ages[is.na(lingcod_male$Ages)] <- lingcod_male$age_pred[is.na(lingcod_male$Ages)]
lingcod_female$Ages[is.na(lingcod_female$Ages)] <- lingcod_female$age_pred[is.na(lingcod_female$Ages)]



### Now let's plot what the weight of sebastes in each lingcod age class

plot(lingcod_male$Ages, lingcod_male$wt..Sebastes.Total, 
     xlim = c(0 , 17), ylim = c(0 , 600),
     xlab = "Ages", ylab = "Sebastes consumed in weight", main = "Male")
plot(lingcod_female$Ages, lingcod_female$wt..Sebastes.Total, 
     xlim = c(0 , 17), ylim = c(0 , 600),
     xlab = "Ages", ylab = "Sebastes consumed in weight", main = "Females")
    # This cuts off one outlier

plot(lingcod_female$Ages, lingcod_female$wt..Sebastes.Total)













