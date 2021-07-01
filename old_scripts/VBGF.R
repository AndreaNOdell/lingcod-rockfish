library(tidyverse)
library(FSA)
library(FSAdata)
library(nlstools)
library(car)

###
### Subsetting data into Age and Length datasets for each sex separately
###

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
sv_schnutes_males <- vbStarts(TL.cm~Ages,data=lingcoddiets_male,type="Schnute", ages2use=c(2,12))
unlist(sv_males)
      ### L1 = 27, L3 = 81.9, K = 0.3291414

sv_francis_males <- vbStarts(TL.cm~Ages,data=lingcoddiets_male,type="Francis",tFrancis=c(2,12))
unlist(sv_francis_males)
      ### L1 = 27, L2 = 67.79846, L3 = 81.90000

sv_gallucci_males <- vbStarts(TL.cm~Ages,data=lingcoddiets_male,type="GallucciQuinn")
unlist(sv_francis_males)

sv_ogle_males <- vbStarts(TL.cm~Ages,data=lingcoddiets_male,type="Ogle", valOgle=c(tr=1))
unlist(sv_ogle_males)
      ### Linf = 80.8869502, K = 0.3291414, Lr = 42.2222342

# Females

sv_schnutes_females <- vbStarts(TL.cm~Ages,data=lingcoddiets_female,type="Schnute", ages2use=c(2,17))
unlist(sv_females)
      ## L1 = 43.9428571, L3 = 112, K = 0.1452901

sv_francis_females <- vbStarts(TL.cm~Ages,data=lingcoddiets_female,type="Francis",tFrancis=c(2,17))
unlist(sv_francis_females)
      ### L1 = 43.94286, L2 = 90.91599, L3 = 112.00000


###
### Call in the Growth function
###

vbschnute <- vbFuns("Schnute")

vbfrancis <- vbFuns("Francis")
ages_francis_male <- c(2,12)
ages_francis_female <- c(2,17)

vbgalluci <- vbFuns("GallucciQuinn")

vbFuns("Typical")

###
### Fit growth model to the data
###

# male

fittedgrowth_male <- nls(TL.cm~vbschnute(Ages, L1,L3,K,t1=2,t3=12),
                        data=lingcoddiets_male,
                        start=sv_schnutes_males)
summary(fittedgrowth_male,correlation=TRUE)
boot_male <- nlsBoot(fittedgrowth_male,niter=200)
confint(boot_male,plot=TRUE)
plot(boot_male)
     # Parameters seem to still be incredibly correlated
    
    # Let's try using the francis parameterization
fittedgrowth_francis_male <- nls(TL.cm~vbfrancis(Ages,L1,L2,L3,t1=ages_francis_male[1],t3=ages_francis_male[2]), 
                                 data=lingcoddiets_male,start=sv_francis_males)
overview(fittedgrowth_francis_male)
bootfrancis_male <- nlsBoot(fittedgrowth_francis_male,niter=1000)
confint(bootfrancis_male,plot=TRUE)
plot(bootfrancis_male)
      # a little better but still somewhat correlated

    # Let's try Ogle to estimate age based off length
#tr = 1
#fittedgrowth_ogle_male <- nls(TL.cm~vbogle(Ages,Linf,K,tr,Lr),data=lingcoddiets_male,start=sv_ogle_males)
#bootogle_male <- nlsBoot(fittedgrowth_ogle_male)
#plot(bootogle_male)
  #Compute tr estimates
#tr_ogle_m <- coef(fittedgrowth_ogle_male)[["tr"]]
    # Set critical length (the length we want to predict)
#Lr <- 62.3
# Set ages and confidence interval at which to make some predictions
#ages_ogle_m <- c(0,14)
#cl_ogle <- 0.95
#bootN2 <- bootogle_male
#bootN2$coefboot <- cbind(bootN2$coefboot,Lr)
#ests_ogle_m=coef(fittedgrowth_ogle_male,Lr)
#confint(bootogle_male,conf.level=cl_ogle)
#vbogle(ages_ogle_m,c(coef(fittedgrowth_ogle_male),Lr))
#predict(bootN2,vbogle,t=ages_ogle_m,conf.level=cl_ogle)[,3:4])



# female

fittedgrowth_female <- nls(TL.cm~vbschnute(Ages, L1,L3,K,t1=2,t3=17),
                         data=lingcoddiets_female,
                         start=sv_schnutes_females)
summary(fittedgrowth_female,correlation=TRUE)
boot_female <- nlsBoot(fittedgrowth_female,niter=1000)
confint(boot_female,plot=TRUE)
plot(boot_female)
    # Parameters seem to still be incredibly correlated

# Let's try using the francis parameterization
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


###
### Predicting age using model
###

# male

AgePredicted_male <- data.frame(TL.cm = 62.3)
predict(fittedgrowth_francis_male, AgePredicted_male)
ests_m <- bootfrancis_male$coefboot
L1_m <- ests_m[,"L1"]
L2_m <- ests_m[,"L2"]
L3_m <- ests_m[,"L3"]
r_m <- (L3_m-L2_m)/(L2_m-L1_m)

confidenceint_male <- matrix(data = NA, nrow=17, ncol = 2)

for(i in AgePredicted$Ages) {
holder <- L1_m + (L3_m-L1_m)*(1-r_m^(2*(AgePredicted$Ages[i]-ages_francis_male[1])/(ages_francis_male[2]-ages_francis_male[1])))/(1-r_m^2)
confidenceint_male[i,] <- quantile(holder,c(0.025,0.975))
}

        # Predicting  using alternative code
male_agerange <- 0:14
LCI_male <- UCI_male <- numeric(length(male_agerange))
for (i in 1:length(male_agerange)) {
  pv_m <- L1_m + (L3_m-L1_m)*(1-r_m^(2*(male_agerange[i]-ages_francis_male[1])/(ages_francis_male[2]-ages_francis_male[1])))/(1-r_m^2)
  LCI_male[i] <- quantile(pv_m,0.025)
  UCI_male[i] <- quantile(pv_m,0.975)
}

        ## Create a dataframe with predicted values and upper/lower confidence intervals
agepredictrange_m <- data.frame(Ages = 0:14)
predict_m <- predict(fittedgrowth_francis_male, agepredictrange_m)
predictions_confint_m <- rbind(LCI_male, predict_m, UCI_male)
colnames(predictions_confint_m) <- c(0:14)
predictions_confint_m<-as.data.frame(t(predictions_confint_m))


#female

predict(fittedgrowth_francis_female, AgePredicted)
ests <- bootfrancis_female$coefboot

confidenceint_female <- matrix(data = NA, nrow=17, ncol = 2)

for(i in AgePredicted$Ages) {
  holder <- L1 + (L3-L1)*(1-r^(2*(AgePredicted$Ages[i]-ages_francis_female[1])/(ages_francis_female[2]-ages_francis_female[1])))/(1-r^2)
  confidenceint_female[i,] <- quantile(holder,c(0.025,0.975))
}


      # Predicting  using alternative code

AgePredicted_female <- data.frame(Ages = 12)
predict(fittedgrowth_francis_female, AgePredicted_female)
ests_f <- bootfrancis_female$coefboot
L1_f <- ests_f[,"L1"]
L2_f <- ests_f[,"L2"]
L3_f <- ests_f[,"L3"]
r_f <- (L3_f-L2_f)/(L2_f-L1_f)

female_agerange <- 0:14
LCI_female <- UCI_female <- numeric(length(female_agerange))
for (i in 1:length(female_agerange)) {
  pv_f <- L1_f+ (L3_f-L1_f)*(1-r_f^(2*(female_agerange[i]-ages_francis_female[1])/(ages_francis_female[2]-ages_francis_female[1])))/(1-r_f^2)
  LCI_female[i] <- quantile(pv_f,0.025)
  UCI_female[i] <- quantile(pv_f,0.975)
}


      ### Welp, the equation is predicting length, not age. So let's do some algebra 
      ### and have it predict Age
LengthPredicted_female <- data.frame(TL.cm = 89.5)
predict(fittedgrowth_francis_female, LengthPredicted_female)

female_lengthrange <- 33:115
lengthLCI_female <- lengthUCI_female <- numeric(length(female_lengthrange))
for (i in 1:length(female_lengthrange)) {
pva_f <- ages_francis_female[1] + (ages_francis_female[2]-ages_francis_female[1])*((log(1-(1-r_f^2)*((female_lengthrange[i]-L1_f)/(L3_f-L1_f))))/(2*log(r)))
}

      ## Create a dataframe with predicted values and upper/lower confidence intervals
agepredictrange_f <- data.frame(Ages = 0:14)
predict_f <- predict(fittedgrowth_francis_female, agepredictrange_f)
predictions_confint_f <- rbind(LCI_female, predict_f, UCI_female)
colnames(predictions_confint_f) <- c(0:14)
predictions_confint_f<-as.data.frame(t(predictions_confint_f))



###
### Plot predictions with confidence intervals
###

# Male
fitPlot(fittedgrowth_francis_male,xlab="Age",ylab="Total Length (mm)",main="Males")
lines(UCI_male~male_agerange,type="l",col="blue",lwd=2,lty=2) 
lines(LCI_male~male_agerange,type="l",col="blue",lwd=2,lty=2)

# Female
fitPlot(fittedgrowth_francis_female,xlab="Age",ylab="Total Length (mm)",main="Females")
lines(UCI_female~female_agerange,type="l",col="blue",lwd=2,lty=2) 
lines(LCI_female~female_agerange,type="l",col="blue",lwd=2,lty=2)


###
### Okay, let's try to fill in the missing values
###

# Males
  # first, extract rows with NA ages for males
NAage_m <- lingcod_rockfish %>%
  filter(Sex.1 == "M" & is.na(Ages)) %>% 
  select(Sample.ID, Ages, TL.cm)
  

#Females

