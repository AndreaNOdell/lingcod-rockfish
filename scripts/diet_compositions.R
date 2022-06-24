library(tidyverse)
library(gridExtra)
library(mgcv)
require(plotrix)
require(stats4)
require(optimx)
require(numDeriv)
require(MASS)
require(fitdistrplus)
library(FSA)
library(gamlss)

# Load in ONE dataset to use
# df <- miceadds::load.Rdata2(filename="cleaned_data/lingcod_rockfish_CA.Rdata")

# We'll go ahead and just use the dataset with all of the ports included.
df <- miceadds::load.Rdata2(filename="cleaned_data/lingcod_rockfish_w_gutcontent.Rdata")

# To calculate the fraction of lingcod diet that is rockfish, I simply divided 
# the total Sebastes weight/number by the total weight/number of all gut contents 
# for each lingcod individual.

df_F = df %>% 
  filter(Sex == "F") %>% 
  filter(Port %in%  c("BRO", "LOS", "SLO", "HMB", "MON", "EME", "MOR", "BOD", "FTB", "EUR", "SBA", 
                      "SDG"))
    # 233 observations

df_M = df %>% 
  filter(Sex == "M") %>% 
  filter(Port %in%  c("BRO", "LOS", "SLO", "HMB", "MON", "EME", "MOR", "BOD", "FTB", "EUR", "SBA", 
                      "SDG"))
    # 409 observations

#### GAMlss approach ####
m1 = gamlss(gut.ratio.sebastes.wt ~ cs(TL.cm), family = BEINF(), data = df_F, method = RS())
m2 = gamlss(gut.ratio.sebastes.wt ~ ps(TL.cm), family = BEINF(), data = df_F, method = RS())
m3 = gamlss(gut.ratio.sebastes.wt ~ fp(TL.cm), family = BEINF(), data = df_F, method = RS())
AIC(m1,m2,m3)

plot(gut.ratio.sebastes.wt ~ TL.cm, data = df_F, xlab = "diet fraction", ylab = "total length", main = "Female")
lines(df_F$TL.cm[order(df_F$TL.cm)], fitted(m1)[order(df_F$TL.cm)], col = "red")
lines(df_F$TL.cm[order(df_F$TL.cm)], fitted(m2)[order(df_F$TL.cm)], col = "blue")
lines(df_F$TL.cm[order(df_F$TL.cm)], fitted(m3)[order(df_F$TL.cm)], col = "green")
legend("left", legend=c("cs", "ps", "fp"), col=c("red", "blue", "green"), lty = 1, box.lty=0, cex = 0.7)

m4 = gamlss(gut.ratio.sebastes.wt ~ cs(TL.cm), family = BEINF(), data = df_M, method = RS())
m5 = gamlss(gut.ratio.sebastes.wt ~ ps(TL.cm), family = BEINF(), data = df_M, method = RS())
m6 = gamlss(gut.ratio.sebastes.wt ~ fp(TL.cm), family = BEINF(), data = df_M, method = RS())

plot(gut.ratio.sebastes.wt ~ TL.cm, data = df_M, xlab = "diet fraction", ylab = "total length", main = "Male")
lines(df_M$TL.cm[order(df_M$TL.cm)], fitted(m4)[order(df_M$TL.cm)], col = "red")
lines(df_M$TL.cm[order(df_M$TL.cm)], fitted(m5)[order(df_M$TL.cm)], col = "blue")
lines(df_M$TL.cm[order(df_M$TL.cm)], fitted(m6)[order(df_M$TL.cm)], col = "green")
legend("left", legend=c("cs", "ps", "fp"), col=c("red", "blue", "green"), lty = 1, box.lty=0, cex = 0.7)
AIC(m4,m5,m6)

#### GAM approach ####

gam_F = gam(gut.ratio.sebastes.wt ~ s(TL.cm), data = df_F, method = "REML")
plot(gam_F,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)

gam_M = gam(gut.ratio.sebastes.wt ~ s(TL.cm), data = df_M, method = "REML")
plot(gam_M,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)

# Now let's predict diet fractions for the length at age for lingcod from pop model
load("cleaned_data/lingcod_parms.Rdata")
lingcod_female = as.data.frame(lingcod$length.at.age["female",])
colnames(lingcod_female) = "TL.cm"
diet_comp_f = predict.gam(gam_F, lingcod_female)
plot(diet_comp_f, ylim = c(0,0.5), main = "Female")

lingcod_male = as.data.frame(lingcod$length.at.age["male",])
colnames(lingcod_male) = "TL.cm"
diet_comp_m = predict.gam(gam_M, lingcod_male)
plot(diet_comp_m, ylim = c(0,0.5), main = "Male")


#### Let's try Pam Moriarty's approach ####

# Check covariance of total stomach content and diet fraction
#females
plot(df_F$gut.ratio.sebastes.wt, df_F$wt..Total)
covariance_f = lm(wt..Total ~ gut.ratio.sebastes.wt, data = df_F)
#males
plot(df_M$gut.ratio.sebastes.wt, df_M$wt..Total)
covariance_m = lm(wt..Total ~ gut.ratio.sebastes.wt, data = df_M)

F_diet_data_for_model = df_F %>% 
  filter(age_useful > 5) %>% 
  dplyr::select(gut.ratio.sebastes.wt, wt..Total)
M_diet_data_for_model = df_M %>% 
  filter(age_useful > 5) %>% 
  dplyr::select(gut.ratio.sebastes.wt, wt..Total)

source('PM_dietfraction_method/run.model.R')#prints the parameter estimates for all 8 mixture model parameters, 
source('PM_dietfraction_method/model.comparison.R')#estimates the prey contribution, c_i, using the mixture model, weighted mean and mean and calculates  error for each estimate

run.model(F_diet_data_for_model[,1],F_diet_data_for_model[,2]) # 0.215
model.comparison(F_diet_data_for_model[,1],F_diet_data_for_model[,2],mat=T)

run.model(M_diet_data_for_model[,1],M_diet_data_for_model[,2]) # 0.139
model.comparison(M_diet_data_for_model[,1],M_diet_data_for_model[,2],mat=T)

load("cleaned_data/lingcod_parms.Rdata")
load("cleaned_data/F_length_variation.Rdata")
load("cleaned_data/M_length_variation.Rdata")

F_length_variation[13:20] = F_length_variation[length(F_length_variation)]
names(F_length_variation) = 1:20

M_length_variation[12:20] = M_length_variation[length(M_length_variation)]
names(M_length_variation) = 1:20

female_length_intervals = as.matrix(bind_cols(lower = lingcod$length.at.age[1,]-F_length_variation, upper = lingcod$length.at.age[1,]+F_length_variation))
male_length_intervals = as.matrix(bind_cols(lower = lingcod$length.at.age[2,]-M_length_variation, upper = lingcod$length.at.age[2,]+M_length_variation))


length_intervals <- array(c(female_length_intervals , male_length_intervals ) , dim = c(20, 2, 2), dimnames = list(c(1:20) , c("lower", "upper"), c("F", "M")))
rm(female_length_intervals, male_length_intervals)


run.model.age = function(sex, age) {
trail = df %>% 
  dplyr:: filter(Sex == sex) %>% 
  dplyr:: filter(TL.cm > length_intervals[age,"lower",sex]) %>% 
  dplyr:: filter(TL.cm < length_intervals[age,"upper",sex]) %>% 
  dplyr:: select(gut.ratio.sebastes.wt, wt..Total)
run.model(trail[,1],trail[,2])
}

F_dietfracs = numeric(20)
for(n in c(3:20)) {
  model = run.model.age(sex = "F", age = n)
  F_dietfracs[n] = model[8,1]$c_i
}

M_dietfracs = numeric(20)
for(n in c(3:20)) {
  model = run.model.age(sex = "M", age = n)
  M_dietfracs[n] = model[8,1]$c_i
}



#### Calculating Average and Weighted Average ####

# Set up bin size
bin_size = 5
round_to = 10

# Create length bins
bins <- seq(plyr::round_any(min(df$TL.cm), round_to, f = floor), plyr::round_any(max(df$TL.cm), round_to, f = ceiling), by= bin_size)
TL.cm <- as.numeric(df$TL.cm)
rangelabels <- paste(head(bins,-1), tail(bins,-1), sep="-")
df$Bin <- cut(TL.cm, bins, rangelabels)

# Let's add a grouped sum by age and sex to do our weightings
gut_summary = df %>% 
  group_by(Sex, Bin) %>% 
  summarise(mean = mean(gut.ratio.sebastes.wt), 
            sd = sd(gut.ratio.sebastes.wt),
            gut_sum_by_sex.lengthbin = sum(wt..Total),
            n = n())
gut_summary$sd[is.na(gut_summary$sd)] <- 0
gut_summary_sum = gut_summary %>% 
  dplyr::select(Sex, Bin, gut_sum_by_sex.lengthbin)
# now add back to dataset
df = left_join(df, gut_summary_sum, by = c("Sex", "Bin"))

df$gut.weighting = df$wt..Total/df$gut_sum_by_sex.lengthbin

# Look at summarized info (mean, weighted mean, sd, n) andn check if weighting adds to 1
gut_summary = df %>% 
  group_by(Sex, Bin) %>% 
  summarise(mean = mean(gut.ratio.sebastes.wt), 
            sd = sd(gut.ratio.sebastes.wt),
            weightedmean = weighted.mean(gut.ratio.sebastes.wt, gut.weighting),
            gutsumcheck = sum(gut.weighting),
            n = n())
gut_summary[is.na(gut_summary)] = 0

# plot out diet fractions
female_dietcomp = as.data.frame(gut_summary) %>% 
  filter(Sex == "F") %>% 
  ggplot() +
   geom_point(aes(x = Bin, y = weightedmean, size = n)) +
    theme_classic() +
    labs(x = "length bin", y = " Weighted mean", title = "Avg. frac. of diet that is Sebastes in female Lingcod") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_point(aes(x = Bin, y = mean, color = "mean", size = 0.05)) +
    ylim(0,1)

male_dietcomp =  as.data.frame(gut_summary) %>% 
  filter(Sex == "M") %>% 
  ggplot() +
    geom_point(aes(x = Bin, y = weightedmean, size = n)) +
    theme_classic() +
    labs(x = "length bin", y = " Weighted mean", title = "Avg. frac. of diet that is Sebastes in male Lingcod") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_point(aes(x = Bin, y = mean, size = 0.05, color = "mean")) +
    ylim(0,1)


### Fitting a GAM ####
# to the weighted mean of diet composition across length bins.

# must convert bin from categorical to numerical since gams need numerical values.
gut_summary$Bin = gsub("\\-.*","\\", gut_summary$Bin)
gut_summary$Bin = as.numeric(gut_summary$Bin) + 0.5*bin_size # Set bin column to be the midpoint of each bin

# female
diet_f = as.data.frame(gut_summary) %>% 
  filter(Sex == "F") 
diet_model_f = gam(mean ~ s(Bin), data = diet_f, method = "REML")
plot(diet_model_f,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
summary(diet_model_f)
gam.check(diet_model_f)

#male
diet_m = as.data.frame(gut_summary) %>% 
  filter(Sex == "M")# %>% 
#  filter(mean < .24) # removed the outlier
diet_model_m = gam(mean ~ s(Bin), data = diet_m, method = "REML")
plot(diet_model_m,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
summary(diet_model_m)
gam.check(diet_model_m)

#### create a vector of diet composition for age ###
load(file = "cleaned_data/lingcod_parms.Rdata")
lingcod_female = as.data.frame(lingcod$length.at.age["female",])
colnames(lingcod_female) = "Bin"
diet_comp_f = predict.gam(diet_model_f, lingcod_female)

lingcod_male = as.data.frame(lingcod$length.at.age["male",])
colnames(lingcod_male) = "Bin"
diet_comp_m = predict.gam(diet_model_m, lingcod_male)

plot(diet_comp_f)
plot(diet_comp_m)









#### Other stuff that isn't super important... ####

# The proportion of stomachs that contained rockfish of all the stomachs
lingcod_rockfish_CA_F <- lingcod_rockfish_CA %>% # filter out only Females
  filter(Sex.1 == "F")
lingcod_rockfish_CA_M <- lingcod_rockfish_CA %>% # filter out only Males
  filter(Sex.1 == "M")
# Divide the number of stomachs with Sebastes by the total number of stomachs
  # Both Sex combined
perc.occurrence <- (sum(lingcod_rockfish_CA$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA))*100 
  # Females
perc.occurrence.F <- (sum(lingcod_rockfish_CA_F$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA_F))*100
  # Males
perc.occurrence.M <- (sum(lingcod_rockfish_CA_M$X..Sebastes.Total > 0)/nrow(lingcod_rockfish_CA_M))*100
# Tried doing percent abundance by number and weight but idk if it's correct...
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


# Now let's look at whether diet compositions are sex-based or size-based by comparing 
# male and females of the same size

ggplot(lingcod_rockfish_CA, aes(x = Bin, y = gut.ratio.sebastes.wt, color = Sex.1)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "length bins (cm)", y = "Diet fraction that is Sebastes")



ggplot(lingcod_rockfish, aes(x = Bin, y = gut.ratio.sebastes.wt, color = Sex.1)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "length bins (cm)", y = "Diet fraction that is Sebastes")







