library(tidyverse)
library(gridExtra)
library(mgcv)


# Load in ONE dataset to use
# df <- miceadds::load.Rdata2(filename="cleaned_data/lingcod_rockfish_CA.Rdata")

# We'll go ahead and just use the dataset with all of the ports included.
df <- miceadds::load.Rdata2(filename="cleaned_data/lingcod_rockfish.Rdata")

# To calculate the fraction of lingcod diet that is rockfish, I simply divided 
# the total Sebastes weight/number by the total weight/number of all gut contents 
# for each lingcod individual.

# Dataframe with the average diet fraction by weight that is sebastes for each age and sex 
dietfrac_by_age_wt <- df %>% 
  group_by(Sex, age_useful) %>% 
  summarise(mean = mean(gut.ratio.sebastes.wt), sd = sd(gut.ratio.sebastes.wt), 
            n = n(), max = max(gut.ratio.sebastes.wt), min = min(gut.ratio.sebastes.wt))
dietfrac_by_age_wt$sd[is.na(dietfrac_by_age_wt$sd)] <- 0

# Dataframe with the average diet fraction by number that is sebastes for each age and sex 
#dietfrac_by_age_X <- df %>% 
#  group_by(Sex.1, age_useful) %>% 
#  summarise(mean = mean(gut.ratio.sebastes.X), sd = sd(gut.ratio.sebastes.X), n = n(),
#            max = max(gut.ratio.sebastes.X), min = min(gut.ratio.sebastes.X))
#dietfrac_by_age_X$sd[is.na(dietfrac_by_age_X$sd)] <- 0

# Set up bin size
bin_size = 8
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
  select(Sex, Bin, gut_sum_by_sex.lengthbin)
# now add back to dataset
df = left_join(df, gut_summary_sum, by = c("Sex", "Bin"))

df$gut.weighting = df$wt..Total/df$gut_sum_by_sex.lengthbin

gut_summary = df %>% 
  group_by(Sex, Bin) %>% 
  summarise(mean = mean(gut.ratio.sebastes.wt), 
            sd = sd(gut.ratio.sebastes.wt),
            weightedmean = weighted.mean(gut.ratio.sebastes.wt, gut.weighting),
            gutsumcheck = sum(gut.weighting),
            n = n())
gut_summary[is.na(gut_summary)] = 0



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


grid.arrange(female_dietcomp, male_dietcomp, ncol=2)

gut_summary$Bin = as.numeric(gut_summary$Bin)
 
# pdf('plots/male_dietcomp_CAonly.pdf')
# male_dietcomp
# dev.off()

diet_f = as.data.frame(gut_summary) %>% 
  filter(Sex == "F") 
diet_model_f = gam(mean ~ s(Bin), data = diet_f, method = "REML")

diet_m = as.data.frame(gut_summary) %>% 
  filter(Sex == "M") %>% 
  filter(mean < .24) # removed the outlier
diet_model_m = gam(mean ~ s(Bin), data = diet_m, method = "ML")

plot(diet_model_f,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)
plot(diet_model_m,pages=1, residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)










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







