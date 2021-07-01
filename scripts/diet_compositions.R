library(tidyverse)
load(file = "cleaned_data/lingcod_rockfish_CA.Rdata")
load(file = "cleaned_data/lingcod_rockfish.Rdata")

# To calculate the fraction of lingcod diet that is rockfish, I simply divided 
# the total Sebastes weight/number by the total weight/number of all gut contents 
# for each lingcod individual.

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

bins <- seq((min(lingcod_rockfish_CA$TL.cm)-2), (max(lingcod_rockfish_CA$TL.cm)+3), by=5)
TL.cm <- as.numeric(lingcod_rockfish_CA$TL.cm)
labels <- bins
rangelabels <- paste(head(labels,-1), tail(labels,-1), sep="-")
lingcod_rockfish_CA$Bin <- cut(TL.cm, bins, rangelabels)

ggplot(lingcod_rockfish_CA, aes(x = Bin, y = gut.ratio.sebastes.wt, color = Sex.1)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "length bins (cm)", y = "Diet fraction that is Sebastes")


bins_full <- seq((min(lingcod_rockfish$TL.cm)-2), (max(lingcod_rockfish$TL.cm)), by=5)
TL.cm_full <- as.numeric(lingcod_rockfish$TL.cm)
labels_full <- bins_full
rangelabels_full <- paste(head(labels_full,-1), tail(labels_full,-1), sep="-")
lingcod_rockfish$Bin <- cut(TL.cm_full, bins_full, rangelabels_full)

ggplot(lingcod_rockfish, aes(x = Bin, y = gut.ratio.sebastes.wt, color = Sex.1)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "length bins (cm)", y = "Diet fraction that is Sebastes")







