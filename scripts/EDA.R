library(tidyverse)
load("cleaned_data/lingcod_rockfish.Rdata")
load("cleaned_data/lingcod_rockfish_CA.Rdata")

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
