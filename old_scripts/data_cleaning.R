library(tidyverse)

# lingcod_working dataset is the dataset with total weight and number of all prey items
#                         and all columns of diet info retained
# lingcod_rockfish dataset is the dataset with a subset of relevant location and fish info
#                          and only columns including rockfish prey info
  

### Read the data in
lingcod <- read.csv("data/Lingcod with Gut Contents_All_BB.csv")
colnames(lingcod)

## Calculate the sum of all the gut content by weight
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

# New table with total food in stomach  
lingcod_working <- left_join(lingcod, lingcod_gut_content_wt, by = "Sample.ID") 
lingcod_working <- left_join(lingcod_working, lingcod_gut_content_X, by = "Sample.ID")
lingcod_working <- lingcod_working %>% # Remove the NA rows at the end
  slice(1:1321)


### There are two columns for sex, Sex and Sex.1. They don't match.
### What is the difference between the two?
lingcod_sexColumns <- lingcod %>% 
  select(Sample.ID, Sex, Sex.1)

lingcod_sexColumns %>% 
  group_by(Sex) %>%
  tally()

lingcod_sexColumns %>% 
  group_by(Sex.1) %>%
  tally()

### Create working dataset with relevant columns and make a
### new column grouping all Sebastes species into one
lingcod_rockfish <- lingcod_working %>%
  select(Sample.ID, Ages, Date..yymmdd., Port, Location, Lat, Long,
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


### Let's take a look at the Unknown Sex values and see if they're even useful
lingcod_rockfish %>% 
  group_by(Sex.1) %>% 
  tally()
  ## So, there are 55 blanks, 1 "O", 630 "F", and 960 "M"

Sex_0 <- lingcod_rockfish %>% 
  filter(Sex.1 == 0)  
## There is 1 observation, with only SampleID and Sebastes wt provided

Sex_NA <- lingcod_rockfish %>% 
  filter(Sex.1 == "")
## All of these observations are just a bunch of 0s

### Replace all 0s in Ages column with NAs
lingcod_rockfish$Ages <- replace(lingcod_rockfish$Ages, lingcod_rockfish$Ages == 0, NA)


### It is important to determine the ages that are NA because we are
### determining Lingcod diet by age

### How did Bonnie et al Age the fish? Otoliths?

### We could possibly determine an Age based off of the length of the fish
### So, let's look at the total length range for each known Age
LengthRangeByAge <- lingcod_rockfish %>% 
  group_by(Ages) %>%
  summarise(
    MaxLengthByAge = max(TL.cm, na.rm = T),
    MinLengthByAge = min(TL.cm, na.rm = T),
    MeanLengthByAge = mean(TL.cm, na.rm = T),
    n = n()
  ) %>%
  arrange(Ages)

189/sum(LengthRangeByAge$n) #14% of data have missing Age values


### Plotting Ages x Total Weight of Sebastes found in stomach content
### with Males/Females as different colors
ggplot(lingcod_rockfish, aes(x = Ages, y = wt..Sebastes.Total, color = Sex.1)) +
  geom_point() +
  theme_classic()


### Plotting Ages x Total Length with males/females as different colors
ggplot(lingcod_rockfish, aes(x = Ages, y = TL.cm, color = Sex.1)) +
  geom_point() +
  theme_classic()
  

  
### Lets look at how the weight distribution of Sebastes found in the 
### stomach varies with age

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


### Look at how Age varies with Depth
ggplot(lingcod_rockfish, aes(x = Ages, y = Depth.ft, color = Sex.1)) +
  geom_point() +
  theme_classic()

### Look at how body weight in kg varies with age between sex
ggplot(lingcod_rockfish, aes(x = Ages, y = Wt.kg, color = Sex.1)) +
  geom_point() +
  theme_classic()

### Lets look at how many observations there are that actually 
### contain sebastes in the stomach
lingcodWithSebastes <- lingcod_rockfish %>%
  filter(wt..Sebastes.Total>0) %>% 
  tally() # there are 165 observations with sebastes in stomach


### How many samples do we have for males and females for each age class
FreqOfSexByAge <- lingcod_rockfish %>% 
  group_by(Ages, Sex.1) %>%
  summarise(
    n = n()
  ) %>%
  arrange(Ages)


  