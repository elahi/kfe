library(tidyverse)
dat <- read_csv("data/seastar_counts - Sheet1.csv")

dat
datL <- dat %>%
  select(DATE, YEAR, SITE, QUADSIZE, rock_percent:pisaster_gig) %>% 
  gather(key = SPECIES, value = COUNT, patiria_min:pisaster_gig) %>% 
  mutate(COUNT10 = COUNT)

datL

write.csv(datL, "data_output/seastar_counts_jim.csv")
