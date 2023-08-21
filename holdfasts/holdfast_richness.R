################################################################################
##' @title Explore species accumulation curves with holdfast data
##' @author Robin Elahi
##' @date 2023-08-18
##' @log 
################################################################################

##### PACKAGES, DATA #####
library(here) 
library(tidyverse) 
library(lubridate)
library(viridis)
library(vegan)
library(iNEXT)

# Theme
theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank()))

#d <- read_csv(here("data", "andrews_henry_combined.csv"))
d <- read_csv(here("data", "andrews_1945_230623.xlsx - andrews_henry_combined.csv"))
names(d)

d <- d %>% 
  mutate(date = lubridate::ymd(paste(year, month, day, sep = "-")))

# Remove rows that were not sampled (ie NAs for species)
d <- d %>% 
  filter(!is.na(Balanophyllia_elegans))

# Order sites north to south
d <- d %>% 
  mutate(site = as.factor(site), 
         site = fct_relevel(site, 
                            "Point_Pinos", 
                            "Otter_Cove", 
                            "Lovers_Point", 
                            "Cabrillo_Point", 
                            "Breakwater", 
                            "Pescadero_Point",
                            "Stillwater_Cove"))

# Carmel vs MBay
d <- d %>% 
  mutate(bay = case_when(site == "Pescadero_Point" ~ "Carmel", 
                         site == "Stillwater_Cove" ~ "Carmel", 
                         TRUE ~ "Monterey"))

##### SPECIES ACCUMULATION CURVES #####

# Create metadata
d_meta <- d %>% 
  select(year:site, date, bay) %>%
  mutate(era = ifelse(year > 2018, "2019", "1930s"))
d_meta

# Create matrix
m <- d %>% 
  select(Balanophyllia_elegans:Gobiesox_eugrammus) %>% 
  as.data.frame()

d_counts <- rowSums(m)
d_counts

d_meta <- d_meta %>% 
  mutate(n_ind = d_counts, 
         n_spp = specnumber(m))

d_meta %>% 
  ggplot(aes(n_ind, n_spp, color = site, shape = era)) + 
  labs(x = "No. of individuals", y = "No. of observed species") + 
  geom_point(size = 3, alpha = 1) + 
  scale_color_viridis_d() 

ggsave(here("figs", "holdfast_richness_abundance.pdf"), height = 5, width = 7)

#### IGNORE BELOW ####

#### VEGAN ####

data("BCI")
BCI
dim(BCI)
accurve<-specaccum(BCI, method="random", permutations=100)
plot(accurve$sites, accurve$richness,
     xlab="Number of Sites",
     ylab="Species Richness",
     main="Barro Colorado Island")

accurve <- specaccum(m, method="random", permutations=100)
plot(accurve$sites, accurve$richness,
     xlab="Number of Sites",
     ylab="Species Richness",
     main="Holdfasts")


#### VEGAN 2 #### 

#total number of species at each site (row of data)
S <- specnumber(BCI)

# Number of INDIVIDULS per site (?)
raremax <- min(rowSums(BCI)) # = 340; 


# rarefy, w/ raremax as input (?)
Srare <- rarefy(BCI, raremax)


#Plot rarefaction results
par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(BCI, raremax))")
abline(0, 1)
rarecurve(BCI, step = 20, 
          sample = raremax, 
          col = "blue", 
          cex = 0.6,
          main = "rarecurve()")

#### VEGAN 3 ####

# https://rpubs.com/brouwern/iNEXTvVEGAN

#total number of species at each site (row of data)
S <- specnumber(m)
S

# Number of INDIVIDULS per site (?)
rowSums(m)
raremax <- min(rowSums(m)) 
raremax

# rarefy, w/ raremax as input (?)
Srare <- rarefy(m, raremax)
Srare

#Plot rarefaction results
par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(m, raremax))")
abline(0, 1)
rarecurve(m, step = 20, 
          sample = raremax, 
          col = "blue", 
          cex = 0.6,
          main = "rarecurve()")

rare_tidy <- rarecurve(m, step = 20, 
          sample = raremax, 
          tidy = TRUE)

#### iNEXT ####

m_transpose <- t(m)
m_transpose

out <- iNEXT(m_transpose, q = 0, datatype = "abundance", size = NULL, endpoint = NULL, 
      knots = 40, se = TRUE, conf = 0.95, nboot = 50)

out
ggiNEXT(out, type = 1, facet.var = "Assemblage")
str(out)

data(spider)
str(spider)


