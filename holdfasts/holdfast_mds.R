################################################################################
##' @title Explore community mds with holdfast data
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

# Theme for mds plots
theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid = element_blank(), 
                  axis.text.x = element_blank(), 
                  axis.text.y = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.title = element_blank(), 
                  strip.background = element_blank()))

d <- read_csv(here("data", "andrews_henry_combined.csv"))

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
                            "Otter_Point", 
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

##### MDS #####

# Create metadata
d_meta <- d %>% 
  select(year:site, date, bay) %>%
  mutate(era = ifelse(year > 2018, "2019", "1930s"))
d_meta

# Create matrix for mds
m <- d %>% 
  select(Balanophyllia_elegans:Gobiesox_eugrammus) %>% 
  as.data.frame()

# Relativize to total
m_rel <- decostand(m, method = "total")
rowSums(m_rel)

# Square-root transformation
m_sqrt <- m^(1/2)

# Create distance matrix
m_dist <- vegdist(m_rel, method = "bray") %>% as.matrix(labels = T)

# Run mds
m_mds <- metaMDS(comm = m_sqrt, distance = "bray", trace = FALSE, autotransform = FALSE)

m_mds$stress

# Join MDS results to metadata
mds_xy <- data.frame(m_mds$points) 
d_meta <- cbind(d_meta, mds_xy)

# Find the convex hull of the points being plotted
hull <- d_meta %>%
  group_by(site) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()
hull

d_meta %>% 
  ggplot(aes(MDS1, MDS2, color = site, shape = era)) + 
  geom_polygon(data = hull, aes(x = MDS1, y = MDS2, color = site),
               inherit.aes = FALSE, alpha = 0.1) +
  scale_color_viridis_d() + 
  geom_point(size = 3, alpha = 0.5) +
  labs(caption = "Transformation: sqrt(abundance)")

ggsave(here("figs", "holdfast_mds_abundance.pdf"), height = 5, width = 7)
