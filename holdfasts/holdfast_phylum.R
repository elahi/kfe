################################################################################
##' @title Explore holdfast abundance data by phylum
##' @author Robin Elahi
##' @date 2023-08-21
##' @log 
################################################################################

##### PACKAGES, DATA #####
library(here) 
library(tidyverse) 
library(lubridate)
library(viridis)
library(vegan)

# Theme 
theme_set(theme_bw(base_size = 12) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank()))

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
                         TRUE ~ "Monterey")) %>% 
  mutate(era = ifelse(year > 2018, "2019", "1930s"))

# Load taxon table
d_taxa <- read_csv(here("data", "andrews_1945_230623.xlsx - taxon_table_re.csv"))
d_taxa

# Pivot to long
d_long <- d %>% 
  pivot_longer(., cols = Balanophyllia_elegans:Unidentified, 
               names_to = "species", values_to = "abundance")

d_long

# Check that species names match
unique(d_long$species)  
unique(d_long$species) %in% d_taxa$species

# Join
d_long <- d_taxa %>% 
  select(species, phylum) %>% 
  left_join(d_long, ., by = "species")

##### SUMMARIZE ABUNDANCE #####
d_abund <- d_long %>% 
  group_by(era, site, date, phylum) %>% 
  summarize(abundance = sum(abundance))

d_abund

abund_summary <- d_abund %>% 
  group_by(era, site, phylum) %>% 
  summarize(mean = mean(abundance), 
            sd = sd(abundance))

abund_summary %>% 
  ggplot(aes(site, mean, color = site, shape = era)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray") +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), lwd = 1) +
  scale_color_viridis_d() +
  facet_wrap(~ phylum, scales = "free_y") + 
  theme(legend.position = "none") + 
  labs(y = "Abundance (mean, SD)", x = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(here("figs", "holdfast_phylum_abundance.pdf"), height = 7, width = 7)


##### SUMMARIZE RICHNESS #####
d_abund <- d_long %>% 
  group_by(era, site, date, phylum) %>% 
  summarize(abundance = sum(abundance))

d_abund

abund_summary <- d_abund %>% 
  group_by(era, site, phylum) %>% 
  summarize(mean = mean(abundance), 
            sd = sd(abundance))

abund_summary %>% 
  ggplot(aes(site, mean, color = site, shape = era)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray") +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), lwd = 1) +
  scale_color_viridis_d() +
  facet_wrap(~ phylum, scales = "free_y") + 
  theme(legend.position = "none") + 
  labs(y = "Abundance (mean, SD)", x = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(here("figs", "holdfast_phylum_abundance.pdf"), height = 7, width = 7)

